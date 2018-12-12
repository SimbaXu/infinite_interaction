"""This script synthesize a 1-dof controller for 1 direction of the
3-dimensional Cartersian admittance control problem.
"""
import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import yaml
from scipy.interpolate import interp1d

# constant
Ts = 0.008  # sampling internval


def plant(gain_E=20, m_int=0.1):
    """ An old plant model. Not in use currently.
    """
    z = co.tf([1, 0], [1], Ts)
    s = co.tf([1, 0], [1])
    zero = co.tf([0], [1], Ts)
    one = co.tf([1], [1], Ts)
    R1 = co.c2d(s / (1 + 0.0437 * s), Ts)
    E = - co.c2d(gain_E / s, Ts)  # exo agent
    R2 = - m_int / Ts ** 2 * (1 - z ** (-1)) ** 2  # internally induced

    P = Ss.tf_blocks([[zero, z**(-1) * R1],  # velocity output
                      # to form T = L / (1 + L)
                      [zero, z**(-1) * (R1 * E + R2)],
                      # to form S = 1 / (1 + L)
                      [one, z**(-1) * (R1 * E + R2)],
                      [z**(-1), z**(-2) * (R1 * E + R2)]])
    return P


def plantMdelta(E_gain=20, m_int=0.1, wI=0.9995, sys_type='22_mult_unt',
                N_in=1, N_out=1, tau_R1=0.0437):
    """ Plant Model with multiplicative uncertainty.

    Diagram is drawn in page 4, notebook.

    Args:
        wI (float): Can be a complex function.  Uncertainty in the exogeneous agent.
        sys_type (str, optional): Type of plant.
            - 22_mult_unt: (2, 2) system with multiplicative uncertainty; See Muji-1811, page4a.
            - 22_inv_mult_unt: (2, 2) system with inverse multiplicative uncertainty;
            - 33_mult_unt: (3, 3) system with multiplicative uncertainty. See Muji-1811, page4.
        N_in (int, optional): Nb. of delay steps in the input path.
        N_out (int, optional): Nb. of delay steps in the output path.
    """
    # constant transfer function
    z = co.tf([1, 0], [1], Ts)
    s = co.tf([1, 0], [1])
    zero = co.tf([0], [1], Ts)
    one = co.tf([1], [1], Ts)

    E_gain = - abs(E_gain)

    # basic blocks
    R1 = co.c2d(s / (1 + tau_R1 * s), Ts)
    sumz = co.c2d(1 / s, Ts)  # discrete-time integrator
    E = E_gain * sumz  # exo agent
    R2 = - m_int / Ts ** 2 * (1 - z ** (-1)) ** 2  # internally induced

    if sys_type == "22_inv_mult_unt":
        P = Ss.tf_blocks([
            [0, 0, R1 * z**(-1)],
            [0, - wI, R1 * E * z**(-1)],
            [z**(-1), - wI * z**(-1), (R1 * E + R2) * z**(-2)]
        ])
    elif sys_type == "22_mult_unt":
        P = Ss.tf_blocks([
            [0, 0, R1 * z**(-N_in)],
            [0, 0, R1 * E * wI * z**(-N_in)],
            [z**(-N_out), z**(-N_out), (R1 * E + R2) * z**(-N_in - N_out)]
        ])
    elif sys_type == "33_mult_unt":
        P = Ss.tf_blocks([
            [0, 0, 0, R1 * z**(- N_in)],
            [0, 0, - E_gain * wI, R1 * sumz * E_gain * wI * z**(-N_in)],
            [0, 1, - E_gain, R1 * sumz * E_gain * z**(-N_in)],
            [z**(-N_in), z**(-N_in), - E_gain * z**(-N_in),
             (R1 * sumz * E_gain + R2) * z ** (-N_in - N_out)]
        ])
        pass
    return P


def analysis(plant, controller, analysis_dict, controller_name='noname'):
    """A basic analysis of a plant/controller pair.

    The plant and controller are first combined using a Linear
    Fractional Transformation to form the closed-loop mapping from
    exogeneous input to exogeneous output. Analysis are then performed
    on this closed-loop mapping.

    Signals are enumerated from 0.

    Args:
        plant: A (3 outputs, 2 inputs) discrete-time LTI system.
        controller: A (1 output, 1 input) discrete-time LTI system.
        analysis_dict: A dictionary that contains descriptions of the
            different analyses to conduct.

    """
    H = Ss.lft(plant, controller)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(
        plant, nu=1, ny=1)
    Ts = plant.dt

    # stability
    w, q = np.linalg.eig(H.A)
    wmax = w[np.argmax(np.abs(w))]
    if np.abs(wmax) > 1:
        print(" -- Closed-loop system UNSTABLE. w_max={:}".format(wmax))
    else:
        print(" -- Closed-loop system STABLE. w_max={:}".format(wmax))

    Nsteps = 800
    freqs = analysis_dict['freqs']
    T_sim = np.arange(Nsteps) * Ts
    nrow, ncol = analysis_dict['row_col']
    fig, axs = plt.subplots(nrow, ncol, figsize=(10, 10))
    vsys = analysis_dict['virtual_sys']
    H_vsys = co.c2d(
        co.tf([vsys['k_E'], 0], [vsys['m'], vsys['b'], vsys['k'] + vsys['k_E']]), Ts)
    for entry in analysis_dict['recipe']:
        i, j = entry[0], entry[1]

        plot_type = entry[2]

        if plot_type == 'q':
            Q = co.feedback(controller, Pyu, sign=1)
            _, Q_imp = co.impulse_response(Q, np.arange(Nsteps) * Ts)
            Q_imp = Q_imp.flatten()
            axs[i, j].plot(T_sim, Q_imp)
            axs[i, j].set_title(r"$Q (\delta[n]) $")

        elif plot_type == 'step_sim':
            # the transfer function of the ideal response of the mapping from
            # position xd to velocity v is:
            # s k_E / (ms^2 + bs + k + k_E)
            _, y_ideal = co.step_response(H_vsys, T_sim)
            axs[i, j].plot(T_sim, y_ideal[0, :], label='ideal resp.')
            axs[i, j].set_title('Step response')

        elif plot_type == 'step_sim_int':
            H_model = co.c2d(
                co.tf([vsys['k_E'], 0], [vsys['m'], vsys['b'], vsys['k'] + vsys['k_E']]), Ts)
            _, y_ideal = co.step_response(H_model, T_sim)
            y_ideal_int = np.cumsum(y_ideal[0, :]) * Ts
            axs[i, j].plot(T_sim, y_ideal_int, label='ideal resp. int.')
            axs[i, j].set_title('Step response')

        elif plot_type == 'step':
            # format: i, j, 'step', (out_idx, in_idx), (out_idx2, in_idx2)
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.step_response(H[output_idx, input_idx], T_sim)
                axs[i, j].plot(T_sim, y_step[0, :],
                               label='step resp. H{:d}{:d}'.format(output_idx, input_idx))
                axs[i, j].set_title('Step response')

        elif plot_type == 'impulse':
            # format: i, j, 'step', (out_idx, in_idx), (out_idx2, in_idx2)
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.impulse_response(
                    H[output_idx, input_idx], T_sim)
                axs[i, j].plot(T_sim, y_step[0, :],
                               label='imp. resp. H{:d}{:d}'.format(output_idx, input_idx))
                axs[i, j].set_title('Impulse response')

        elif plot_type == 'step_int':
            # format: i, j, 'step_int', (out_idx, in_idx), (out_idx2, in_idx2)
            # the integral of the step response
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.step_response(H[output_idx, input_idx], T_sim)
                y_step_int = np.cumsum(y_step[0, :]) * Ts
                axs[i, j].plot(T_sim, y_step_int,
                               label='step resp. int. H{:d}{:d}'.format(output_idx, input_idx))
                axs[i, j].set_title('Step response integral')

        elif plot_type == 'nyquist':
            # format: i, j, 'nyquist', (out_idx, in_idx), (w0, w1, w3) [these
            # are interested frequencies]
            output_idx, input_idx = entry[3]
            mag, phase, freqs = H[output_idx, input_idx].freqresp(freqs)
            H = mag[0, 0] * np.exp(1j * phase[0, 0])
            axs[i, j].plot(H.real, H.imag, '-',
                           label='H{:d}{:d}'.format(output_idx, input_idx))

            if len(entry) > 4:
                toplot_idx = []
                for omega in entry[4]:
                    idx = np.argmin(np.abs(freqs - omega))
                    axs[i, j].text(H[idx].real, H[idx].imag,
                                   "{:.3f} rad/s".format(freqs[idx]))
                    toplot_idx.append(idx)
                axs[i, j].scatter(H[toplot_idx].real, H[toplot_idx].imag)

            axs[i, j].set_aspect('equal')

        elif plot_type == 'bode_mag':
            # format: i, j, 'bode_mag', (out_idx, in_idx), (out_idx2, in_idx2)
            for k in range(3, len(entry)):
                if len(entry[k]) == 2:
                    output_idx, input_idx = entry[k]
                    mag, phase, freqs = H[output_idx,
                                          input_idx].freqresp(freqs)
                    label = 'H{:d}{:d}'.format(output_idx, input_idx)
                elif entry[k] == 'vsys':
                    H_model = co.c2d(
                        co.tf([vsys['k_E'], 0], [vsys['m'], vsys['b'], vsys['k'] + vsys['k_E']]), Ts)
                    mag, phase, freqs = H_model.freqresp(freqs)
                    label = 'vsys'
                axs[i, j].plot(freqs, mag[0, 0], label=label)

            axs[i, j].set_xscale('log')
            axs[i, j].set_yscale('log')
            axs[i, j].set_xlabel('Freq(rad/s)')
            axs[i, j].set_title("Frequency Response (non-db)")

        elif plot_type == 'bode_phs':
            # format: i, j, 'bode_mag', (out_idx, in_idx), (out_idx2, in_idx2)

            for k in range(3, len(entry)):
                if len(entry[k]) == 2:
                    output_idx, input_idx = entry[k]
                    mag, phase, freqs = H[output_idx,
                                          input_idx].freqresp(freqs)
                    label = 'H{:d}{:d}'.format(output_idx, input_idx)
                elif entry[k] == 'vsys':
                    H_model = co.c2d(
                        co.tf([vsys['k_E'], 0], [vsys['m'], vsys['b'], vsys['k'] + vsys['k_E']]), Ts)
                    mag, phase, freqs = H_model.freqresp(freqs)
                    label = 'vsys'
                axs[i, j].plot(freqs, np.rad2deg(phase[0, 0]), label=label)

            for mult in range(-2, 2):
                axs[i, j].plot([freqs[0], freqs[-1]],
                               [90 * mult, 90 * mult], '--', c='red')

            axs[i, j].set_xscale('log')
            axs[i, j].set_xlabel('Freq(rad/s)')
            axs[i, j].set_title("Phase lag")

    if 'sharex' in analysis_dict.keys():
        for (i1, j1), (i2, j2) in analysis_dict['sharex']:
            axs[i1, j1].get_shared_x_axes().join(axs[i1, j1], axs[i2, j2])

    for i in range(nrow):
        for j in range(ncol):
            axs[i, j].grid()
            axs[i, j].legend()

    fig.suptitle('Analysis plots: {:}'.format(controller_name))
    plt.tight_layout()
    plt.show()


class Qsyn:
    """ Collection of sub-routines used in Q-parametriation based synthesis.
    """
    @staticmethod
    def obtain_time_response_var(
            weight, io_idx, input_kind, T_sim, Pzw, Pzu, Pyw):
        """ Return a time-reponse variable.

        A basis of delayed impulses is assumed. That is:

            Q = weight[0] * 1 + weight[1] * z^-1 + ... + weight[n - 1] * z^-(n - 1)
        """
        print(" -- generate time reponse {:}".format(io_idx))
        i, j = io_idx
        if input_kind == 'step' or input_kind == 'step_int':
            out = co.step_response(Pzw[i, j], T_sim)
        elif input_kind == 'impulse':
            out = co.impulse_response(Pzw[i, j], T_sim)
        Pzw_resp = out[1][0, :]
        resp = Pzw_resp
        Nsteps = len(T_sim)

        PzuPyu_resp = None
        resp_mat = []
        for k in range(weight.shape[0]):
            if k == 0:
                if input_kind == 'step' or input_kind == 'step_int':
                    out = co.step_response((Pzu * Pyw)[i, j], T_sim)
                elif input_kind == 'impulse':
                    out = co.impulse_response((Pzu * Pyw)[i, j], T_sim)
                else:
                    raise(
                        ValueError(
                            "Unknown input kind {:}".format(input_kind)))
                # assign base response
                PzuPyu_resp = out[1][0, :] * Ts
                resp_mat.append(PzuPyu_resp)
                # resp = resp + weight[0] * PzuPyu_resp
            else:
                resp_k = np.zeros_like(PzuPyu_resp)
                resp_k[k:] = PzuPyu_resp[:(Nsteps - k)]
                resp_mat.append(resp_k)
                # resp = resp + weight[k] * resp_k
        resp_mat = np.array(resp_mat).T
        resp = resp + resp_mat * weight
        if input_kind == 'step_int':
            resp = cvx.cumsum(resp) * Ts
        return resp

    def obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw):
        """Return a frequency-response variable.

        The frequency-response variable is the sum of basis
        responses. The basis Q(s) are delayed impulse with magnitude
        Ts.

        """
        print(" -- generate frequency reponse {:}".format(io_idx))
        i, j = io_idx
        mag, phase, _ = Pzw[i, j].freqresp(freqs)
        freqresp = mag[0, 0] * np.exp(1j * phase[0, 0])
        delay_operator = np.exp(- 1j * freqs * Ts)

        PzuPyw_resp = None
        PzuPyw_mat = []
        for k in range(weight.shape[0]):
            if k == 0:
                mag, phase, _ = (Pzu * Pyw)[i, j].freqresp(freqs)
                PzuPyw_resp = mag[0, 0] * np.exp(1j * phase[0, 0]) * Ts
                PzuPyw_mat.append(PzuPyw_resp)
            else:
                delayed_resp = PzuPyw_resp * delay_operator ** k
                PzuPyw_mat.append(delayed_resp)
        PzuPyw_mat = np.array(PzuPyw_mat).T

        # check if the frequency condition is for DC gain
        if len(freqs) == 1 and freqs[0] == 1e-5:
            import ipdb
            ipdb.set_trace()
        freqresp += PzuPyw_mat * weight
        return freqresp

    def lambda_log_interpolate(pts, preview=False):
        """Return a log-log interpolated function.

        Args:
            pts: List of points.
            preview (bool, optional): Plot the interpolated function
                over a wide interval. Only for visualization.

        Returns:
            a lambda function.
        """
        x, y = zip(*pts)
        xlog = np.log10(x)
        ylog = np.log10(y)
        log_interpolator = interp1d(xlog, ylog, fill_value='extrapolate')

        def interpolate_func(xi):
            return np.power(10, log_interpolator(np.log10(xi)))

        if preview:
            xdata = np.logspace(-3, 5, 100)
            ydata = interpolate_func(xdata)
            plt.plot(xdata, ydata)
            plt.scatter(x, y)
            plt.xlim([0, 300])
            plt.show()

        return interpolate_func


def Q_synthesis(Pz_design, specs):
    """Synthesize a controller for the given plant.

    A dictionary is used to supply specification to the synthesizer.

    [recipe]: A list of lists.
    -

    """
    # setup
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(
        Pz_design, nu=1, ny=1)

    Ts = Pz_design.dt
    T_sim = np.arange(specs['Nsteps']) * Ts

    z = co.tf([1, 0], [1], Ts)
    freqs = specs['freqs']

    # synthesis
    weight = cvx.Variable(specs['Ntaps'])
    constraints = []

    # objective: time-domain based objective
    input_kind, io_idx, sys_desired, obj_weight = specs['objective']
    imp_var = Qsyn.obtain_time_response_var(
        weight, io_idx, input_kind, T_sim, Pzw, Pzu, Pyw)
    if input_kind == 'step':
        out = co.step_response(sys_desired, T_sim)
    elif input_kind == 'impulse':
        out = co.impulse_response(sys_desired, T_sim)
    elif input_kind == 'step_int':
        out = co.step_response(sys_desired, T_sim)
        out_1 = np.cumsum(out[1], axis=1) * Ts
        out = (0, out_1)  # TODO: bad hack, improve this
    imp_desired = out[1][0, :]
    imp_desired[specs['resp_delay']:] = (imp_desired[:specs['Nsteps'] - specs['resp_delay']])
    imp_desired[:specs['resp_delay']] = 0
    cost1 = obj_weight * cvx.norm(imp_var - imp_desired)

    # objective: time-domain quadratic regulator
    if 'objective-qreg' in specs:
        input_kind, io_idx, target, obj_weight = specs['objective-qreg']
        imp_var = Qsyn.obtain_time_response_var(
            weight, io_idx, input_kind, T_sim, Pzw, Pzu, Pyw)
        if input_kind == 'step':
            out = co.step_response(sys_desired, T_sim)
        elif input_kind == 'impulse':
            out = co.impulse_response(sys_desired, T_sim)
        elif input_kind == 'step_int':
            out = co.step_response(sys_desired, T_sim)
            out_1 = np.cumsum(out[1], axis=1) * Ts
            out = (0, out_1)  # TODO: bad hack, improve this
        cost1 += obj_weight * cvx.norm(imp_var - target)

        # no overshoot constraint
        constraints.append(imp_var <= target)

    # frequency-domain constraints
    freq_vars = {}
    for io_idx, func in specs['freq-constraints']:
        freq_var = Qsyn.obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw)
        constraints.append(
            cvx.abs(freq_var) <= func(freqs)
        )
        freq_vars[io_idx] = freq_var

    # min H-infinity norm
    cost_inf = 0
    for io_idx, func, obj_weight in specs['objective-Hinf']:
        if io_idx in freq_vars:
            freq_var = freq_vars[io_idx]
        else:
            freq_var = Qsyn.obtain_freq_var(
                weight, io_idx, freqs, Pzw, Pzu, Pyw)
            freq_vars[io_idx] = freq_var
        cost_inf += obj_weight * cvx.norm(freq_var - func(freqs), 'inf')

    # passivity
    for io_idx, wc in specs['passivity']:
        if io_idx in freq_vars:
            freq_var = freq_vars[io_idx]
        else:
            freq_var = Qsyn.obtain_freq_var(
                weight, io_idx, freqs, Pzw, Pzu, Pyw)
            freq_vars[io_idx] = freq_var
        # freqs[idx_wc] is the smallest rotationl velocity that is greater than
        # wc
        idx_wc = 0
        while freqs[idx_wc] < wc:
            idx_wc += 1
        assert(freqs[idx_wc] > wc)
        constraints.append(
            cvx.real(freq_var[:idx_wc]) >= 0
        )
        print(
            "Passivity constraint accounted for with wc={:f}".format(
                freqs[idx_wc]))

    # dc-gain
    for io_idx, gain in specs['dc-gain']:
        # only w=0 is needed
        imp_var_ = Qsyn.obtain_time_response_var(
            weight, io_idx, 'impulse', T_sim, Pzw, Pzu, Pyw)

        constraints.append(
            cvx.sum(imp_var_) == gain
        )

    # # noise attenuation:
    # constraints.append(
    #     cvx.abs(H00_freqrp) <= 2e-2
    # )

    # distance moved should never be negative (not effective, hence, not in use)
    # dist_mat = np.zeros((Nstep, Nstep))
    # for i in range(Nstep):
        # dist_mat[i, :i] = 1.0
    # constraints.append(dist_mat * imp_var >= 0)

    # zero distance travelled  (not effective) (this constraint is satisfied by definition)
    # constraints.append(cvx.sum(imp_var) == 0)

    # # good (positive phase lag) passive until w = 4 (not effective)
    # omega_ths_idx = np.argmin(np.abs(freqs - 8.0))
    # constraints.append(
    #     cvx.imag(H00_freqrp[:omega_ths_idx]) >= 0
    # )

    # penalize rapid changes in Q
    diffQ = np.zeros((specs['Ntaps'] - 1, specs['Ntaps']))
    for i in range(specs['Ntaps'] - 1):
        diffQ[i, i: i + 2] = [1.0, -1.0]

    diffY = np.zeros((specs['Nsteps'] - 1, specs['Nsteps']))
    for i in range(specs['Nsteps'] - 1):
        diffY[i, i: i + 2] = [1.0, -1.0]

    imp_weight = np.ones(specs['Nsteps'])
    cost2 = specs['reg'] * cvx.norm1(weight * Ts)
    cost3 = specs['reg_diff_Q'] * cvx.norm(diffQ * weight * Ts, p=2)
    cost4 = specs['reg_diff_Y'] * cvx.norm(diffY * imp_var)
    cost5 = specs['reg'] * cvx.norm1(imp_var - imp_desired)
    cost = cost1 + cost2 + cost3 + cost4 + cost5 + cost_inf
    print(" -- Form cvxpy Problem instance")
    prob = cvx.Problem(cvx.Minimize(cost), constraints)
    print(" -- Start solving")
    prob.solve(verbose=True)
    assert(prob.status == 'optimal')

    print("cost: cost1={:f}, cost2={:f},\n      cost3={:f} cost4={:f}"
          "\n     cost5={:f} cost-inf={:f}" .format(
              cost1.value, cost2.value, cost3.value, cost4.value, cost5.value,
              cost_inf.value))
    print("Explanations:\n"
          "cost1: |y[n] - y_ideal[n]|_2\n"
          "cost2: |Q|_1\n"
          "cost3: |DIFF * Q|_2\n"
          "cost4: |DIFF * y[n]|_2 (for smooth movement)\n"
          "where DIFF is the difference matrix.")

    # debug
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(imp_var.value, label="imp_var")
    axs[0, 0].plot(imp_desired, label="imp_desired")
    axs[0, 0].plot(imp_weight * 0.004, label="scaled weight")
    axs[0, 0].legend()
    axs[0, 1].plot(weight.value, label="weight")
    _i = 0
    for io_idx, func in specs['freq-constraints']:
        freq_var = freq_vars[io_idx]
        axs[1, 0].plot(
            freqs, np.abs(freq_var.value), label="H{:}".format(io_idx), c='C{:d}'.format(_i))
        axs[1, 0].plot(
            freqs, func(freqs), '--', c='C{:d}'.format(_i))
        _i += 1
    for i, j in [(1, 0), (1, 1)]:
        axs[i, j].grid()
        axs[i, j].legend()
        axs[i, j].set_xscale('log')
        axs[i, j].set_yscale('log')
    plt.show()

    # report result
    # form individual filter
    taps = weight.value * Ts
    Q_fir = co.tf(taps, [1] + [0] * (specs['Ntaps'] - 1), Ts)
    Q = Ss.tf2ss(Q_fir, minreal=True, via_matlab=True)
    K = co.feedback(Q, Pyu, sign=-1)

    return K, {'Qtaps': taps, 'Pyu': Pyu, 'zPyu': z * Pyu}


def desired_impulse(m=1, b=1, k=1, Nstep=600, Ts=0.008):
    # impulse response, comparison with an ideal mass/spring/damper
    H_model = co.c2d(co.tf([1, 0], [m, b, k]), Ts)
    T_imp = np.arange(Nstep) * Ts
    _, y_imp = co.impulse_response(H_model, T_imp)
    return y_imp.flatten()


def print_controller(relative_file_path, Qtaps, zPyu_tf):
    """ Print a 3-dimensional DelayedFeedback controller.

    A delayed-feedback controller is defined as

         feedback(Q, z^-1 * zPyu, sign=-1)

    Every dimension shares the same feedback law.
    """
    num = [float(x) for x in zPyu_tf.num[0][0]]
    den = [float(x) for x in zPyu_tf.den[0][0]]
    num = num + [0] * (len(den) - len(num))
    taps = [float(x) for x in np.array(Qtaps).flatten()]
    ctrl_dict = {
        'xff_taps': taps,
        'xfb_b': num,
        'xfb_a': den,
        'yff_taps': taps,
        'yfb_b': num,
        'yfb_a': den,
        'zff_taps': taps,
        'zfb_b': num,
        'zfb_a': den,
    }
    output_dict = {
        'diagonal3D': {
            'filter_dev': ctrl_dict
        }
    }
    with open(relative_file_path, 'w') as f:
        yaml.dump(output_dict, f, default_flow_style=False)
    print("--> Wrote a controller to {:}".format(relative_file_path))


def main():
    # # impulse response and weight
    # imp_desired = desired_impulse(m=2, b=18, k=5, Nstep=1000)
    # imp_weight = np.zeros(1000)
    # imp_weight[:] = 1
    # plt.plot(imp_desired); plt.show()

    # desired response from w3 to z1
    Pz_design = plantMdelta(
        E_gain=50, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)

    noise_atten_func = Qsyn.lambda_log_interpolate(
        [[0.1, 0.1], [25, 0.1], [25, 0.015], [50, 0.004], [200, 0.002]], preview=True)

    desired_sys = co.c2d(co.tf([50, 0], [2.5, 12, 0 + 50]), Ts)

    def desired_sys_up(freqs, desired_sys=desired_sys):
        return desired_sys.freqresp(freqs)[0][0, 0]

    design_dict = {
        'Ntaps': 800,
        'Nsteps': 1000,
        # 'freqs': np.logspace(-2, np.log10(np.pi / Ts), 1000),
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),
        'resp_delay': 3,  # number of delayed time step

        # different objective
        'objective': ['step_int', (0, 2), desired_sys, 1.5],
        'objective-Hinf': [
            [(0, 2), desired_sys_up, 1.0]
        ],

        # a quadratic regulator objective:
        # (input_kind, input indices, target level, objective weight)
        # 'objective-qreg': ['step_int', (0, 2), 1.0, 1.0],

        # 'objective2': ['inf', (0, 2), co.c2d(co.tf([50, 0], [2.5, 21, 0 + 50]), Ts)]
        # 'objective': ['step', (0, 2), co.c2d(co.tf([50, 0], [3.5, 21, 0 + 50]), Ts)],
        # 'objective': ['step', (0, 2), co.c2d(co.tf([50, 0], [2., 25, 0 + 50]), Ts)],
        # frequencies to apply constraint on

        # 'objective': ['impulse', (0, 0), co.c2d(co.tf([1, 0], [2, 18, 5]), Ts)],
        # the impulse-based objective is dropped in favour of the step objective.
        # The reason is that it is found that the step objectively
        # leads to better performance.

        # cost different paramerter

        # regulation parameter
        'reg': 0e-2,
        'reg_diff_Q': 1e-2,
        'reg_diff_Y': 1,

        # list of ((idxi, idxj), function)
        # frequency-domain constraint: H_ij(e^jwT) <= func(w)
        'freq-constraints': [
            [(1, 1), lambda omega: np.ones_like(omega)],
            [(0, 0), noise_atten_func],
        ],

        # passivity: Re[H_ij(e^jwT)] >= 0, for all w <= w_c
        'passivity': [
            [(0, 0), 26]
        ],

        # DC gain
        'dc-gain': [
            [(2, 2), 0]
        ]
    }

    if input("Design controller?") == 'y':
        K_Qparam, data = Q_synthesis(Pz_design, design_dict)

    if input("Write controller? y/[n]") == 'y':
        print_controller(
            "../config/Nov21_Cart_synthesis_Qparam_synres.yaml",
            data['Qtaps'],
            data['zPyu'])

    analysis_dict = {
        'row_col': (3, 2),
        'freqs': np.logspace(-2, np.log10(np.pi / Ts), 200),
        'sharex': ([(0, 1), (1, 1)], [(0, 0), (1, 0)]),
        'virtual_sys': {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 50},
        'recipe': [
            (0, 0, 'impulse', (0, 2)),
            (2, 0, 'step_int', (0, 2)),
            (2, 0, 'step_sim_int'),
            (1, 0, 'step', (0, 2)),
            (1, 0, 'step_sim'),
            (0, 1, 'bode_mag', (0, 0), (2, 2), (0, 2), 'vsys'),
            (1, 1, 'bode_phs', (0, 0), (2, 2), (0, 2), 'vsys'),
            (2, 1, 'nyquist', (0, 0), [1, 10, 24, 50])
        ]
    }

    Pz_contracted = plantMdelta(
        E_gain=100, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)
    Pz_relaxed = plantMdelta(
        E_gain=20, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)

    if input("Analyze stuffs? y/[n]") == 'y':
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 50}
        analysis(Pz_design, K_Qparam, analysis_dict, controller_name='Qparam')
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 100}
        analysis(
            Pz_contracted,
            K_Qparam,
            analysis_dict,
            controller_name='Qparam')
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 20}
        analysis(Pz_relaxed, K_Qparam, analysis_dict, controller_name='Qparam')

    K_ = {
        'ad_light': co.c2d(co.tf([1], [2.5, 12, 0]), Ts),
        'ad_heavy': co.c2d(co.tf([1], [6, 18, 0]), Ts)
    }

    if input("Analyze Admittance controller") == 'y':
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 50}
        analysis(Pz_design, K_['ad_light'],
                 analysis_dict, controller_name="ad_light")
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 100}
        analysis(Pz_contracted, K_['ad_light'],
                 analysis_dict, controller_name="ad_light")
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 20}
        analysis(Pz_relaxed, K_['ad_light'],
                 analysis_dict, controller_name="ad_light")

    if input("Write controller? y/[n]") == 'y':
        print_controller(
            "../config/Nov21_Cart_synthesis_Qparam_synres.yaml",
            data['Qtaps'],
            data['zPyu'])

    if input("Dump Q controller with pickle? y/[n]") == 'y':
        import pickle
        with open('Nov21_ctrl.pkl', 'wb') as f:
            pickle.dump(K_Qparam, f, -1)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    main()
