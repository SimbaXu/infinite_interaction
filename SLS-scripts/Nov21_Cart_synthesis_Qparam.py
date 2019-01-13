"""This script synthesize a 1-dof controller for 1 direction of the
3-dimensional Cartersian admittance control problem.
"""
import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import yaml

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

    k_nom = 30
    m_des = 1.5
    b_des = 8

    # desired response from w3 to z1
    Pz_design = plantMdelta(
        E_gain=k_nom, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=2)

    noise_atten_func = Ss.Qsyn.lambda_log_interpolate(
        [[0.1, 0.1], [10, 0.1], [20, 0.02], [25, 0.01], [30, 0.01], [50, 0.005], [200, 0.002]], preview=True)

    desired_sys = co.c2d(co.tf([m_des * k_nom, b_des * k_nom, 0], [m_des, b_des, 0 + k_nom]), Ts)

    def desired_sys_up(freqs, desired_sys=desired_sys):
        return desired_sys.freqresp(freqs)[0][0, 0]

    design_dict = {
        'ny': 1,
        'nu': 1,
        'Ntaps': 800,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),

        # different objective
        'shape-time-delay': 2,  # number of delayed time step
        'shape-time': [
            ['step', (2, 2), desired_sys, 5]
        ],
        'shape-freq': [
            [(2, 2), desired_sys, 'inf', 1.0]
        ],
        'constraint-freq': [
            # [(1, 1), lambda omega: np.ones_like(omega)],
            [(0, 0), noise_atten_func],
        ],

        # regulation: reg2 * ||Q||_2 + reginf * ||Q||_inf
        'reg2': 0,
        'reginf': 0,

        # list of ((idxi, idxj), function)
        # frequency-domain constraint: H_ij(e^jwT) <= func(w)


        # passivity: Re[H_ij(e^jwT)] >= 0, for all w <= w_c
        'passivity': [
            [(0, 0), 0.001, 20]
        ],

        # DC gain
        'dc-gain': [
            [(2, 2), 0]
        ]
    }

    if input("Design controller?") == 'y':
        syn_data = Ss.Qsyn.Q_synthesis(Pz_design, design_dict)

    if input("Write controller? y/[n]") == 'y':
        print_controller(
            "../config/Nov21_Cart_synthesis_Qparam_synres.yaml",
            syn_data['Qtaps'],
            syn_data['zPyu'])

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

    m_test = 0.1
    Pz_nom = plantMdelta(
        E_gain=k_nom, wI=1.0, sys_type='33_mult_unt', m_int=m_test, N_in=1, N_out=2,
        tau_R1=0.0437)

    Pz_contracted = plantMdelta(
        E_gain=100, wI=1.0, sys_type='33_mult_unt', m_int=m_test, N_in=1, N_out=2)
    Pz_relaxed = plantMdelta(
        E_gain=20, wI=1.0, sys_type='33_mult_unt', m_int=m_test, N_in=1, N_out=2)

    if input("Analyze stuffs? y/[n]") == 'y':
        # form Q param controller (state-space)
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            syn_data['Qtaps'], syn_data['Pyu'])

        # start analysis
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 50}
        analysis(Pz_nom, K_Qparam_ss, analysis_dict, controller_name='Qparam')

        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 100}
        analysis(
            Pz_contracted,
            K_Qparam_ss,
            analysis_dict,
            controller_name='Qparam')

        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 20}
        analysis(Pz_relaxed, K_Qparam_ss, analysis_dict, controller_name='Qparam')

    K_ = {
        'ad_light': co.c2d(co.tf([1], [2.5, 12, 0]), Ts),
        'ad_heavy': co.c2d(co.tf([1], [6, 18, 0]), Ts)
    }

    if input("Analyze Admittance controller") == 'y':
        analysis_dict['virtual_sys'] = {'m': 2.5, 'b': 12, 'k': 0, 'k_E': 50}
        analysis(Pz_nom, K_['ad_light'],
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
            syn_data['Qtaps'],
            syn_data['zPyu'])

    if input("Dump Q controller with pickle? y/[n]") == 'y':
        import pickle
        with open('Nov21_ctrl.pkl', 'wb') as f:
            pickle.dump(K_Qparam_ss, f, -1)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    main()
