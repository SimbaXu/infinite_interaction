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
dt = 0.008


def plant(gain_E=20, m_int=0.1):
    """ An old plant model. Not in use currently.
    """
    z = co.tf([1, 0], [1], dt)
    s = co.tf([1, 0], [1])
    zero = co.tf([0], [1], dt)
    one = co.tf([1], [1], dt)
    R1 = co.c2d(s / (1 + 0.0437 * s), dt)
    E = - co.c2d(gain_E / s, dt)  # exo agent
    R2 = - m_int / dt ** 2 * (1 - z ** (-1)) ** 2  # internally induced

    P = Ss.tf_blocks([[zero, z**(-1) * R1],  # velocity output
                      [zero, z**(-1) * (R1 * E + R2)],  # to form T = L / (1 + L)
                      [one, z**(-1) * (R1 * E + R2)],  # to form S = 1 / (1 + L)
                      [z**(-1), z**(-2) * (R1 * E + R2)]])
    return P


def plantMdelta(E_gain=20, m_int=0.1, wI=0.9995, sys_type='22_mult_unt',
                N_in=1, N_out=1):
    """ Plant Model with multiplicative uncertainty.

    Diagram is drawn in page 4, notebook.

    Args:
        wI (float): Can be a complex function.  Uncertainty in the exogeneous agent.
        sys_type (str, optional): Type of plant.
            - 22_mult_unt: (2, 2) system with multiplicative uncertainty;
            - 22_inv_mult_unt: (2, 2) system with inverse multiplicative uncertainty;
            - 33_mult_unt: (3, 3) system with multiplicative uncertainty.
        N_in (int, optional): Nb. of delay steps in the input path.
        N_out (int, optional): Nb. of delay steps in the output path.
    """
    # constant transfer function
    z = co.tf([1, 0], [1], dt)
    s = co.tf([1, 0], [1])
    zero = co.tf([0], [1], dt)
    one = co.tf([1], [1], dt)

    E_gain = - abs(E_gain)

    # basic blocks
    R1 = co.c2d(s / (1 + 0.0437 * s), dt)
    sumz = co.c2d(1 / s, dt)  # discrete-time integrator
    E = E_gain * sumz  # exo agent
    R2 = - m_int / dt ** 2 * (1 - z ** (-1)) ** 2  # internally induced

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
            [0, 0, E_gain * wI, R1 * sumz * E_gain * wI * z**(-N_in)],
            [0, 1, E_gain, R1 * sumz * E_gain * z**(-N_in)],
            [z**(-N_in), z**(-N_in), E_gain * z**(-N_in), (R1 * sumz * E_gain + R2) * z **(-N_in - N_out)]
        ])
        pass
    return P


def analysis(plant, controller, controller_name='noname',
             m=0.5, b=10, k=80,
             freqs_bnd_yn=[1e-2, 255], mag_bnd_yn=[-10, -10],
             freqs_bnd_T=[1e-2, 357], mag_bnd_T=[6, 6]):
    """A basic analysis of a plant/controller pair.

    Simulate the system and draw several plots.

      1 |  3
    ----|----
      2 |  4

    1): input output 1
    2): input output 2
    3): frequency responses
    4): nyquist

    Args:
        plant: A (3 outputs, 2 inputs) discrete-time LTI system.
        controller: A (1 output, 1 input) discrete-time LTI system.
        internal_data: The internal responses {R, M, N, L, H} in that
                       order. Or can be None.
    """
    H = Ss.lft(plant, controller)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(plant, nu=1, ny=1)
    dT = plant.dt
    Nsteps = 800

    # stability
    w, q = np.linalg.eig(H.A)
    wmax = w[np.argmax(np.abs(w))]
    if np.abs(wmax) > 1:
        print(" -- Closed-loop system UNSTABLE. w_max={:}".format(wmax))
    else:
        print(" -- Closed-loop system STABLE. w_max={:}".format(wmax))

    # spectral analysis
    w_nyquist = np.pi / dT
    freqs = np.logspace(-2, np.log10(w_nyquist) - 1e-2, 1000)
    mag, phase, omega = H.freqresp(freqs)

    # step and impulse response
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Q-params
    Q = co.feedback(controller, Pyu, sign=1)
    _, Q_imp = co.impulse_response(Q, np.arange(Nsteps) * dT)
    Q_imp = Q_imp.flatten()
    axs[0, 0].plot(Q_imp)

    # impulse response, comparison with an ideal mass/spring/damper
    H_model = co.c2d(co.tf([1, 0], [m, b, k]), dT)
    T_imp = np.arange(Nsteps) * dT
    _, y_imp = co.impulse_response(H[0, 0], T_imp)
    _, y_imp_model = co.impulse_response(H_model, T_imp)

    # step response
    _, y_step = co.step_response(H[0, 0], T_imp)

    axs[1, 0].plot(T_imp, y_step[0, :], label='step response')
    axs[1, 0].plot(T_imp, y_imp[0, :], label='h1[n]')
    axs[1, 0].plot(T_imp, y_imp_model[0, :],
                   label='vel-impulse[sys(m={:},b={:},k={:})]'.format(m, b, k))
    axs[1, 0].legend()
    axs[1, 0].text(5, 0.0003, 'model(m={:},b={:},k={:})'.format(
        m, b, k), horizontalalignment='center')

    # frequency responses
    mag_yn = 20 * np.log10(mag[0, 0])
    mag_M = 20 * np.log10(mag[1, 1])
    try:
        mag_33 = 20 * np.log10(mag[2, 2])  # mapping from intention displacement to force feedback
        axs[0, 1].plot(freqs, mag_33, label='H33(e^jw)', c='C2')
    except:
        pass
    # bounds on H_yn and H_T
    axs[0, 1].plot(freqs, mag_yn, label='|H11(e^jw)|', c='C0')
    axs[0, 1].plot(freqs, mag_M, label='|H22(e^jw)|', c='C1')
    axs[0, 1].plot([w_nyquist, w_nyquist], [
                   np.min(mag_yn), np.max(mag_yn)], '--', c='red')
    axs[0, 1].plot(freqs_bnd_yn, mag_bnd_yn, 'x--', c='C0', label='wN(w)^-1')
    axs[0, 1].plot(freqs_bnd_T, mag_bnd_T, 'x--', c='C1', label='wT(w)^-1')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_ylabel('Mag(dB)')
    axs[0, 1].set_xlabel('Freq(rad/s)')
    axs[0, 1].set_title("Frequency Response")
    axs[0, 1].legend()
    axs[0, 1].grid()

    # # nyquist plot of H_yn (noise to output)
    # H_yn = mag[0, 0] * np.exp(1j * phase[0, 0])
    # axs[1, 1].plot(H_yn.real, H_yn.imag, '-')
    # axs[1, 1].set_title("Nyquist plot of H_yn(s)")
    # axs[1, 1].grid()

    # int_omegas = [6, 10, 50, 80]  # interested angular velocity
    # for int_omega_ in int_omegas:
    #     idx = np.argmin(np.abs(freqs - int_omega_))
    #     axs[1, 1].text(H_yn[idx].real, H_yn[idx].imag, "{:.3f} rad/s".format(freqs[idx]))

    # Axs[1,1]: phase lag plot
    axs[1, 1].set_title("phase lag")
    axs[1, 1].plot(freqs, phase[2, 2], label='angle{H22}')
    axs[1, 1].set_xscale('log')
    axs[1, 1].grid()
    axs[1, 1].legend()

    fig.suptitle('Analysis plots: {:}'.format(controller_name))
    plt.tight_layout()
    plt.show()


def Q_synthesis(Pz_design, imp_weight=np.ones(500), Ntaps=300, Nstep=500, imp_desired=None, reg=1e-5):
    """
    """
    # setup
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(Pz_design, nu=1, ny=1)

    dt = Pz_design.dt
    z = co.tf([1, 0], [1], dt)
    zero = co.tf([0], [1], dt)

    # check and change impulse desired shape
    if imp_desired is None:
        raise(ValueError("No desired impulse supplied!"))
    if imp_desired.shape[0] < Nstep:
        imp_desired = np.concatenate(
            (imp_desired, np.zeros(Nstep - imp_desired.shape[0])))
    else:
        imp_desired = imp_desired[:Nstep]

    # basis Q elements
    print("--> Form basis closed-loop output mappings")
    # basis_Qs = [zero, dt / (1 - z**(-1))]  # (0, accumulator)
    basis_Qs = [zero, zero]  # (0, accumulator)
    # delay
    for i in range(Ntaps):
        den = [1] + [0] * i
        basis_Qs.append(co.tf([dt], den, dt))

    # basis output mappings H
    print("--> Form basis closed-loop output mappings")
    basis_Hs = []
    for Q in basis_Qs:
        H = Pzw + Pzu * Q * Pyw
        basis_Hs.append(H)

    # basis impulse response
    print("--> Compute basis impulse responses")
    basis_imps = []
    # zero and accumulator
    for i in range(2):
        H = basis_Hs[i]
        _, imp = co.impulse_response(H[0, 0], np.arange(Nstep) * dt)
        basis_imps.append(imp.flatten())

    # impulse response: This code assume that  the last impulse are all finite
    imp = None
    for i in range(Ntaps):
        if i == 0:
            H = basis_Hs[2]
            _, imp = co.impulse_response(H[0, 0], np.arange(Nstep) * dt)
            imp = imp.flatten()
            basis_imps.append(imp)
        else:
            imp_ = np.zeros(Nstep)
            imp_[i:] = imp[:(Nstep - i)]
            basis_imps.append(imp_)

    # basis frequency response
    print("--> Compute basis frequency responses")
    omega_nyquist = np.pi / dt
    omegas = np.linspace(20, omega_nyquist, 100)
    omegas = np.concatenate((np.linspace(1e-2, 19.9, 100), omegas))
    basis_H00_freqrp = []
    basis_H11_freqrp = []
    basis_H22_freqrp = []
    for H in basis_Hs:
        # H00
        mag, phase, _ = H[0, 0].freqresp(omegas)
        freqresp = mag[0, 0] * np.exp(1j * phase[0, 0])
        basis_H00_freqrp.append(freqresp)
        # H11
        mag, phase, _ = H[1, 1].freqresp(omegas)
        freqresp = mag[0, 0] * np.exp(1j * phase[0, 0])
        basis_H11_freqrp.append(freqresp)
        # H22
        mag, phase, _ = H[2, 2].freqresp(omegas)
        freqresp = mag[0, 0] * np.exp(1j * phase[0, 0])
        basis_H22_freqrp.append(freqresp)

    # ## debug impulse respnse
    # for i in range(5):
    #     plt.plot(basis_impresps[i], 'x-', label='{:}'.format(i))
    # for i in range(0, Ntaps, 50):
    #     plt.plot(basis_impresps[i], 'x-', label='{:}'.format(i))
    # plt.legend()
    # plt.show()

    # synthesis
    Nvar = len(basis_Qs)
    weight = cvx.Variable(Nvar)
    constraints = [
        cvx.sum(weight) == 1
    ]

    imp_var = 0
    H00_freqrp = 0
    H11_freqrp = 0
    for i in range(Nvar):
        imp_var += weight[i] * basis_imps[i]
        H00_freqrp += weight[i] * basis_H00_freqrp[i]
        H11_freqrp += weight[i] * basis_H11_freqrp[i]

    try:
        H22_freqrp = 0
        for i in range(Nvar):
            H22_freqrp += weight[i] * basis_H22_freqrp[i]
        # constraints.append(cvx.abs(H22_freqrp) <= 55)
    except Exception as e:
        print("Unable to constraint H22, error: {:}".format(e))

    # robust stability
    constraints.append(
        cvx.abs(H11_freqrp) <= 1
    )

    # # noise attenuation:
    # constraints.append(
    #     cvx.abs(H00_freqrp) <= 2e-2
    # )

    # constraint disp -> force

    # distance moved should never be negative (not effective, hence, not in use)
    # dist_mat = np.zeros((Nstep, Nstep))
    # for i in range(Nstep):
        # dist_mat[i, :i] = 1.0
    # constraints.append(dist_mat * imp_var >= 0)

    # zero distance travelled  (not effective)
    # constraints.append(cvx.sum(imp_var) == 0)

    # penalize rapid changes in impulse response
    imp_diff_mat = np.zeros((Nstep - 1, Nstep))
    for i in range(Nstep - 1):
        imp_diff_mat[i, i: i + 2] = [1.0, -1.0]

    imp_weight_mat = np.diag(imp_weight)
    cost1 = cvx.norm(imp_weight_mat * (imp_var - imp_desired))
    cost2 = reg * cvx.norm1(weight)
    cost3 = cvx.norm(imp_diff_mat * imp_var)
    cost = cost1 + cost2 + cost3
    prob = cvx.Problem(cvx.Minimize(cost), constraints)
    prob.solve(verbose=True)
    assert(prob.status == 'optimal')

    # debug
    print("cost: mse={:f}, reg={:f}, mse diff={:f}".format(cost1.value, cost2.value, cost3.value))
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(imp_var.value, label="imp_var")
    axs[0, 0].plot(imp_desired, label="imp_desired")
    axs[0, 0].plot(imp_weight * 0.004, label="scaled weight")
    axs[0, 0].legend()
    axs[0, 1].plot(weight.value, label="weight")
    axs[1, 0].plot(omegas, np.abs(H00_freqrp.value), label="H00")
    axs[1, 1].plot(omegas, np.abs(H11_freqrp.value), label="H11")
    axs[1, 1].plot(omegas, np.abs(H22_freqrp.value), label="H22")
    for i, j in [(1, 0), (1, 1)]:
        axs[i, j].grid()
        axs[i, j].legend()
        axs[i, j].set_xscale('log')
        axs[i, j].set_yscale('log')
    plt.show()

    # report result
    # form individual filter
    taps = weight.value[2:] * dt
    Q_fir = co.tf(taps, [1] + [0] * (Ntaps - 1), dt)
    Q = Ss.mtf2ss(Q_fir, minreal=True)
    K = co.feedback(Q, Pyu, sign=-1)

    return K, {'Qtaps': taps, 'Pyu': Pyu, 'zPyu': z * Pyu}


def desired_impulse(m=1, b=1, k=1, Nstep=600, dT=0.008):
    # impulse response, comparison with an ideal mass/spring/damper
    H_model = co.c2d(co.tf([1, 0], [m, b, k]), dT)
    T_imp = np.arange(Nstep) * dT
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

    K_ = {
        'admittance': co.c2d(co.tf([1], [6, 18, 15]), dt)
    }

    # impulse response and weight
    imp_desired = desired_impulse(m=2, b=20, k=5, Nstep=1000)
    imp_weight = np.zeros(1000)
    imp_weight[:] = 1
    plt.plot(imp_desired); plt.show()

    Pz_design = plantMdelta(
        E_gain=70, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)
    # Pz_design = plantMdelta(E_gain=50, wI=1.0, sys_type='22_mult_unt', m_int=0.15, N_in=2, N_out=1)

    K_Qparam, data = Q_synthesis(
        Pz_design, imp_desired=imp_desired, imp_weight=imp_weight,
        Ntaps=800, Nstep=1000)

    # analysis Qparam
    if input("Analyze stuffs? y/[n]") == 'y':
        analysis(Pz_design, K_Qparam, m=4, b=12, k=0, controller_name='design plant')

        Pz_contract = plantMdelta(E_gain=50, N_in=2, N_out=2, sys_type='33_mult_unt')
        analysis(Pz_contract, K_Qparam, m=4, b=12, k=0, controller_name='contacting (high stiffness)')

        Pz_open = plantMdelta(E_gain=60)
        analysis(Pz_open, K_Qparam, m=4, b=12, k=0, controller_name='open (low stiffness)')


        analysis(Pz_design, K_['admittance'], m=4, b=12, k=0, controller_name="admittance")

    if input("Write controller? y/[n]") == 'y':
        print_controller("../config/Nov21_Cart_synthesis_Qparam_synres.yaml", data['Qtaps'], data['zPyu'])

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()

if __name__ == '__main__':
    main()
