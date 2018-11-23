import numpy as np
import control as co
import matplotlib.pyplot as plt
import cvxpy as cvx
import SLSsyn as Ss
import yaml
import os.path
from scipy.linalg import block_diag
import scipy.sparse as sparse

dT = 0.008
s = co.tf([1, 0], [1])


class OldControllers:
    """ Controllers that I designed sometime ago using different methods.
    """
    A1 = co.c2d(co.tf([1], [1, 6, 5]), dT)
    Alp = co.c2d(1.0 / 5 / (0.3 * s + 1) / (0.5 * s + 1)
                 * (0.032 * s + 1) * (0.0021 * s + 1), dT)
    Ainf = co.tf([0, 0.000234247581937, -
                  0.001519129809213, 0.004117723048895, -
                  0.005971926120458, 0.004887991561889, -
                  0.002140996863836, 0.000392090600785],
                 [1.0, - 7.692092610286219, 25.432665781511435, -
                  46.858394653166798, 51.963159172009036, -
                  34.686122493591625, 12.905678916635351, -
                  2.064894113111176], dT)

    Ac14 = co.tf([0, 0.000780043976221, -
                  0.001879759125377, 0.001518988906942, -
                  0.000411374804692, 0],
                 [1.0, -
                  2.915259796026691, 2.835076504158775, -
                  0.919776846642308, 0.000000000000039, 0], dT)


def plant(Larm=0.4, Hgain=20, Hdelay=0.01):
    """Nominal Model as defined in Oct02_SISO_model.m

    Return a linear model of the admittance control robot.
    The model has the following form:

    [d  ] = [P11 P12] [fH]
    [fHT]   [P21 P22] [qu]
    [mm ]   [P31 P32]

    d:   position of the end-effector
    fHT: force output, does not have physical meaning, only use to
    mm: effective torque measured,

    fH: force exerted by human,
    qu: control generated by the admittance controller,

    Consider a controller that maps from [mm] to [qu], the closed-loop
    response would be a mapping from [fH] -> [d, fHT].

    Check org_imgs/Sysdiagram0ax.pdf

    """
    s = co.tf([1, 0], [1])
    R1 = (-s + 55.56) / (s + 55.56) / (0.0437 * s + 1)
    # NOTE: the near last two terms is the first-order pade
    # approximation of the transport delay exp(-sHdelay).
    # NOTE: the last two terms: 10 / (10 + s) encodes the fact that
    # human is very bad at following anything faster than 10
    # rad/s. These terms essentially mean unit gain at low frequency
    # and quick drop afterward.
    H = Hgain * (1 - Hdelay * s) / (1 + Hdelay * s) * 10 / (10 + s)
    Se = 2 * np.pi * 73 / (s + 2 * np.pi * 73)
    R2 = 0.475 * s ** 2 / (s / 100 + 1) ** 2 * 35 / \
        (s + 35) * (-s + 66) / (s + 66)

    P = Ss.tf_blocks([[0, R1 * Larm],
                      [0, (R1 * Larm * H + R2)],
                      [Larm * Se, - (R1 * Larm * H + R2) * Se * Larm]])
    return P


def analysis(plant, controller, Mp=1.05, Tr=0.9, controller_name='noname',
             internal_data=None, m=0.5, b=10, k=80,
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
    w, q = np.linalg.eig(H.A)
    wmax = w[np.argmax(np.abs(w))]
    if np.abs(wmax) > 1:
        print(" -- Closed-loop system UNSTABLE. w_max={:}".format(wmax))
    else:
        print(" -- Closed-loop system STABLE. w_max={:}".format(wmax))

    # spectral analysis
    w_nyquist = np.pi / dT
    freqs = np.logspace(-2, np.log10(w_nyquist) - 1e-2, 100)
    mag, phase, omega = H.freqresp(freqs)
    # step and impulse response
    T, y_step = co.step_response(H, [0, 10])
    # NOTE1: if T = [0, 10], the algorithm would not compute the
    # impulse response, but the response to a linearly decreasing
    # input from 1 at time 0 to 0 at time 10. To address this issue,
    # input the whole time vector.

    # NOTE2: co.impulse is not identical to matlab impulse. The
    # difference is that co.impulse compute the reponse w.r.t to a
    # input of the form [1, 0, ...], while matlab returns the reponse
    # to [1 / Ts, 0, ...]. Matlab's convention would produce similar
    # impulse response for two corresponding c-time and d-time. co's
    # convention is more theoretical.

    T_step = np.arange(y_step.shape[1]) * dT  # get correct T
    dss = y_step[0, -1]
    # const-jerk input: accelerate to 5 N in 0.3 sec, then deccelerate to 1 in
    # 0.3 sec
    T_forced, F_forced = const_jerk_input(
        duration=10, accel_duration=0.6, f_max=4, f_ss=1.0, dT=dT)
    T_forced, Y_forced, _ = co.forced_response(H, T_forced, F_forced)

    # step response
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    axs[0, 0].plot(T_step, y_step[0, :], label='Hd(t)*u(t)')
    axs[0, 0].plot(T_forced, Y_forced[0, :], label='resp to const-jerk')
    axs[0, 0].plot([0, Tr, Tr, 10], [0, 0, 0.98 *
                                     dss, 0.98 * dss], '--', c='gray')
    axs[0, 0].plot([0, 10], [Mp * dss, Mp * dss], '--',
                   c='gray', label='dss={:.5f}'.format(dss))
    axs[0, 0].legend(loc=1)  # upper right
    axs[0, 0].set_xlabel('Time(sec)')
    axs[0, 0].set_ylabel('y(m)')

    # impulse response, comparison with an ideal mass/spring/damper
    H_model = co.c2d(co.tf([1], [m, b, k]), dT)
    T_imp = np.arange(0, 10, 0.008)
    _, y_imp = co.impulse_response(H, T_imp)
    _, y_imp_model = co.impulse_response(H_model, T_imp)

    axs[1, 0].plot(T_imp, y_imp[0, :], label='H[n]')
    axs[1, 0].plot(T_imp, y_imp_model[0, :], label='ref-model')
    axs[1, 0].legend()
    axs[1, 0].text(5, 0.0003, 'model(m={:},b={:},k={:})'.format(
        m, b, k), horizontalalignment='center')

    # frequency responses
    mag_yn = 20 * np.log10(mag[0, 0])
    mag_T = 20 * np.log10(mag[1, 0])
    # bounds on H_yn and H_T
    if internal_data is not None:
        T = internal_data['internal responses'][0].shape[0]
    else:
        T = 256
    omegas = np.arange(int(T / 2)) * 2 * np.pi / T / 0.008
    omegas[0] = 1e-2
    wS_inv = np.ones(T) * 100  # infinity
    wS_inv[:int(T / 2)] = np.power(10, np.interp(np.log10(omegas),
                                                 np.log10(freqs_bnd_yn), mag_bnd_yn) / 20)
    mag_wS = 20 * np.log10(wS_inv)
    axs[0, 1].scatter(omegas[:int(T / 2)],
                      mag_wS[:int(T / 2)], label='1/wN', c='C2')
    axs[0, 1].plot(freqs, mag_yn, label='H_yn(z)', c='C0')
    axs[0, 1].plot(freqs, mag_T, label='T(z)', c='C1')
    axs[0, 1].plot([w_nyquist, w_nyquist], [
                   np.min(mag_yn), np.max(mag_yn)], '--', c='red')
    axs[0, 1].plot(freqs_bnd_yn, mag_bnd_yn, 'x--', c='C0', label='1/wN')
    axs[0, 1].plot(freqs_bnd_T, mag_bnd_T, 'x--', c='C1', label='1/wT')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_ylabel('Mag(dB)')
    axs[0, 1].set_xlabel('Freq(rad/s)')
    axs[0, 1].legend()

    # nyquist plot of H_yn (noise to output)
    H_yn = mag[0, 0] * np.exp(1j * phase[0, 0])
    axs[1, 1].plot(H_yn.real, H_yn.imag, 'o-')
    axs[1, 1].set_title("Nyquist plot of H_yn(s)")

    for i in [0, 30, 40, 50, 60, 63, 65, 69, 73, 77, 80, 99]:
        axs[1, 1].text(H_yn[i].real, H_yn[i].imag,
                       "{:.3f} rad/s".format(freqs[i]))

    fig.suptitle('Analysis plots: {:}'.format(controller_name))
    plt.show()

    if internal_data is not None:
        fig, axs = plt.subplots(2, 1)
        (Rval, Nval, Mval, Lval) = internal_data['internal responses']
        T = Rval.shape[0]
        T_Half = int(T / 2)
        Rdft = np.fft.fft(Rval[:, 0, 0], axis=0)
        axs[0].vlines(np.arange(T_Half) * 2 * np.pi /
                      T, 0, np.abs(Rdft[:T_Half]))
        axs[0].set_title('DFT{{R[0, 0]}} T={:d}'.format(T))
        axs[1].plot(internal_data['L'].flatten(), label='L[n]')
        axs[1].plot(internal_data['MB2'].flatten(), label='MB2[n]')
        axs[1].legend()
        plt.show()


def const_jerk_input(
        duration=10,
        accel_duration=0.6,
        f_max=4,
        f_ss=1.0,
        dT=dT):
    """
    """
    t_arr = np.arange(0, duration, dT)
    f_arr = []
    for t in t_arr:
        if t < accel_duration / 2:
            f_arr.append(t * f_max / (accel_duration / 2))
        elif t < accel_duration:
            f_arr.append(f_max - (t - accel_duration / 2) *
                         (f_max - f_ss) / (accel_duration / 2))
        else:
            f_arr.append(f_ss)
    return t_arr, np.array(f_arr)


def SLS_synthesis_p1(Pssd, T, regularization=-1, test_signal='step',
                     Mp=1.01, m=0.5, b=10, k=80,
                     freqs_bnd_yn=[1e-2, 255], mag_bnd_yn=[-10, -10],
                     freqs_bnd_T=[1e-2, 357], mag_bnd_T=[6, 6], T_delay=11):
    """Synthesize a controller using SLS.

    Procedure p1

        Constraints
        - 20c, 20a, 20b (achievability constraints),
        - lower bound on impulse response to prevent negative response,
        - steady-state displacement given constant acting force;
        - noise attenuation: upper bound on the mapping from noise to displacement;
        - robust stability: upper bound on the complementary sensitivity transfer function;

        Objectives
        - regularization using scaled l1 norm on the impulse responses L, MB2
        - distance to a desired mass/spring/damper model (m,b,k)
    """
    print("-- Starting SLS_synthesis_p1")
    # parameters
    nu = 1  # 1 dof controller
    ny = 1

    # form response matrices
    R, N, M, L, H, constraints = Ss.SLS.form_SLS_response_matrices(
        Pssd, nu, ny, T)

    # # constants
    nx = Pssd.states
    A, B1, B2, C1, C2, D11, D12, D21, D22 = Ss.get_partitioned_mats(
        Pssd, nu, ny)

    # select component of H that correspond to the mapping from acting
    # force (and noise) to actual robot displacement (H_yn) and the
    # mapping that corresponds to the complementary transfer function T.
    H_yn = H[0::2, :]
    H_T = H[1::2, :]

    # objective: match the impulse response of a given system
    sys_model = co.c2d(co.tf([1], [m, b, k]), dT)
    _, imp_model = co.impulse_response(sys_model, np.arange(T) * dT)

    # NOTE: have to use norm, sum_of_squares does not work. The reason
    # is the magnitude of the objective function must not be too
    # small. The optimizer seems to get confuse and simply stop
    # working.
    imp_diff = H_yn[T_delay:, 0] - imp_model[0, :T - T_delay]
    weight = np.diag(1 + 2 * (1.0 / T) * np.arange(T - T_delay))
    objective = 1e6 * cvx.norm(weight * imp_diff)

    # try some regularization
    if regularization > 0:
        reg = regularization * (cvx.norm1(H_yn))
    else:
        reg = cvx.abs(cvx.Variable())

    # constraint in frequency-domain, if specified
    W_dft = Ss.dft_matrix(T)
    Hz_yn = W_dft * H_yn
    Hz_T = W_dft * H_T
    omegas = np.arange(int(T / 2)) * 2 * np.pi / T / dT
    omegas[0] = 1e-2

    # upper bound for noise attenuation
    wN_inv = np.ones((T, 1)) * 100  # infinity
    wN_inv[:int(T / 2), 0] = np.power(10, np.interp(np.log10(omegas),
                                                    np.log10(freqs_bnd_yn), mag_bnd_yn) / 20)

    # upper bound for complementary sensitivity transfer function
    wT_inv = np.ones((T, 1)) * 100  # infinity
    wT_inv[:int(T / 2), 0] = np.power(
        10, np.interp(np.log10(omegas), np.log10(freqs_bnd_T), mag_bnd_T) / 20)

    # add both frequency-domian constraints
    constraints.append(cvx.abs(Hz_yn) <= wN_inv)
    constraints.append(cvx.abs(Hz_T) <= wT_inv)

    # optimize
    obj = cvx.Minimize(objective + reg)
    prob = cvx.Problem(obj, constraints)
    print("-- [SLS_synthesis_p1] Preparing problem with cvxpy!")
    prob.solve(verbose=True, solver=cvx.MOSEK)
    print("-- [SLS_synthesis_p1] optimization status: {:}".format(prob.status))
    print(
        "-- [SLS_synthesis_p1] obj = {:}, reg = {:}".format(objective.value, reg.value))

    if prob.status != "optimal":
        return None, None

    print("-- [SLS_synthesis_p1] Forming controllers!")
    # form controllers (Structure 1, Figure 4b, Wang 2018)
    L_value = np.array(L.value).reshape(ny, nu, -1)
    MB2_value = np.array((M * B2).value).reshape(nu, nu, -1)

    # since ny=nu=1, we have
    fir_den = [1] + [0 for n in range(T - 1)]
    MB2_tf = co.tf(MB2_value[0, 0], fir_den, dT)
    L_tf = co.tf(L_value[0, 0], fir_den, dT)
    K = co.feedback(1, MB2_tf, sign=-1) * L_tf
    K = Ss.tf2ss(K, minreal=True)

    # response mapping
    Rval = np.array(
        [R[n * nx: (n + 1) * nx, :].value for n in range(T)]).reshape(-1, nx, nx)
    Nval = np.array(
        [N[n * nx: (n + 1) * nx, :].value for n in range(T)]).reshape(-1, nx, ny)
    Mval = np.array(
        [M[n * nu: (n + 1) * nu, :].value for n in range(T)]).reshape(-1, nu, nx)
    Lval = np.array(
        [L[n * nu: (n + 1) * nu, :].value for n in range(T)]).reshape(-1, nu, ny)
    Hval_yn = np.array(H_yn.value)
    Hval_T = np.array(H_T.value)

    return K, {'internal responses': (Rval, Nval, Mval, Lval),
               'output impulse': (Hval_yn, Hval_T),
               'L': L_value, 'MB2': MB2_value}


def main():
    Ptf_design = plant(Hdelay=0.05, Hgain=50)
    Pss_design = Ss.tf2ss(Ptf_design, minreal=True)
    Pssd_design = co.c2d(Pss_design, dT)

    # synthesize controller
    freqs_bnd_T = [1e-2, 2.3, 7.3, 25, 357]
    mag_bnd_T = [-3, -3, -3, -10, -40]
    # freqs_bnd_yn = [1e-2, 3.0, 30, 80, 255]  # rad
    # mag_bnd_yn = [-10, -10, -20, -74, -100]  # db
    freqs_bnd_yn = [1e-2, 3.0, 20, 50, 255]  # rad
    mag_bnd_yn = [-20, -20, -20, -94, -130]  # db
    Asls, internal_data = SLS_synthesis_p1(Pssd_design, 256, regularization=1,
                                           freqs_bnd_T=freqs_bnd_T, mag_bnd_T=mag_bnd_T,
                                           freqs_bnd_yn=freqs_bnd_yn, mag_bnd_yn=mag_bnd_yn,
                                           m=1.5, b=24, k=60, T_delay=7)

    # # test/analysis
    Ptf_test = plant(Hgain=50, Hdelay=0.05)
    Pssd_test = co.c2d(Ss.tf2ss(Ptf_test, minreal=True), dT)
    if Asls is not None:
        analysis(Pssd_test, Asls,
                 internal_data=internal_data, Tr=1.0, controller_name='SLS',
                 freqs_bnd_T=freqs_bnd_T, mag_bnd_T=mag_bnd_T, freqs_bnd_yn=freqs_bnd_yn,
                 mag_bnd_yn=mag_bnd_yn, m=1.5, b=24, k=60)

    analysis(Pssd_test, OldControllers.A1, Tr=1.0, controller_name='admittance',
             freqs_bnd_T=freqs_bnd_T, mag_bnd_T=mag_bnd_T,
             freqs_bnd_yn=freqs_bnd_yn, mag_bnd_yn=mag_bnd_yn,
             m=1.5, b=28, k=65)
    analysis(Pssd_test, OldControllers.Ac14, Tr=1.0, controller_name='Hinf',
             freqs_bnd_T=freqs_bnd_T, mag_bnd_T=mag_bnd_T,
             freqs_bnd_yn=freqs_bnd_yn, mag_bnd_yn=mag_bnd_yn,
             m=1.5, b=28, k=65)

    # to print controller
    Ss.SLS.print_controller(internal_data['L'], internal_data['MB2'])

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    main()
