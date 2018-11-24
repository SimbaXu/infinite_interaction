"""This script dessigns a controller using SLS for rendering a simple
robot model with time-delay approximated bya first-order pade
approximation passive. While at the same time, feel good to the
exogeneous agent.

Method: SYN1

status, as of 18/11/2018: Unable to synthesize with MOSEK

"""
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import SLSsyn as Ss
import scipy.sparse as sparse


def main():
    # auxilliary variables
    s = co.tf([1, 0], [1])
    z = co.tf([1, 0], [1], 0.008)
    dT = 0.008

    # plant
    R = s / (1 + 0.0437 * s)
    Rd = co.c2d(R, 0.008) / z**2

    # I parametrize $Q(z)$ as a FIR filter with, say 10 taps. The
    # output mapping $H(z)$ is therefore a linear combination of
    # delayed impulse responses of $R(z)$
    N = 734
    Td = 500  # 100 steps
    # impulse
    Harr = np.zeros((N, Td + N))  # number of tap + number of step each impulse
    Harr_dft = 1j * np.zeros((N, Td + N))
    # compute
    _, H0_ = co.impulse_response(Rd, np.arange(Td) * dT, transpose=True)
    Harr[0, :Td] = H0_[:, 0]
    # other taps
    for i in range(1, N):
        Harr[i, i: Td + i] = H0_[:, 0]
    # dft
    W_dft = Ss.dft_matrix(Td + N)
    for i in range(N):
        Harr_dft[i] = W_dft @ Harr[i]

    # an ideal controller
    A = co.c2d(co.tf([1], [3, 10, 10]), dT)
    H = Rd * A
    _, H_actual = co.impulse_response(H, np.arange(Td + N) * dT, transpose=True)
    _, taps = co.impulse_response(A, np.arange(N) * dT, transpose=True)
    # response as computed using Q-param
    H_Qparam = 0
    for i in range(N):
        H_Qparam += Harr[i, :] * taps[i, 0]

    H_actual = H_actual[:, 0]
    plt.plot(H_actual, label="output")
    plt.plot(taps, label="tap")
    plt.plot(H_Qparam, label="output as linear combination Qparam")
    plt.legend()
    plt.grid()
    plt.show()

    # a single tap
    fig, axs = plt.subplots(1, 2, figsize=[10, 4])
    axs[0].plot(Harr[0], label='tap0')
    axs[0].plot(Harr[5], label='tap5')
    axs[0].set_title("A single impulse response")
    axs[0].legend()
    axs[1].plot(Harr_dft[0].real, Harr_dft[0].imag, 'o-', label='tap0')
    axs[1].plot(Harr_dft[5].real, Harr_dft[5].imag, 'o-', label='tap5')
    axs[1].set_title("Nyquist Plot")
    axs[1].grid()
    axs[1].legend()
    plt.show()

    # the desired impulse response
    sys_model = co.c2d(co.tf([1, 0], [3, 10, 10]), dT)
    _, H_desired = co.impulse_response(sys_model, np.arange(Td + N) * dT, transpose=True)
    H_desired = H_desired[:, 0]
    H_desired[3:] = H_desired[:Td + N - 3]
    H_desired[:3] = 0
    plt.plot(H_desired)
    plt.grid()
    plt.show()

    # optimize with cvxpy
    weights = cvx.Variable(N)
    weight0 = cvx.Variable()
    H = weight0 * 0
    H_dft = weight0 * np.zeros(Td + N) * 1j
    for i in range(N):
        H += weights[i] * Harr[i]
        H_dft += weights[i] * Harr_dft[i]
    obj = cvx.norm(H - H_desired)
    constraints = [cvx.sum(weights) + weight0 == 1,
                   cvx.real(H_dft[:100]) >= 0,
                   cvx.abs(H_dft[100:]) <= 1e-2]
    prob = cvx.Problem(cvx.Minimize(obj), constraints)
    prob.solve(verbose=True)
    if prob.status == 'optimal':
        Qval = np.array(weights.value).flatten()
        plt.plot(H.value, label='actual')
        plt.plot(H_desired, label='desired')
        plt.legend()
        plt.show()

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()
    

if __name__ == '__main__':
    main()
