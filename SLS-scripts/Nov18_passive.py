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


def plant(dT=0.008):
    """ Discrete-time transfer function.

    Diagram:
                          +--------+    v
                   +----->|  R(s)  |-------->
                   |      +--------+
                   |      +--------+
                   +<-----|  K(s)  |<------- f
                       u  +--------+ m


    Transfer functions:

    [v] = [0 R(s)] [f]
    [m]   [1 0   ] [u]

    Signals:
        v: output velocity
        m: measured force
        f: acting force
        u: position command
    """
    s = co.tf([1, 0], [1])
    # NOTE: This mapping is from commanding position to velocity
    R = (-s + 55.56) / (s + 55.56) / (0.0437 * s + 1) * s
    P = Ss.tf_blocks([[0, R], [1, 0]])
    Pss = Ss.tf2ss(P, minreal=True)
    Pdss = co.c2d(Pss, dT)
    return Pdss


def SLS_synthesis_passive(Pdss, T, m=1.5, b=24, k=60, T_delay=7,
                          freqs_bnd=None, mags_bnd=None, dT=0.008):
    """
    """
    # parameters
    nu = 1  # 1 dof controller
    ny = 1

    # form response matrices
    # R, N, M, L, H, constraints = Ss.SLS.form_SLS_response_matrices(
        # Pdss, nu, ny, T)
    nx = Pdss.states
    A, B1, B2, C1, C2, D11, D12, D21, D22 = Ss.get_partitioned_mats(
        Pdss, nu, ny)

    ny_out, nu_exo = D11.shape  # control output dimension, also known as nz, nw

    # variables
    R = cvx.Variable((T * nx, nx))
    N = cvx.Variable((T * nx, ny))
    M = cvx.Variable((T * nu, nx))
    L = cvx.Variable((T * nu, ny))

    # constraints
    # 20c
    constraints = [
        R[:nx, :] == 0, M[:nu, :] == 0, N[:nx, :] == 0,
    ]
    # 20a: t20a_1 R - t20a_2 R - t20a_3 M = t20a_4
    t20a_1 = np.zeros((T * nx, T * nx))
    t20a_2 = np.zeros((T * nx, T * nx))
    t20a_3 = np.zeros((T * nx, T * nu))
    t20a_4 = np.zeros((T * nx, nx))
    for n in range(T):
        if n != T - 1:
            t20a_1[n * nx: (n + 1) * nx,
                   (n + 1) * nx: (n + 2) * nx] = np.eye(nx)
        t20a_2[n * nx: (n + 1) * nx, n * nx: (n + 1) * nx] = A
        t20a_3[n * nx: (n + 1) * nx, n * nu: (n + 1) * nu] = B2
        if n == 0:
            t20a_4[:nx, :nx] = np.eye(nx)
    constraints.extend(
        [t20a_1 * R - t20a_2 * R - t20a_3 * M == t20a_4,
         t20a_1 * N - t20a_2 * N - t20a_3 * L == 0]
    )

    # 20b: t20a_1 R - R * A - N C2 == t20a_4
    # 20b-lower: t20b_1 M - M * A - L * C2 == 0
    t20b_1 = np.zeros((T * nu, T * nu))
    for n in range(T):
        if n != T - 1:
            t20b_1[n * nu: n * nu + nu,
                   (n + 1) * nu: (n + 2) * nu] = np.eye(nu)

    constraints.extend(
        [
            t20a_1 * R - R * A - N * C2 == t20a_4,
            t20b_1 * M - M * A - L * C2 == 0
        ]
    )

    # mapping from exo input to control output
    from scipy.linalg import block_diag
    C1_blk = block_diag(* [C1] * T)
    D12_blk = block_diag(* [D12] * T)
    D11_blk = np.zeros((T * ny_out, nu_exo))
    D11_blk[:ny_out, :nu_exo] = D11
    H = C1_blk * R * B1 + D12_blk * M * B1 + \
        C1_blk * N * D21 + D12_blk * L * D21 + D11_blk

    # fourier transform of H
    W_dft = Ss.dft_matrix(T)
    H_dft = W_dft @ H

    # constants
    nx = Pdss.states
    A, B1, B2, C1, C2, D11, D12, D21, D22 = Ss.get_partitioned_mats(
        Pdss, nu, ny)

    # objective: match the impulse response of a given system
    # min ||H - resp_desired||
    sys_model = co.c2d(co.tf([1, 0], [m, b, k]), dT)
    _, resp_desired = co.impulse_response(sys_model, np.arange(T) * dT, transpose=True)
    resp_diff = H[T_delay:] - resp_desired[:(T - T_delay)]
    weight = np.diag(1 + 2 * (1.0 / T) * np.arange(T - T_delay)) * 1e3  # normalize
    objective = cvx.norm(weight * resp_diff)

    # constraint: passivity
    constraints.append(cvx.real(H_dft) >= 1e-3)

    # constraint: frequency response magnitude
    # list of angular velocities corresponding to the
    # discrete(digital) frequencies 2pi / T * k, k=[0...N-1]
    omegas = np.arange(int(T / 2)) * 2 * np.pi / T / dT
    omegas[0] = 1e-2
    mags = np.ones((int(T / 2), 1))
    mags[:, 0] = np.power(
        10, np.interp(np.log10(omegas), np.log10(freqs_bnd), mags_bnd) / 20)
    # constraints.append(cvx.abs(H_dft[:int(T / 2)]) <= mags)

    prob = cvx.Problem(cvx.Minimize(objective), constraints)
    prob.solve(verbose=True, solver=cvx.SCS)
    import ipdb; ipdb.set_trace()
    a = 1


def main():
    # define constant
    dT = 0.008
    # plant
    Pdss_design = plant()
    # bound on loop gain magnitude
    lg_freqs = [1e-2, 2.3, 7.3, 25, 357]
    lg_mag = [-3, -3, -3, -10, -40]
    # synthesize controller
    Asls, internal_data = SLS_synthesis_passive(
        Pdss_design, 256, freqs_bnd=lg_freqs, mags_bnd=lg_mag)

if __name__ == '__main__':
    main()
