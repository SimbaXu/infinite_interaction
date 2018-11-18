# previous synthesis script
from Oct23_1dof_admittance import plant, analysis
# standard libs
import SLSsyn as Ss
import numpy as np
import control as co
import matplotlib.pyplot as plt
import cvxpy as cvx


def load_atomic_pool(Pssd_design, T=5):
    """A pool of basis atomic controller.

    Will be used for finding an optimal atomic controller via affine
    combinations.

    Args:
        Pssd_design(tf.ss): A discrete-time state-space system.
    """
    s = co.tf([1, 0], [1])
    ctrl_pool = [
        # co.tf([1], [0.1, 1, 1]),  # this controller should be unstable
        1 / (s**2 + 4 * s + 6),
        1 / (0.2 * s**2 + 3 * s + 2),
        1 / (2 * s**2 + 5 * s + 10),
        1 / (0.2 * s**2 + 4 * s + 6),
        1 / (0.1 * s**2 + 6 * s + 6),
        co.tf([0], [1]),  # zero controller
        5e-3,
        5e-3 * (1 + s) / (1 + 1e-3 * s),
        5e-3 * (1 + s) / (1 - s),
        # co.tf([0.001], [1]),
        # co.tf([0.01], [1, 1])
    ]
    report_string = "Controllers pool contain: {:d} controllers\n".format(
        len(ctrl_pool))
    for i, c in enumerate(ctrl_pool):
        report_string += "> [c{:d}] {:}".format(i, repr(c))

    # Now, generate the responses for each controller in the pool. Those that are not
    # stabilizing are excluded.
    atomic_pool = {}
    for i, ctrl in enumerate(ctrl_pool):
        atomic = AtomicFIRController.init_via_simulation(
            Pssd_design, proc_ctrl(ctrl), T=T)
        if atomic is not None:
            atomic.name = "c{:d}".format(i)
            atomic_pool[atomic.name] = atomic
    report_string = "-- {:d}/{:d} controllers in controllers pool are stabilizing".format(
        len(atomic_pool), len(ctrl_pool))
    print(report_string)
    return atomic_pool


def main():
    # constants
    dT = 0.008

    # setup
    Ptf_design = plant(Hdelay=0.05, Hgain=50)
    Pss_design = Ss.tf2ss(Ptf_design, minreal=True)
    Pssd_design = co.c2d(Pss_design, dT)
    atomic_pool = load_atomic_pool(Pssd_design, T=10)

    # frequencies bounds
    freqs_bnd_T = [1e-2, 2.3, 7.3, 25, 357]
    mag_bnd_T = [-3, -3, -3, -10, -40]
    freqs_bnd_yn = [1e-2, 3.0, 20, 50, 255]  # rad
    mag_bnd_yn = [-20, -20, -20, -94, -130]  # db
    m = 1.5
    b = 28
    k = 65
    atomic_sls, internal_data = SLS_synthesis2(
        atomic_pool, m, b, k, freqs_bnd_yn, mag_bnd_yn, freqs_bnd_T, mag_bnd_T)
    # just for fun, analyze a controller
    freqs_bnd_T = [1e-2, 2.3, 7.3, 25, 357]
    mag_bnd_T = [-3, -3, -3, -10, -40]
    # freqs_bnd_yn = [1e-2, 3.0, 30, 80, 255]  # rad
    # mag_bnd_yn = [-10, -10, -20, -74, -100]  # db
    freqs_bnd_yn = [1e-2, 3.0, 20, 50, 255]  # rad
    mag_bnd_yn = [-20, -20, -20, -94, -130]  # db
    analysis(
        Pssd_design,
        # atomic_pool['c4'].get_executable_controller(),
        atomic_sls.get_executable_controller(),
        Tr=1.0,
        controller_name='admittance',
        freqs_bnd_T=freqs_bnd_T,
        mag_bnd_T=mag_bnd_T,
        freqs_bnd_yn=freqs_bnd_yn,
        mag_bnd_yn=mag_bnd_yn,
        m=1.5,
        b=28,
        k=65
    )
    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


class AtomicFIRController(object):
    """A atomic controller.

    An atomic controller can be generated by simulation of a plant and
    controller, or as an affine combination of several atomic
    controllers.

    """
    nb_atoms = 0

    def __init__(self):
        # do nothing
        self.name = "atom nb.{:d}".format(AtomicFIRController.nb_atoms)
        AtomicFIRController.nb_atoms += 1

    @staticmethod
    def init_as_convex_comb(weights, controllers):
        self = AtomicFIRController()
        self.name = "nb.{:d}(affcomb)".format(AtomicFIRController.nb_atoms - 1)
        self.exec_ctrl = None
        assert len(weights) == len(controllers)
        # normalize
        self.weights = np.abs(weights) / np.sum(np.abs(weights))

        # basic parameters
        self.T = controllers[0].T
        self.dT = controllers[0].dT
        self.NT = controllers[0].NT
        self.nx = controllers[0].nx
        self.ny = controllers[0].ny
        self.nxc = controllers[0].nxc
        self.nu = controllers[0].nu
        self.ny = controllers[0].ny

        # responses
        self.R = 0
        self.N = 0
        self.M = 0
        self.L = 0
        self.MB2 = 0
        self.H = 0
        for wi, ctrli in zip(self.weights, controllers):
            self.R += wi * ctrli.R
            self.N += wi * ctrli.N
            self.M += wi * ctrli.M
            self.L += wi * ctrli.L
            self.MB2 += wi * ctrli.MB2
            self.H += wi * ctrli.H
        return self

    @staticmethod
    def init_via_simulation(Pssd_design, A1dss, T=5):
        self = AtomicFIRController()
        nu = A1dss.outputs
        ny = A1dss.inputs
        nxc = A1dss.states
        nx = Pssd_design.states
        dT = Pssd_design.dt
        self.T = T
        self.dT = dT
        self.NT = int(T / dT)
        self.nx = nx
        self.ny = ny
        self.nxc = nxc
        self.nu = nu
        self.ny = ny
        self.nw = Pssd_design.inputs - A1dss.outputs
        self.nz = Pssd_design.outputs - A1dss.inputs

        # closed-loop response, check stability
        Pssd_cl = Pssd_design.lft(A1dss)
        w, q = np.linalg.eig(Pssd_cl.A)
        wmax = w[np.argmax(np.abs(w))]
        if np.abs(wmax) > 1:
            print(" -- Closed-loop system UNSTABLE. w_max={:}"
                  " ----> Return None".format(wmax))
            return None
        else:
            print(" -- Closed-loop system STABLE. w_max={:}".format(wmax))

        # control systems
        self.plant_ol = Pssd_design
        self.plant_cl = Pssd_cl
        self.exec_ctrl = A1dss

        # matrix impulse response
        A, B1, B2, C1, C2, D11, D12, D21, D22 = map(
            np.array, Ss.get_partitioned_mats(Pssd_design, nu, ny))
        Ak, Bk, Ck, Dk = A1dss.A, A1dss.B, A1dss.C, A1dss.D
        assert np.all(D22 == 0)

        # Concatenate system matrices to find the matrices of the ss
        # rep of the impulse response mapping
        # [R N] =   A_cb  |  B_cb
        # [M L]     ------|-------
        #           C_cb  |  D_cb
        A_cb = np.block([[A + B2.dot(Dk).dot(C2), B2.dot(Ck)],
                         [Bk.dot(C2), Ak]])
        B_cb = np.block([[np.eye(nx), B2.dot(Dk)],
                         [np.zeros((nxc, nx)), Bk]])
        C_cb = np.block([[np.eye(nx), np.zeros((nx, nxc))],
                         [Dk.dot(C2), Ck]])
        D_cb = np.block([[np.zeros((nx, nx)), np.zeros((nx, ny))],
                         [np.zeros((nu, nx)), Dk]])
        P_resp_dss = co.ss(A_cb, B_cb, C_cb, D_cb, dT)

        # Compute impulse responses of the matrices R, N, M, L
        Tarr = np.arange(0, T, dT)  # 5 sec horizon
        NT = Tarr.shape[0]
        mtrx_impulses = []
        for i in range(self.nx + self.ny):
            _, yarr = co.impulse_response(
                P_resp_dss, Tarr, input=i, transpose=True)
            if np.linalg.norm(yarr[-1]) > 1e-10:
                raise(ValueError("Matrix Impulse response has non-zero norm at the last time-step! ({:}) "
                                 "Consider increasing horizont T. (Current val={:})".format(
                                     np.linalg.norm(yarr[-1]), T)))
            mtrx_impulses.append(yarr)
        mtrx_impulses = np.stack(mtrx_impulses, axis=2)

        # Individual responses
        self.R = mtrx_impulses[:, :self.nx, :self.nx]
        self.N = mtrx_impulses[:, :self.nx, self.nx:self.nx + self.ny]
        self.M = mtrx_impulses[:, self.nx:self.nx + self.nu, :self.nx]
        self.L = mtrx_impulses[:, self.nx:self.nx + self.nu, self.nx:self.nx + self.ny]
        self.MB2 = self.M@B2

        # Output Impulse responses H (from control input w to output z)
        output_impulses = []
        for i in range(self.nw):
            _, yarr = co.impulse_response(
                Pssd_cl, Tarr, input=i, transpose=True
            )
            # no need to check for finite horizon since H is a linear
            # combination of the matrix responses
            output_impulses.append(yarr)
        self.H = np.stack(output_impulses, axis=2)

        # frequency response of the mapping
        W = Ss.dft_matrix(self.NT)
        self.H_fft = np.zeros_like(self.H) * 1j
        for i in range(self.nw):
            for j in range(self.nz):
                self.H_fft[:, i, j] = W@self.H[:, i, j]
        return self

    def get_executable_controller(self):
        if self.exec_ctrl:
            return self.exec_ctrl
        else:
            fir_den = [1] + [0 for n in range(self.NT - 1)]
            MB2_tf = co.tf(self.MB2[:, 0, 0], fir_den, self.dT)
            L_tf = co.tf(self.L[:, 0, 0], fir_den, self.dT)
            A1dss_FIR = co.feedback(1, MB2_tf, sign=-1) * L_tf
            A1dss_FIR = Ss.mtf2ss(A1dss_FIR, minreal=False)
            self.exec_ctrl = A1dss_FIR
            return self.exec_ctrl


def SLS_synthesis2(atoms, m, b, k, freqs_bnd_yn=[1e-2, 255], mag_bnd_yn=[-10, -10],
                   freqs_bnd_T=[1e-2, 357], mag_bnd_T=[6, 6]):
    """ The new (version 2) SLS synthesis procedure.
    """
    keys = list(atoms.keys())
    key0 = keys[0]
    NT = atoms[key0].NT
    dT = atoms[key0].dT

    # optimization parameters
    theta = cvx.Variable(len(atoms))
    regularization = cvx.norm1(theta)
    constraints = [cvx.sum(theta) == 1]

    # desired impulse response
    Tarr = np.arange(NT) * dT
    desired_model = co.c2d(co.tf([1], [m, b, k]), dT)
    _, H_desired = co.impulse_response(desired_model, Tarr, transpose=True)
    H_desired = np.array(H_desired).flatten()

    # affine controller's properties
    H_affcomb = 0
    H_fft_affcomb = 0
    for i, key in enumerate(keys):
        H_affcomb += theta[i] * atoms[key].H[:, 0]
        H_fft_affcomb += theta[i] * atoms[key].H_fft

    # list of frequencies for frequency-domain constraints
    omegas = np.arange(int(NT / 2)) * 2 * np.pi / NT / dT
    omegas[0] = 1e-2

    # frequency-based upper bound for noise attenuation
    wN_inv = np.ones(NT) * 100  # infinity
    wN_inv[:int(NT / 2)] = np.power(
        10, np.interp(np.log10(omegas), np.log10(freqs_bnd_yn), mag_bnd_yn) / 20)

    # frequency-based upper bound on the complementary sensitivity
    # function for robust stability
    wT_inv = np.ones(NT) * 100  # infinity
    wT_inv[:int(NT / 2)] = np.power(
        10, np.interp(np.log10(omegas), np.log10(freqs_bnd_T), mag_bnd_T) / 20)

    # add frequency-based constraints
    constraints.append(cvx.abs(H_fft_affcomb[:, 0]) <= wN_inv)
    constraints.append(cvx.abs(H_fft_affcomb[:, 1]) <= wT_inv)

    # mimic desired system response
    mse = cvx.norm(H_affcomb - H_desired)
    obj = cvx.Minimize(mse + regularization)
    constraints.append(
        H_affcomb >= -1e-4
    )

    # opt problem
    prob = cvx.Problem(obj, constraints)
    prob.solve(verbose=True)

    if prob.status == "optimal":
        report_string = " -- Optimization successful!\n     controller key | weight\n"
        for i, key in enumerate(keys):
            report_string += "    {:10} | {:.3f}\n".format(key, theta[i].value)
        print(report_string)
        # form the affinely combined controller
        atoms_list = [atoms[key] for key in keys]
        res_ctrl = AtomicFIRController.init_as_convex_comb(theta.value, atoms_list)
        return res_ctrl, None
    else:
        print(" -- Optimization fails! Returning None!")
        return None, None


def proc_ctrl(A1c, dT=0.008, use_matlab=True):
    if isinstance(A1c, float) or isinstance(A1c, int):
        A1c = co.tf([1], [1]) * A1c
    A1d = co.c2d(A1c, dT)
    A1dss = Ss.tf2ss(A1d, minreal=True, via_matlab=use_matlab)
    return A1dss


if __name__ == '__main__':
    main()
