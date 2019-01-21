import numpy as np
try:
    import matlab.engine
    FOUND_MATLAB = True
except ImportError:
    # unable to find matlab-engine
    FOUND_MATLAB = False
import control as co
import yaml
import os
import cvxpy as cvx
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class MatlabEngine(object):
    """ Matlab Engine (Singleton)

    """
    matlab_engine_instance = None
    loaded_matlab = False

    def __init__(self):
        if MatlabEngine.matlab_engine_instance is None:
            try:
                print("-- Trying to start matlab engine!")
                MatlabEngine.matlab_engine_instance = matlab.engine.start_matlab()
                print("-- Matlab Engine started successfully!")
            except Exception as e:
                print(e)
                MatlabEngine.matlab_engine_instance = None

    def __getattr__(self, name):
        """ Delegate calls to MatlabEngine to the internal instance.
        """
        return getattr(self.matlab_engine_instance, name)


def impulse2tf(M, dT):
    """ Convert a impulse reponse to a tranfer matrix.

    Args:
        M (ndarray): Shape (T, ny, nu). Impulse reponse of the filter.
        dT (float): Time step.
    """
    M = np.array(M)
    T, ny, nu = M.shape
    common_den = [1] + [0 for i in range(T - 1)]

    fir_num = np.zeros((ny, nu, T))
    fir_den = np.zeros((ny, nu, T))
    for i in range(ny):
        for j in range(nu):
            fir_den[i, j, :] = common_den
            fir_num[i, j, :] = M[:, i, j]
    return co.tf(fir_num, fir_den, dT)


def get_partitioned_mats(P, nu, ny):
    """ Return the partitioned matrices of P.

    Args:
        P (control.ss): A state-space transfer function.

    Returns:
        ([np.matrix]): A, B1, B2, C1, C2, D11, D12, D21, D22

    """
    if not isinstance(P, co.StateSpace):
        raise ValueError('P needs to be a StateSpace system')
    A = P.A
    B1 = P.B[:, :P.inputs - nu]
    B2 = P.B[:, P.inputs - nu:]
    C1 = P.C[:P.outputs - ny, :]
    C2 = P.C[P.outputs - ny:, :]
    D11 = P.D[:P.outputs - ny, :P.inputs - nu]
    D12 = P.D[:P.outputs - ny, P.inputs - nu:]
    D21 = P.D[P.outputs - ny:, :P.inputs - nu]
    D22 = P.D[P.outputs - ny:, P.inputs - nu:]

    return A, B1, B2, C1, C2, D11, D12, D21, D22


def dft_matrix(N):
    """ Return W s.t. DFT{x} = W x
    """
    W0 = np.exp(- 2j * np.pi / N)
    W = np.ones((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            W[i, j] = W0 ** (i * j)
    return W


def get_num_dens(Plist, dt=None):
    """Extract the numerators and denominators of a list of SISO TFs.

    The list of lists of SISO TFs is assumed to have the following
    form:

    [[P11, P12, ...],
     [P21, P22, ...],
     [P21, P22, ...],
     ...]

    The tfs are also summed to have consistent sampling time.
    """
    Nout = len(Plist)
    Nin = len(Plist[0])
    nums = []
    dens = []
    for i in range(Nout):
        nums.append([])
        dens.append([])
        for j in range(Nin):
            if isinstance(Plist[i][j], co.xferfcn.TransferFunction):
                nums[i].append(Plist[i][j].num[0][0])
                dens[i].append(Plist[i][j].den[0][0])
                if co.isdtime(Plist[i][j]):
                    if dt is None:
                        dt = Plist[i][j].dt
                    elif dt != Plist[i][j].dt:
                        raise(ValueError("Inconsistent sampling time"))
                    else:  # okay
                        pass

            # other than tf function, must be a float
            elif isinstance(Plist[i][j], float) or isinstance(Plist[i][j], int):
                nums[i].append([float(Plist[i][j])])
                dens[i].append([1])
            else:
                raise(ValueError(
                    "Input: {:} not understood. Must be a transfer function or a numeric value."))

    return nums, dens, dt


def tf_blocks(Plist, dt=None):
    """Create a MIMO TF from multiple SISO TFs.

    See :func:`get_num_dens` for details.
    """
    nums, dens, dt = get_num_dens(Plist, dt=dt)

    if dt is None:  # continuous time
        return co.tf(nums, dens)
    else:  # discrete time
        return co.tf(nums, dens, dt)


def lft(P, K, nu=-1, ny=-1):
    """ Linear Fractional Transformation

    https://www.mathworks.com/help/control/ref/lft.html

    """
    return co.ss(P).lft(K, nu, ny)


def mminreal(P):
    """ Matlab Engine based minreal.
    """
    raise NotImplementedError


def tf2ss(P, minreal=False, via_matlab=True):
    if via_matlab:
        return mtf2ss(P, minreal)
    else:
        Pss = co.tf2ss(P)
        Pss = Pss.minreal()
        return Pss


def mtf2ss(P, minreal=False):
    """Convert a Transfer Function to State-Space.

    The conversion is done using Matlab via Matlab Engine.

    If `minreal` is true, also perform a minimal realization.

    """
    eng = MatlabEngine()
    maxlen = 0
    num = []
    den = []
    num_cmd = '{'
    den_cmd = '{'
    for i in range(P.outputs):
        for j in range(P.inputs):
            eng.workspace['num{:d}{:d}'.format(
                i, j)] = matlab.double(P.num[i][j].tolist())
            eng.workspace['den{:d}{:d}'.format(
                i, j)] = matlab.double(P.den[i][j].tolist())
            num_cmd += 'num{:d}{:d}, '.format(i, j)
            den_cmd += 'den{:d}{:d}, '.format(i, j)
        num_cmd += ';'
        den_cmd += ';'
    num_cmd += '}'
    den_cmd += '}'

    eng.eval('num=' + num_cmd + ';', nargout=0)
    eng.eval('den=' + den_cmd + ';', nargout=0)
    if co.isctime(P):
        eng.eval('P=tf(num, den);', nargout=0)
    else:
        eng.eval('P=tf(num, den, {:f});'.format(P.dt), nargout=0)
    eng.eval('Pss=ss(P);', nargout=0)
    if minreal:
        eng.eval('Pss=minreal(Pss);', nargout=0)
    A = np.array(eng.eval('Pss.A'))
    B = np.array(eng.eval('Pss.B'))
    C = np.array(eng.eval('Pss.C'))
    D = np.array(eng.eval('Pss.D'))
    if co.isctime(P):
        return co.ss(A, B, C, D)
    else:
        return co.ss(A, B, C, D, P.dt)


class SLS:
    """ A collection of functions used in SLS synthesis.
    """

    def print_controller(L, MB2, file_name='super.yaml',
                         controller_name='super'):
        """ Print the coefficients to yaml format.
        """
        CONFIG_DIR = "/home/hung/catkin_ws/src/infinite_interaction/config"
        file_dir = os.path.join(CONFIG_DIR, file_name)
        # convert to simple python float before dumping
        L = list([float(L_) for L_ in np.array(L).flatten()])
        MB2 = list([float(L_) for L_ in np.array(MB2).flatten()])
        ny = 1
        nu = 1
        T = len(L)
        ctrl_dict = {'T': T, 'ny': ny, 'nu': nu, 'L': L, 'MB2': MB2}
        output_dict = {
            'fir_siso': {
                controller_name: ctrl_dict
            }
        }
        with open(file_dir, 'w') as f:
            yaml.dump(output_dict, f, default_flow_style=False)
        print("-- Wrote controller to {:}".format(file_dir))

    def form_SLS_response_matrices(Pdss, nu, ny, T):
        """Form response matrices and output matrices.

        Return matrices {R, N, M, L} respectively and output matrices
        {H}. Also return all achievability constraints that guarantee
        the mappings are valid.

        Return:
            R, N, M, L (cvxpy.Variables): Closed-loop reponse matrices.
            H (cvxpy.Variable): Closed-loop output from exogeneous input to exogeneous output.
            constraints (list of cvxpy.Constraint): Achievability constraints.

        """
        nx = Pdss.states
        A, B1, B2, C1, C2, D11, D12, D21, D22 = get_partitioned_mats(
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
        t20a_1 = sparse.lil_matrix((T * nx, T * nx))
        t20a_2 = sparse.lil_matrix((T * nx, T * nx))
        t20a_3 = sparse.lil_matrix((T * nx, T * nu))
        t20a_4 = sparse.lil_matrix((T * nx, nx))
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
        t20b_1 = sparse.lil_matrix((T * nu, T * nu))
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
        C1_blk = sparse.block_diag([C1] * T)
        D12_blk = sparse.block_diag([D12] * T)
        D11_blk = sparse.lil_matrix((T * ny_out, nu_exo))
        D11_blk[:ny_out, :nu_exo] = D11
        H = C1_blk * R * B1 + D12_blk * M * B1 + \
            C1_blk * N * D21 + D12_blk * L * D21 + D11_blk

        # debug
        ctrl = co.c2d(tf2ss(co.tf([5e-3], [1]), minreal=True), 0.008)
        Rnp, Nnp, Mnp, Lnp = SLS.compute_cl_matrix_responses(Pdss, ctrl, 0.008 * 256)
        Rnp = Rnp.reshape(*R.shape)
        Nnp = Nnp.reshape(*N.shape)
        Mnp = Mnp.reshape(*M.shape)
        Lnp = Lnp.reshape(*L.shape)
        # check constraint

        tmp1 = np.matmul(t20a_1, Rnp) - np.matmul(t20a_2, Rnp) - np.matmul(t20a_3, Mnp)
        tmp2 = np.matmul(t20a_1, Nnp) - np.matmul(t20a_2, Nnp) - np.matmul(t20a_3, Lnp)
        tmp3 = np.matmul(t20a_1, Rnp) - np.matmul(Rnp, A) - np.matmul(Nnp, C2)
        tmp4 = np.matmul(t20b_1, Mnp) - np.matmul(Mnp, A) - np.matmul(Lnp, C2)

        return R, N, M, L, H, constraints

    def form_convolutional_matrix(input_signal, T):
        """Find the convolutional matrix for computing convolution.

        In particular, consider the convolution operation, which can
        be re-written as matrix multiplication.

        u * w = U w

        where U is the convolutional matrix.

        Consider a signal x, the output signal y = input_signal * x can be
        computed as the matrix multiplication y = [conv matrix of input_signal] x.

        Let input_signal = u[0],...,u[T-1], the convolutional matrix is
        given by:

        [u[0]                ]
        [u[1], u[0],         ]
        [u[2], u[1], u[0], ..]
        ...
        [u[T-1], ..,     u[0]]
        ...
        [u[N-1], ...         ]

        """
        N = input_signal.shape[0]
        conv_mat = np.zeros((N, T))
        for i in range(T):
            conv_mat[i:, i] = input_signal[:N - i]
        return conv_mat

    def compute_cl_matrix_responses(plant, ctrl, Tsim):
        """ Compute the response mapping.

        Returns:
            R, N, M, L (array): Each is a 3-dimensionl array: (T, nx, nu)
        """
        # Input check
        if plant.dt != ctrl.dt:
            raise(ValueError("Plant and controller do not have same sampling time!"))
        else:
            dT = plant.dt
        if not (isinstance(plant, co.StateSpace) or isinstance(ctrl, co.StateSpace)):
            raise(ValueError("Plant and controller both must be StateSpace!"))

        nu = ctrl.outputs
        ny = ctrl.inputs
        nxc = ctrl.states
        nx = plant.states
        # form the state-space form of the combined closed-loop mappings transfer function
        A, B1, B2, C1, C2, D11, D12, D21, D22 = map(
            np.array, get_partitioned_mats(plant, nu, ny))
        Ak, Bk, Ck, Dk = ctrl.A, ctrl.B, ctrl.C, ctrl.D
        assert np.all(D22 == 0)

        # Concatenate system matrices to find the matrices of the ss rep of the impulse response mapping
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
        resp_dss = co.ss(A_cb, B_cb, C_cb, D_cb, dT)
        # Compute impulse responses of the matrices
        Tarr = np.arange(0, Tsim, dT)  # 5 sec horizon 
        NT = Tarr.shape[0]
        impulses = []
        for i in range(nx + ny):
            _, yarr = co.impulse_response(resp_dss, Tarr, input=i, transpose=True)
            impulses.append(yarr)
        impulse_full = np.stack(impulses, axis=2)

        # Individual responses
        R = impulse_full[:, :nx, :nx]
        N = impulse_full[:, :nx, nx:nx + ny]
        M = impulse_full[:, nx:nx + nu, :nx]
        L = impulse_full[:, nx:nx + nu, nx:nx + ny]

        return R, N, M, L


def get_partitioned_transfer_matrices(Pz_design, nu=1, ny=1):
    """ Return the partitioned transfer matrices.

    Returns:
        Pzw: control.TransferFunction
        Pzu: control.TransferFunction
        Pyw: control.TransferFunction
        Pyu: control.TransferFunction
    """
    nw = Pz_design.inputs - nu
    nz = Pz_design.outputs - ny
    # block matrices
    Pzw = Pz_design[:nz, :nw]
    Pzu = Pz_design[:nz, nw:]
    Pyw = Pz_design[nz:, :nw]
    Pyu = Pz_design[nz:, nw:]
    return Pzw, Pzu, Pyw, Pyu


def plot_step_responses(P: co.TransferFunction, Tmax: float, input_descriptions=range(10), output_descriptions=range(10)):
    """ Simulate and plot all elements of the transfer matrix.

    Args:
        P: MIMO plant
        Tmax: Simulation duration.
    """
    Ts = P.dt
    tsim = np.arange(int(Tmax / Ts)) * Ts
    fig, axs = plt.subplots(P.outputs, P.inputs, sharex=True)
    # plot
    for i in range(P.outputs):
        axs[i, 0].set_ylabel(output_descriptions[i])
        for j in range(P.inputs):
            t, y = co.step_response(P[i, j], tsim)
            axs[i, j].plot(tsim, y[0, :])
    for j in range(P.inputs):
        axs[P.outputs - 1, j].set_xlabel(input_descriptions[j])
    plt.show()


class Qsyn:
    """ Collection of sub-routines used in Q-parametriation based synthesis.
    """
    @staticmethod
    def obtain_time_response_var(weight: cvx.Variable,  io_idx: tuple, T_sim: np.ndarray, Pzw: co.TransferFunction, Pzu: co.TransferFunction, Pyw: co.TransferFunction) -> dict:
        """Return the (i, j) indexed response of the closed-loop response in the time-domain.

        Weight Convention: The basis contains delayed impulses. 

        In particular, let ny and nu be the dimensions of the measure
        output and control input respectively, Q(z) is a nu times ny
        MIMO filter. There are in total of nu * ny filters. Each
        filter is represented by a Finite Impulse Reponse with `Ntaps`
        taps:

        Q_ij = t0 + t1 z^-1 + ... t(N-1) z^(-N+1)

        To represent the complete MIMO filter Q, the taps of the
        filters are concatenated in a row-major order. For example, if
        Q is a 2-by-2 transfer matrix, the coefficients are the
        concatenation of the followings:

        [taps of Q00], [taps of Q01], [taps of Q10], [taps of Q11]

        Hence, the number of weights equals the product of the number
        of taps, ny and nu.

        Args:
            weight: cvxpy Variable containing the taps.
            io_idx: (output, input) indices.
            T_sim: Simulation time indices.
            Pzw: Plant block transfer functions.
            Pzu: Plant block transfer functions.
            Pyw: Plant block transfer functions.

        Returns:
            dictionary of result
        """
        # Post checks
        if not np.allclose(np.diff(T_sim), Pzw.dt):
            raise(ValueError("Inconsistent simulation time `Tsim`"))

        if weight.shape[0] % (Pzu.inputs * Pyw.outputs) != 0:
            raise(ValueError("Inconsistent weight dimension"))

        print(" -- generate time response {:}".format(io_idx))
        i, j = io_idx
        Ts = Pzw.dt
        out = co.impulse_response(Pzw[i, j], T_sim)
        Pzw_resp = out[1][0, :]
        resp = Pzw_resp
        Nsteps = len(T_sim)

        # shape info
        nu = Pzu.inputs
        ny = Pyw.outputs
        Ntaps = int(weight.shape[0] / (nu * ny))

        basic_resp = None
        resp_mat = []
        for k in range(weight.shape[0]):
            iQ = int((k / Ntaps) / ny)  # current row on Q matrix
            jQ = int((k / Ntaps) % ny)  # current col on Q matrix
            delay_steps = k - (ny * iQ + jQ) * Ntaps

            if delay_steps == 0:
                out = co.impulse_response(Pzu[i, iQ] * Pyw[jQ, j], T_sim)
                # assign base response
                basic_resp = out[1][0, :] * Ts
                resp_mat.append(basic_resp)
            else:
                resp_k = np.zeros_like(basic_resp)
                resp_k[delay_steps:] = basic_resp[:(Nsteps - delay_steps)]
                resp_mat.append(resp_k)

        resp_mat = np.array(resp_mat).T
        if type(weight) == cvx.Variable:
            resp = resp + resp_mat * weight
        elif type(weight) == np.ndarray:
            resp = resp + np.dot(resp_mat, weight)
        else:
            raise ValueError("Unknown weight type")
        return resp

    @staticmethod
    def obtain_freq_var(weight: cvx.Variable, io_idx: tuple, freqs: np.ndarray, Pzw: co.TransferFunction, Pzu: co.TransferFunction, Pyw: co.TransferFunction) -> dict:
        """Return a frequency-response variable.

        For convention on weight, see :func:`obtain_time_response_var`.

        Args:
            weight: cvxpy Variable containing the taps.
            io_idx: (output, input) indices.
            freqs: Frequencies at which to evaluate the transfer function.
            Pzw: Plant block transfer functions.
            Pzu: Plant block transfer functions.
            Pyw: Plant block transfer functions.
        """
        # Post checks
        if weight.shape[0] % (Pzu.inputs * Pyw.outputs) != 0:
            raise(ValueError("Inconsistent weight dimension"))
        print(" -- generate frequency reponse {:}".format(io_idx))

        Ts = Pzw.dt
        i, j = io_idx
        mag, phase, _ = Pzw[i, j].freqresp(freqs)
        resp = mag[0, 0] * np.exp(1j * phase[0, 0])
        delay_operator = np.exp(- 1j * freqs * Ts)

        # shape info
        nu = Pzu.inputs
        ny = Pyw.outputs
        Ntaps = int(weight.shape[0] / (nu * ny))

        basic_resp = None
        resp_mat = []
        delay_steps = 0
        for k in range(weight.shape[0]):
            iQ = int((k / Ntaps) / ny)  # current row on Q matrix
            jQ = int((k / Ntaps) % ny)  # current col on Q matrix
            delay_steps = k - (ny * iQ + jQ) * Ntaps

            if delay_steps == 0:
                mag, phase, _ = (Pzu[i, iQ] * Pyw[jQ, j]).freqresp(freqs)
                basic_resp = mag[0, 0] * np.exp(1j * phase[0, 0]) * Ts
                resp_mat.append(basic_resp)
            else:
                delayed_resp = basic_resp * delay_operator ** delay_steps
                resp_mat.append(delayed_resp)
        resp_mat = np.array(resp_mat).T

        if type(weight) == cvx.Variable:
            resp = resp + resp_mat * weight
        elif type(weight) == np.ndarray:
            resp = resp + np.dot(resp_mat, weight)
        else:
            raise ValueError("Unknown weight type")
        return resp

    @staticmethod
    def form_Q_feedback_controller_ss(taps: np.ndarray, Pyu: co.TransferFunction) -> co.StateSpace:
        """Return a state-space realization of the Q controller.

        In Q-parametrization theory, the optimal controller given
        paramtrization Q is given as the feedback interconnection
        between Q and Pyu.

        Args:
            taps: All taps of Q, which are FIRs arranged in array.
            Pyu: Block transfer function of the plant.

        Returns:
            K_ss: Resulting state-space controller.
        """
        nu = Pyu.inputs
        ny = Pyu.outputs
        Ts = Pyu.dt
        Ntaps = int(taps.shape[0] / (nu * ny))

        # form Q transfer function from the taps
        Q_fir_list = []
        for i in range(nu):
            Q_fir_list.append([])
            for j in range(ny):
                current_index = (i * ny + j) * Ntaps
                Q_fir_ij = co.tf(
                    taps[current_index: current_index + Ntaps], [1] + [0] * (Ntaps - 1), Ts)
                Q_fir_list[-1].append(Q_fir_ij)
        Q_fir = tf_blocks(Q_fir_list)

        # compute state-space realizationof Q_fir with matlab
        Q_ss = tf2ss(Q_fir, minreal=True, via_matlab=True)

        K_ss = co.feedback(Q_ss, Pyu, sign=-1)
        return K_ss

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

    @staticmethod
    def Q_synthesis(Pz_design: co.TransferFunction, specs: dict) -> dict:
        """Synthesize a controller for the given plant.

        A dictionary is used to supply specification to the synthesizer.

        Args:
            Pz_design: Plant for controller design.
            specs: Specifications as a dictionary.

        Returns:
            Dictionary containing optimization result.
        """
        nu = specs['nu']
        ny = specs['ny']
        # setup
        Pzw, Pzu, Pyw, Pyu = get_partitioned_transfer_matrices(
            Pz_design, nu=nu, ny=ny)

        Ts = Pz_design.dt
        T_sim = np.arange(specs['Nsteps']) * Ts

        z = co.tf([1, 0], [1], Ts)
        freqs = specs['freqs']
        all_costs = {'shape-time': [], 'shape-freq': [], 'reg2': None, 'reginf': None}
        # for book-keeping only
        freq_vars = {}
        time_vars = {}
        misc = {}

        weight = cvx.Variable(specs['Ntaps'])
        constraints = []
        if 'additional-time-vars' in specs:
            for input_kind, io_idx in specs['additional-time-vars']:
                # time-var
                if (io_idx, input_kind) not in time_vars:
                    imp_var = Qsyn.obtain_time_response_var(
                        weight, io_idx, T_sim, Pzw, Pzu, Pyw)
                    if input_kind == 'step':
                        imp_var = cvx.cumsum(imp_var)
                    elif input_kind == 'step_int':
                        imp_var = cvx.cumsum(cvx.cumsum(imp_var))

                    time_vars[(io_idx, input_kind)] = imp_var

        if 'additional-freq-vars' in specs:
            # func: a function that returns magnitude
            for io_idx in specs['additional-freq-vars']:
                if io_idx in freq_vars:
                    freq_var = freq_vars[io_idx]
                else:
                    freq_var = Qsyn.obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw)
                    freq_vars[io_idx] = freq_var

        # frequency response shaping
        if 'shape-time' in specs:
            for input_kind, io_idx, sys_desired, obj_weight in specs['shape-time']:
                # time-var
                if (io_idx, input_kind) in time_vars:
                    imp_var = time_vars[(io_idx, input_kind)]
                else:
                    imp_var = Qsyn.obtain_time_response_var(weight, io_idx, T_sim, Pzw, Pzu, Pyw)
                    if input_kind == 'step':
                        imp_var = cvx.cumsum(imp_var)
                    elif input_kind == 'step_int':
                        imp_var = cvx.cumsum(cvx.cumsum(imp_var))

                    time_vars[(io_idx, input_kind)] = imp_var

                # desired time-var
                if input_kind == 'step':
                    out = co.step_response(sys_desired, T_sim)
                elif input_kind == 'impulse':
                    out = co.impulse_response(sys_desired, T_sim)
                elif input_kind == 'step_int':
                    out = co.step_response(sys_desired, T_sim)
                    out_1 = np.cumsum(out[1], axis=1) * Ts
                    out = (0, out_1)  # TODO: bad hack, improve this

                imp_desired = out[1][0, :]
                imp_desired[specs['shape-time-delay']:] = (
                    imp_desired[:specs['Nsteps'] - specs['shape-time-delay']])
                imp_desired[:specs['shape-time-delay']] = 0
                all_costs['shape-time'].append(obj_weight * cvx.norm(imp_var - imp_desired))

                # bookkeeping
                time_vars[(io_idx, input_kind, 'desired')] = imp_desired

        # frequency response shaping
        if 'shape-freq' in specs:
            for io_idx, sys_desired, norm, obj_weight in specs['shape-freq']:
                if io_idx in freq_vars:
                    freq_var = freq_vars[io_idx]
                else:
                    freq_var = Qsyn.obtain_freq_var(
                        weight, io_idx, freqs, Pzw, Pzu, Pyw)
                    freq_vars[io_idx] = freq_var
                freq_desired = sys_desired.freqresp(freqs)[0][0, 0]
                all_costs['shape-freq'].append(
                    obj_weight * cvx.norm(freq_var - freq_desired, norm))

        if 'reg2' in specs:
            all_costs['reg2'] = cvx.norm(weight, 2) * specs['reg2']

        if 'reginf' in specs:
            all_costs['reginf'] = cvx.norm(weight, 'inf') * specs['reginf']

        # frequency-domain constraints
        if 'constraint-freq' in specs:
            # func: a function that returns magnitude
            for io_idx, func in specs['constraint-freq']:
                if io_idx in freq_vars:
                    freq_var = freq_vars[io_idx]
                else:
                    freq_var = Qsyn.obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw)
                    freq_vars[io_idx] = freq_var

                # obtain magnitude upper-bound
                if type(func) == co.TransferFunction:
                    mag_upbnd = func.freqresp(freqs)[0][0, 0]
                else:
                    mag_upbnd = func(freqs)

                constraints.append(
                    cvx.abs(freq_var) <= mag_upbnd
                )
                freq_vars[io_idx] = freq_var

        # passivity
        if 'passivity' in specs:
            for io_idx, delta, wc in specs['passivity']:
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
                constraints.extend([
                    cvx.real(freq_var[:idx_wc]) >= 0,
                    cvx.real(freq_var[:idx_wc]) + cvx.imag(freq_var[:idx_wc]) * delta >= 0
                ]
                )
                print(
                    "Passivity constraint accounted for with wc={:f}".format(
                        freqs[idx_wc]))

        # dc-gain
        for io_idx, gain in specs['dc-gain']:
            # only w=0 is needed
            imp_var = Qsyn.obtain_time_response_var(weight, io_idx, T_sim, Pzw, Pzu, Pyw)
            time_vars[(io_idx, 'impulse')] = imp_var
            constraints.append(
                cvx.sum(imp_var) == gain
            )

        # sum all previous found cost
        cost = 0
        for key in all_costs:
            if type(all_costs[key]) == list:
                for cost_ in all_costs[key]:
                    cost += cost_
            elif all_costs[key] is not None:  # scalar cost
                cost += all_costs[key]

        print(" -- Form cvxpy Problem instance")
        prob = cvx.Problem(cvx.Minimize(cost), constraints)

        print(" -- Solve instance")
        prob.solve(verbose=True)
        assert(prob.status == 'optimal')

        print(" -- Cost report\n")
        keys_sorted = sorted(all_costs.keys())
        for key in keys_sorted:
            if type(all_costs[key]) == list:
                for i, cost_ in enumerate(all_costs[key]):
                    print(" > {:} nb. {:d}: {:f}".format(key, i, cost_.value))
            elif all_costs[key] is not None:
                print(" > {:}         : {:f}".format(key, all_costs[key].value))
            else:
                print(" > {:}         : None".format(key))
        print("\n -- Finish reporting cost")

        # return the weight found
        taps = weight.value * Ts

        # form identity matrix Iz
        Iz_list = np.zeros((ny, ny)).tolist()
        for i in range(ny):
            Iz_list[i][i] = z
        Iz = tf_blocks(Iz_list, Ts)
        zPyu = Iz * Pyu
        res = {'Qtaps': taps, 'Pyu': Pyu,
               'zPyu': zPyu,
               'time-vars': time_vars,
               'freq-vars': freq_vars,
               'misc': misc,
               'Tsim': T_sim}
        return res


def analysis(plant, controller, analysis_dict, controller_name='noname',
             in_idxname=range(100), out_idxname=range(100), Nsteps=250):
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
    H = lft(plant, controller)
    nu = controller.outputs
    ny = controller.inputs
    Pzw, Pzu, Pyw, Pyu = get_partitioned_transfer_matrices(
        plant, nu=nu, ny=ny)
    Ts = plant.dt

    # stability
    w, q = np.linalg.eig(H.A)
    wmax = w[np.argmax(np.abs(w))]
    if np.abs(wmax) > 1:
        print(" -- Closed-loop system UNSTABLE. w_max={:}".format(wmax))
    else:
        print(" -- Closed-loop system STABLE. w_max={:}".format(wmax))

    freqs = analysis_dict['freqs']
    T_sim = np.arange(Nsteps) * Ts
    nrow, ncol = analysis_dict['row_col']
    fig, axs = plt.subplots(nrow, ncol, figsize=(10, 10))

    for entry in analysis_dict['recipe']:
        i, j = entry[0], entry[1]

        plot_type = entry[2]

        if plot_type == 'q':
            Q = co.feedback(controller, Pyu, sign=1)
            _, Q_imp = co.impulse_response(Q, np.arange(Nsteps) * Ts)
            Q_imp = Q_imp.flatten()
            axs[i, j].plot(T_sim, Q_imp)
            axs[i, j].set_title(r"$Q (\delta[n]) $")

        elif plot_type == 'step':
            # format: i, j, 'step', (out_idx, in_idx), (out_idx2, in_idx2)
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.step_response(H[output_idx, input_idx], T_sim)
                axs[i, j].plot(T_sim, y_step[0, :],
                               label='step resp. {:}, {:}'.format(
                                   out_idxname[output_idx],
                                   in_idxname[input_idx]))
                axs[i, j].set_title('Step response')

        elif plot_type == 'impulse':
            # format: i, j, 'step', (out_idx, in_idx), (out_idx2, in_idx2)
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.impulse_response(
                    H[output_idx, input_idx], T_sim)
                axs[i, j].plot(T_sim, y_step[0, :],
                               label='imp. resp. {:}, {:}'.format(
                                   out_idxname[output_idx], in_idxname[input_idx]
                               ))
                axs[i, j].set_title('Impulse response')

        elif plot_type == 'step_int':
            # format: i, j, 'step_int', (out_idx, in_idx), (out_idx2, in_idx2)
            # the integral of the step response
            for k in range(3, len(entry)):
                output_idx, input_idx = entry[k]
                _, y_step = co.step_response(H[output_idx, input_idx], T_sim)
                y_step_int = np.cumsum(y_step[0, :]) * Ts
                axs[i, j].plot(T_sim, y_step_int,
                               label='step resp. int. {:}, {:}'.format(
                                   out_idxname[output_idx], in_idxname[input_idx]
                               ))
                axs[i, j].set_title('Step response integral')

        elif plot_type == 'nyquist':
            # format: i, j, 'nyquist', (out_idx, in_idx), (w0, w1, w3) [these
            # are interested frequencies]
            output_idx, input_idx = entry[3]
            mag, phase, freqs = H[output_idx, input_idx].freqresp(freqs)
            H = mag[0, 0] * np.exp(1j * phase[0, 0])
            axs[i, j].plot(H.real, H.imag, '-',
                           label='H{:}{:}'.format(
                               out_idxname[output_idx], in_idxname[input_idx]
                           ))

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
                    label = '{:}, {:}'.format(
                        out_idxname[output_idx], in_idxname[input_idx])
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
                    label = 'H{:}{:}'.format(
                        out_idxname[output_idx], in_idxname[input_idx])
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

