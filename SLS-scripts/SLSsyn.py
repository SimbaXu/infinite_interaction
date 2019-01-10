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


def step_responses(P, Tmax):
    Ts = P.dt
    tsim = np.arange(int(Tmax / Ts)) * Ts
    fig, axs = plt.subplots(P.outputs, P.inputs, sharex=True)
    for i in range(P.outputs):
        for j in range(P.inputs):
            t, y = co.step_response(P[i, j], tsim)
            axs[i, j].plot(tsim, y[0, :])
    plt.show()


class Qsyn:
    """ Collection of sub-routines used in Q-parametriation based synthesis.
    """
    @staticmethod
    def obtain_time_response_var(weight, io_idx, T_sim, Pzw, Pzu, Pyw):
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
        delay_steps = 0
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

    def obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw):
        """Return a frequency-response variable.

        For convention on weight, see :func:`obtain_time_response_var`.


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

    def Q_synthesis(Pz_design, specs):
        """Synthesize a controller for the given plant.

        A dictionary is used to supply specification to the synthesizer.

        [recipe]: A list of lists.
        -

        """
        # shape info
        if 'ny' in specs:
            ny = specs['ny']
        else:
            ny = 1
        if 'nu' in specs:
            nu = specs['nu']
        else:
            nu = 1

        Pzw, Pzu, Pyw, Pyu = get_partitioned_transfer_matrices(
            Pz_design, nu=nu, ny=ny)

        Ts = Pz_design.dt
        T_sim = np.arange(specs['Nsteps']) * Ts

        z = co.tf([1, 0], [1], Ts)
        freqs = specs['freqs']

        # synthesis
        weight = cvx.Variable(specs['Ntaps'] * ny * nu)
        constraints = []

        if 'objective' in specs:
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
        else:
            cost1 = 0

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
        if 'objective-Hinf' in specs:
            for io_idx, func, obj_weight in specs['objective-Hinf']:
                if io_idx in freq_vars:
                    freq_var = freq_vars[io_idx]
                else:
                    freq_var = Qsyn.obtain_freq_var(
                        weight, io_idx, freqs, Pzw, Pzu, Pyw)
                    freq_vars[io_idx] = freq_var
                cost_inf += obj_weight * cvx.norm(freq_var - func(freqs), 'inf')

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
        if 'dc-gain' in specs:
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
        NQ = weight.shape[0]
        diffQ = np.zeros((specs['Ntaps'] - 1, specs['Ntaps']))
        for i in range(specs['Ntaps'] - 1):
            diffQ[i, i: i + 2] = [1.0, -1.0]

        diffY = np.zeros((specs['Nsteps'] - 1, specs['Nsteps']))
        for i in range(specs['Nsteps'] - 1):
            diffY[i, i: i + 2] = [1.0, -1.0]

        imp_weight = np.ones(specs['Nsteps'])
        cost2 = specs['reg'] * cvx.norm1(weight * Ts)
        # cost3 = specs['reg_diff_Q'] * cvx.norm(diffQ * weight * Ts, p=2)
        # cost4 = specs['reg_diff_Y'] * cvx.norm(diffY * imp_var)
        # cost5 = specs['reg'] * cvx.norm1(imp_var - imp_desired)
        cost = cost1 + cost2
        print(" -- Form cvxpy Problem instance")
        prob = cvx.Problem(cvx.Minimize(cost), constraints)
        print(" -- Start solving")
        prob.solve(verbose=True)
        assert(prob.status == 'optimal')

        print("cost: cost1={:f}, cost2={:f},\n      cost3={:f} cost4={:f}"
              "\n     cost5={:f} cost-inf={:f}" .format(
                  cost1.value, cost2.value, 0, 0, 0,
                  cost_inf.value))

        print("Explanations:\n"
              "cost1: |y[n] - y_ideal[n]|_2 * weight + |y[n] - 1|_2 * weight_reg\n"
              "cost2: |Q|_1\n"
              "cost3: |DIFF * Q|_2\n"
              "cost4: |DIFF * y[n]|_2 (for smooth movement)\n"
              "cost5: |y[n] - y_des[n]|_1\n"
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
        Q = tf2ss(Q_fir, minreal=True, via_matlab=True)
        K = co.feedback(Q, Pyu, sign=-1)

        return {'Qtaps': taps, 'Pyu': Pyu, 'zPyu': z * Pyu}

