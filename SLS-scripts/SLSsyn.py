import numpy as np
import matlab.engine
import control as co
import yaml
import os
import cvxpy as cvx
import scipy.sparse as sparse


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


def get_num_dens(Plist):
    """Extract the numerators and denominators of a list of SISO TFs.

    The list of lists of SISO TFs is assumed to have the following
    form:

    [[P11, P12, ...],
     [P21, P22, ...],
     [P21, P22, ...],
     ...]
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

            # other than tf function, must be a float
            elif isinstance(Plist[i][j], float) or isinstance(Plist[i][j], int):
                nums[i].append([float(Plist[i][j])])
                dens[i].append([1])
            else:
                raise(ValueError(
                    "Input: {:} not understood. Must be a transfer function or a numeric value."))

    return nums, dens


def tf_blocks(Plist):
    """Create a MIMO TF from multiple SISO TFs.

    See :func:`get_num_dens` for details.
    """
    nums, dens = get_num_dens(Plist)
    return co.tf(nums, dens)


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

    def form_SLS_response_matrices(Pssd, nu, ny, T):
        """Form response matrices and output matrices.

        Return matrices {R, N, M, L} respectively and output matrices
        {H}. Also return all achievability constraints that guarantee
        the mappings are valid.

        Return:
            R, N, M, L (cvxpy.Variables): Closed-loop reponse matrices.
            H (cvxpy.Variable): Closed-loop output from exogeneous input to exogeneous output.
            constraints (list of cvxpy.Constraint): Achievability constraints.

        """
        nx = Pssd.states
        A, B1, B2, C1, C2, D11, D12, D21, D22 = get_partitioned_mats(
            Pssd, nu, ny)

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
