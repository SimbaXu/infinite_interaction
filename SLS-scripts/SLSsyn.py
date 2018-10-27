import numpy as np
import matlab.engine
import control as co


class MatlabEngine(object):
    """A Singleton for Matlab Engine.

    """
    matlab_engine_instance = None
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


def get_partitioned_mats(P, nu, ny):
    """ Return the partitioned matrices of P.
    """
    if type(P) != co.StateSpace:
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
            if type(Plist[i][j]) == co.xferfcn.TransferFunction:
                nums[i].append(Plist[i][j].num[0][0])
                dens[i].append(Plist[i][j].den[0][0])
            else:
                nums[i].append([0])
                dens[i].append([1])
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
            eng.workspace['num{:d}{:d}'.format(i, j)] = matlab.double(P.num[i][j].tolist())
            eng.workspace['den{:d}{:d}'.format(i, j)] = matlab.double(P.den[i][j].tolist())
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



