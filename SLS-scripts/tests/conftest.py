import numpy as np
import control as co


def random_dtime_sys(Nout, Nin, Ts, Nstate=5):
    np.random.seed(0)
    sys_ = co.drss(Nstate, Nin, Nout)
    sys_.dt = Ts
    sys = co.tf(sys_)
    return sys
