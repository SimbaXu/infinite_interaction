import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt


Ts = 0.008

def plant(kind='v1', tau_R1 = 0.0437, ke=10e3, br=10, kr=500, mr=0.05):
    # constant transfer function
    z = co.tf([1, 0], [1], Ts)
    s = co.tf([1, 0], [1])
    zero = co.tf([0], [1], Ts)
    one = co.tf([1], [1], Ts)

    R1 = 1 / (1 + tau_R1 * s)
    if kind == 'v1':
        """
        outputs = [z1, x, y1, y2]
        inputs = [fd, xe, n, u1]
        """
        P = Ss.tf_blocks([
            [0, z ** (-2) * co.c2d(ke*(br*s + kr)/(br*s + ke + kr + mr*s**2), Ts), 0, z ** (-4) * co.c2d(-R1*(br*ke*s + br*mr*s**3 + ke*kr + kr*mr*s**2)/(br*s + ke + kr + mr*s**2), Ts)],
            [0, 0, 0, z ** (-2) * co.c2d(R1, Ts)],
            [1,                                                    0, 0,                                                                                     0],
            [0, z ** (-2) * co.c2d(ke*(br*s + kr)/(br*s + ke + kr + mr*s**2), Ts), 1, z ** (-4) * co.c2d(-R1*(br*ke*s + br*mr*s**3 + ke*kr + kr*mr*s**2)/(br*s + ke + kr + mr*s**2), Ts)]
        ])
    return P


class Qsyn:
    """ Collection of sub-routines used in Q-parametriation based synthesis.
    """
    @staticmethod
    def obtain_time_response_var(
            weight, io_idx, input_kind, T_sim, Pzw, Pzu, Pyw):
        """Form an indexed response in the time-domain.

        A basis contains of delayed impulses is used. 

        For i in [0,...,N-1]:    Q_i = [Ts * z^i, 0],
        For i in [N, 2N -1]:     Q_i = [0, Ts * z^(i - N)].

        The weights are therefore understood as the two FIR taps
        placed in parallel.

        """
        print(" -- generate time response {:}".format(io_idx))
        i, j = io_idx
        if input_kind == 'step' or input_kind == 'step_int':
            out = co.step_response(Pzw[i, j], T_sim)
        elif input_kind == 'impulse':
            out = co.impulse_response(Pzw[i, j], T_sim)
        Pzw_resp = out[1][0, :]
        resp = Pzw_resp
        Nsteps = len(T_sim)
        Ntaps = weight.shape[0] / 2

        basic_resp_1 = None
        basic_resp_2 = None
        resp_mat = []
        for k in range(weight.shape[0]):
            if k == 0:
                if input_kind == 'step' or input_kind == 'step_int':
                    out = co.step_response(Pzu[i, 0] * Pyw[0, j], T_sim)
                elif input_kind == 'impulse':
                    out = co.impulse_response(Pzu[i, 0] * Pyw[0, j], T_sim)
                else:
                    raise(
                        ValueError(
                            "Unknown input kind {:}".format(input_kind)))
                # assign base response
                basic_resp_1 = out[1][0, :] * Ts
                resp_mat.append(basic_resp_1)
            elif k == Ntaps:
                if input_kind == 'step' or input_kind == 'step_int':
                    out = co.step_response(Pzu[i, 0] * Pyw[1, j], T_sim)
                elif input_kind == 'impulse':
                    out = co.impulse_response(Pzu[i, 0] * Pyw[1, j], T_sim)
                else:
                    raise(
                        ValueError(
                            "Unknown input kind {:}".format(input_kind)))
                # assign base response
                basic_resp_2 = out[1][0, :] * Ts
                resp_mat.append(basic_resp_2)

            elif 0 < k and k < Ntaps:
                resp_k = np.zeros_like(basic_resp_1)
                resp_k[k:] = basic_resp_1[:(Nsteps - k)]
                resp_mat.append(resp_k)

            elif k > Ntaps and k < 2 * Ntaps:
                resp_k = np.zeros_like(basic_resp_2)
                resp_k[k:] = basic_resp_2[:(Nsteps - k)]
                resp_mat.append(resp_k)

            else:
                assert False, "Error occurs"

        resp_mat = np.array(resp_mat).T
        resp = resp + resp_mat * weight
        if input_kind == 'step_int':
            resp = cvx.cumsum(resp) * Ts
        return resp

    def obtain_freq_var(weight, io_idx, freqs, Pzw, Pzu, Pyw):
        """Return a frequency-response variable.

        The frequency-response variable is the sum of basis
        responses. The basis Q(s) are delayed impulse with magnitude
        Ts.

        """
        print(" -- generate frequency reponse {:}".format(io_idx))
        i, j = io_idx
        mag, phase, _ = Pzw[i, j].freqresp(freqs)
        freqresp = mag[0, 0] * np.exp(1j * phase[0, 0])
        delay_operator = np.exp(- 1j * freqs * Ts)
        Ntaps = int(weight.shape[0] / 2)

        basic_resp_1 = None
        basic_resp_2 = None
        resp_mat = []
        for k in range(weight.shape[0]):
            if k == 0:
                mag, phase, _ = (Pzu[i, 0] * Pyw[0, j]).freqresp(freqs)
                basic_resp_1 = mag[0, 0] * np.exp(1j * phase[0, 0]) * Ts
                resp_mat.append(basic_resp_1)
            elif k < Ntaps:
                delayed_resp = basic_resp_1 * delay_operator ** k
                resp_mat.append(delayed_resp)
            elif k == Ntaps:
                mag, phase, _ = (Pzu[i, 0] * Pyw[1, j]).freqresp(freqs)
                basic_resp_2 = mag[0, 0] * np.exp(1j * phase[0, 0]) * Ts
                resp_mat.append(basic_resp_2)
            elif k > Ntaps and k < 2 * Ntaps:
                delayed_resp = basic_resp_2 * delay_operator ** k
                resp_mat.append(delayed_resp)
        resp_mat = np.array(resp_mat).T

        # check if the frequency condition is for DC gain
        if len(freqs) == 1 and freqs[0] == 1e-5:
            import ipdb
            ipdb.set_trace()
        freqresp += resp_mat * weight
        return freqresp

def step_responses(P, Tmax):
    tsim = np.arange(int(Tmax / Ts)) * Ts
    fig, axs = plt.subplots(P.outputs, P.inputs, sharex=True)
    for i in range(P.outputs):
        for j in range(P.inputs):
            t, y = co.step_response(P[i, j], tsim)
            axs[i, j].plot(tsim, y[0, :])
    plt.show()

def Q_synthesis(P, design_dict):
    """
    """
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(
        Pz_design, nu=1, ny=2)
    T_sim = np.arange(specs['Nsteps']) * Ts
    z = co.tf([1, 0], [1], Ts)
    freqs = specs['freqs']

    # synthesis
    weight = cvx.Variable(2 * specs['Ntaps'])
    constraints = []

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


if __name__ == '__main__':
    P = plant()
    # step_responses(P, 0.5)  # preview

    design_dict = {
        'Ntaps': 800,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),
        'resp_delay': 2,  # number of delayed time step
        'objective': ['step', (2, 2), desired_sys, 5],
       
    }
    if input("Design controller?") == 'y':
        K_Qparam, data = Q_synthesis(Pz_design, design_dict)



    
import IPython
if IPython.get_ipython() is None:
    IPython.embed()

