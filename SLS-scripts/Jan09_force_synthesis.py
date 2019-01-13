import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import Jan09_model_derivation as mo


def gen_description(*key):
    z_idx_map = {
        0: 'fm', 1: 'x'
    }
    w_idx_map = {
        0: 'fd', 1: 'xe', 2: 'n'
    }
    out = "{:},{:} > ".format(z_idx_map[key[0][0]], w_idx_map[key[0][1]])
    for i in range(1, len(key)):
        out += str(key[i])
    return out


if __name__ == '__main__':
    Ts = 0.008
    s = co.tf([1, 0], [1])

    Pz_design = mo.PlantV1.plant()
    # step_responses(Pz_design, 0.5)  # preview

    desired_resp = co.c2d(1 / (1 + 0.2 * s), Ts)

    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 800,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),
        'resp_delay': 2,  # number of delayed time step

        'shape-time-delay': 1,
        'shape-time': [
            ['step', (0, 0), desired_resp, 1],
        ],

        'dc-gain': [
            [(0, 0), 1],  # (fm, fd)
            [(0, 1), 0],  # (fm, xe)
        ],

        'constraint-freq': [
            [(0, 2), lambda omegas: np.ones_like(omegas)],
            [(1, 1), lambda omegas: np.ones_like(omegas)],
        ],

        'additional-time-vars': [
            ['step', (1, 0)]
        ],
        'additional-freq-vars': [
            (0, 1)
        ]
    }
    if input("Design controller?") == 'y':
        data = Ss.Qsyn.Q_synthesis(Pz_design, design_dict)

        # preview result
        imp00 = data['time-vars'][((0, 0), 'impulse')].value

        # scale by xe magnitude to a step of 1mm
        imp01 = data['time-vars'][((0, 1), 'impulse')].value * 0.001
        step_des00 = data['time-vars'][((0, 0), 'step', 'desired')]

        plt.plot(data['Tsim'], np.cumsum(imp00), label=gen_description((0, 0), 'step'))
        plt.plot(data['Tsim'], np.cumsum(imp01), label=gen_description((0, 1), 'step'))
        plt.plot(data['Tsim'], step_des00, label=gen_description((0, 0), 'step', 'desired'))
        plt.plot(data['Tsim'], data['time-vars'][((1, 0), 'step')].value, label=gen_description((1, 0), 'step'))
        plt.legend()
        plt.show()

        io2plot = [(0, 1)]
        for io in io2plot:
            mag = np.abs(data['freq-vars'][io].value)
            plt.plot(design_dict['freqs'], mag, label=gen_description(io))
        plt.legend()
        plt.show()

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()
