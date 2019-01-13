import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_model_derivation as mo


def gen_description(Plant, *key):
    out = "{:},{:} > ".format(Plant.z_idx_map[key[0][0]], Plant.w_idx_map[key[0][1]])
    for i in range(1, len(key)):
        out += str(key[i])
    return out


if __name__ == '__main__':
    Ts = 0.008
    s = co.tf([1, 0], [1])

    Pz_design = mo.PlantV1.plant()

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
            ['step', (1, 0)],
            ['step', (1, 1)]
        ],
        'additional-freq-vars': [
            (0, 1), (2, 3)
        ]
    }

    if input("Design controller?") == 'y':
        data = Ss.Qsyn.Q_synthesis(Pz_design, design_dict)

    if input("View immediate synthesis result?") == 'y':
        # debug: preview result
        imp00 = data['time-vars'][((0, 0), 'impulse')].value

        # scale by xe magnitude to a step of 1mm
        imp01 = data['time-vars'][((0, 1), 'impulse')].value * 0.001
        step_des00 = data['time-vars'][((0, 0), 'step', 'desired')]

        plt.plot(data['Tsim'], np.cumsum(imp00), label=gen_description(mo.PlantV1, (0, 0), 'step'))
        plt.plot(data['Tsim'], np.cumsum(imp01), label=gen_description(mo.PlantV1, (0, 1), 'step'))
        plt.plot(data['Tsim'], step_des00, label=gen_description(mo.PlantV1, (0, 0), 'step', 'desired'))
        plt.plot(data['Tsim'], data['time-vars'][((1, 0), 'step')].value, label=gen_description(mo.PlantV1, (1, 0), 'step'))
        plt.plot(data['Tsim'], data['time-vars'][((1, 1), 'step')].value, label=gen_description(mo.PlantV1, (1, 1), 'step'))
        plt.legend()
        plt.show()

        io2plot = [(1, 1)]
        for io in io2plot:
            mag = np.abs(data['freq-vars'][io].value)
            plt.plot(design_dict['freqs'], mag, label=gen_description(mo.PlantV1, io))
        plt.legend()
        plt.show()

    if input("Save controller for later analysis?") == "y":
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            data['Qtaps'], data['Pyu'])
        np.savez("Jan09_K_ss_params.npz", A=K_Qparam_ss.A, B=K_Qparam_ss.B,
                 C=K_Qparam_ss.C, D=K_Qparam_ss.D, dt=K_Qparam_ss.dt)
        

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()
