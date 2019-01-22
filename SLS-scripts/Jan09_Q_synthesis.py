import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_plant_pool as mo


def gen_description(Plant, *key):
    out = "{:},{:} > ".format(Plant.z_idx_map[key[0][0]], Plant.w_idx_map[key[0][1]])
    for i in range(1, len(key)):
        out += str(key[i])
    return out


def synthesize_controller_general_configuration():
    """ Synthesize a controller in the general controller configuration.

    See drawings/PlanV1.svg for an illustration of this configuration.
    """
    Ts = 0.008
    s = co.tf([1, 0], [1])

    plant_nominal = mo.PlantV1.plant(K_env=3.8e3, omega_add=-30, k_link=60e3)
    desired_resp = co.c2d(1 / (1 + 0.5 * s), Ts)

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
            # [(0, 2), lambda omegas: np.ones_like(omegas)],
            # [(1, 1), lambda omegas: np.ones_like(omegas)],
        ],

        'additional-time-vars': [
            ['step', (1, 0)],
            ['step', (1, 1)]
        ],
        'additional-freq-vars': [
            (4, 4)
        ]
    }

    if input("Design controller?") == 'y':
        data = Ss.Qsyn.Q_synthesis(plant_nominal, design_dict)

    if input("View immediate synthesis result?") == 'y':
        omega_interested = [1, 5, 10, 20, 40]
        analysis_dict = {
            'row_col': (2, 2),
            'freqs': np.logspace(-3, 2.2, 500),
            'recipe': [
                (0, 0, "step", (0, 0)),
                (0, 1, "step", (1, 0)),
                (1, 0, "nyquist", (4, 4), omega_interested),
                (1, 1, "bode_mag", (4, 4), omega_interested),
            ]
        }

        Ss.Qsyn.Q_synthesis_analysis(data, analysis_dict)

    if input("Save controller for later analysis?") == "y":
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            data['Qtaps'], data['Pyu'])
        np.savez("Jan09_controller_statespace_general_configuration.npz", A=K_Qparam_ss.A, B=K_Qparam_ss.B,
                 C=K_Qparam_ss.C, D=K_Qparam_ss.D, dt=K_Qparam_ss.dt)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    synthesize_controller_general_configuration()
