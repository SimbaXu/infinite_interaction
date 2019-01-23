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

    plant_nominal = mo.PlantV2.plant(K_env=3.8e3, omega_add=-20, K_env_aug=60e3)
    response_time = 0.5
    tracking_desired_resp = co.c2d(1 / (1 + response_time * s), Ts)
    rejection_desired_resp = co.c2d(3500 * response_time * s / (1 + response_time * s), Ts)

    def env_aug_shaping_function(freqs):
        centers = np.ones_like(freqs) * np.complex(49, -50)
        radius = np.ones_like(freqs) * 200
        return centers, radius
        # return radius


    env_aug_shaping_function = Ss.Qsyn.lambda_log_interpolate(
        [[0.1, 30], [10, 10], [100, 0.15], [200, 0.05]])


    x_cmd_fdesired_shaping = Ss.Qsyn.lambda_log_interpolate(
        [[10, 2e-4], [100, 2e-5]])


    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 1400,
        'Nsteps': 1000,
        'freqs': np.linspace(2, np.pi / Ts, 1000),
        'resp_delay': 3,  # number of delayed time step

        'shape-time-delay': 1,
        'shape-time': [
            ['step', (0, 0), tracking_desired_resp, 1],
            ['step', (0, 1), rejection_desired_resp, 0.1 / 3500],
        ],

        # 'shape-freq': [
        #     [(1, 0), co.tf([0], [1]), 2, 100]
        # ],

        'dc-gain': [
            # [(0, 0), 1],  # (fm, fd)
            # [(0, 1), 0],  # (fm, xe)
        ],

        'reg2': 1,

        'constraint-freq': [
            [(3, 3), env_aug_shaping_function, False],
            # [(1, 0), x_cmd_fdesired_shaping]
            # [(4, 4), lambda omegas: np.ones_like(omegas) * 50],
        ],

        'constraint-nyquist-stability': [
            [(3, 3), 0.7, 0.4, (18, 40)],
            [(3, 3), 0.8, 1.6, (38, 280)]
        ],

        'additional-time-vars': [
            ['step', (1, 0)],
            ['step', (1, 1)]
        ],
        'additional-freq-vars': [
            (4, 4), (0, 0), (1, 0), (3, 3)
        ]
    }

    if input("Design controller?") == 'y':
        data = Ss.Qsyn.Q_synthesis(plant_nominal, design_dict)

    if input("View immediate synthesis result?") == 'y':
        omega_interested = [1, 5, 10, 20, 40, 60, 80, 100, 200, 300]
        analysis_dict = {
            'row_col': (3, 2),
            'freqs': design_dict['freqs'],
            'recipe': [
                (0, 0, "step", (0, 0)),
                (0, 1, "step", (0, 1)),
                (1, 0, "nyquist", (3, 3), omega_interested),

                (1, 1, "bode_mag", (3, 3), omega_interested),
                (1, 1, "bode_mag", env_aug_shaping_function, omega_interested),

                (2, 1, "bode_mag", (1, 0), omega_interested),
                # (2, 0, 'q')
                (2, 0, 'bode_phs', (3, 3))
            ]
        }
        Ss.Qsyn.Q_synthesis_analysis(data, analysis_dict, output_descriptions=mo.PlantV1.output_descriptions, input_descriptions=mo.PlantV1.input_descriptions)

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
