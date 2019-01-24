import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_plant_pool as mo
import Jan09_print_controllers as print_controllers


def gen_description(Plant, *key):
    out = "{:},{:} > ".format(Plant.z_idx_map[key[0][0]], Plant.w_idx_map[key[0][1]])
    for i in range(1, len(key)):
        out += str(key[i])
    return out


def synthesize_controller_scaled_1dof_configuration():
    """Synthesize a controller in the general controller configuration.

    See drawings/PlanV1.svg for an illustration of this configuration.

    """
    Ts = 0.008
    s = co.tf([1, 0], [1])

    plant_nominal = mo.PlantV2.plant(K_env=3.8e3, K_env_aug=50e3,
                                     x_cmd_scale=1e-4, x_env_scale=1e-4)
    response_time = 0.5

    tracking_desired_resp = co.c2d(1 / (1 + response_time * s), Ts)
    rejection_desired_resp = co.c2d(response_time * s / (1 + response_time * s), Ts)

    def env_aug_shaping_function(freqs):
        centers = np.ones_like(freqs) * np.complex(49, -50)
        radius = np.ones_like(freqs) * 200
        return centers, radius
        # return radius

    env_aug_shaping_function = Ss.Qsyn.lambda_log_interpolate(
        [[0.1, 30], [10, 10], [100, 0.15], [200, 0.05]])


    x_cmd_fdesired_shaping = Ss.Qsyn.lambda_log_interpolate(
        [[10, 1], [100, 0.4]])


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
            ['step', (0, 1), rejection_desired_resp, 1.0],
        ],

        # 'shape-freq': [
        #     [(1, 0), co.tf([0], [1]), 2, 100]
        # ],

        'dc-gain': [
            # [(0, 0), 1],  # (fm, fd)
            # [(0, 1), 0],  # (fm, xe)
        ],

        'reg2': 1e-3,

        'constraint-freq': [
            [(3, 3), env_aug_shaping_function, False],
            [(1, 0), x_cmd_fdesired_shaping, False]
            # [(4, 4), lambda omegas: np.ones_like(omegas) * 50],
        ],

        'constraint-nyquist-stability': [
            # [(3, 3), 0.7, 0.4, (18, 40)],
            # [(3, 3), 0.8, 1.6, (38, 280)]
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
                (2, 1, "bode_mag", x_cmd_fdesired_shaping, omega_interested),
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



def synthesize_controller_general_configuration():
    """ Synthesize a controller in the general controller configuration.

    See drawings/PlanV1.svg for an illustration of this configuration.
    """
    Ts = 0.008
    s = co.tf([1, 0], [1])

    plant_nominal = mo.PlantV2.plant(K_env=5, omega_add=-20, K_env_aug=60)
    response_time = 0.17
    tracking_desired_resp = co.c2d(1 / (1 + response_time * s), Ts)
    rejection_desired_resp = co.c2d(3 * response_time * s / (1 + response_time * s), Ts)
    rejection_force_resp = co.c2d(response_time * s / (1 + response_time * s), Ts)

    def env_aug_shaping_function(freqs):
        centers = np.ones_like(freqs) * np.complex(49, -50)
        radius = np.ones_like(freqs) * 200
        return centers, radius

    env_aug_shaping_function = Ss.Qsyn.lambda_log_interpolate(
        [[0.1, 50], [20, 25], [50, 2], [100, 2], [200, 2]])

    x_cmd_fdesired_shaping = Ss.Qsyn.lambda_log_interpolate(
        # [[10, 0.35], [15, 0.25], [25, 0.15], [100, 0.04], [200, 0.02]])
        [[10, 1], [38, 0.5], [100, 0.1]])
        # [[10, 1e-1], [100, 0.03]])

    x_cmd_fnoise_shaping = Ss.Qsyn.lambda_log_interpolate(
        [[10, 1], [38, 0.5], [100, 0.1]])

    # def env_aug_shaping_function(omegas):
    #     s = co.tf([1, 0], [1])
    #     return 1.05 * (1 / (1 + 0.03 * s) ** 2).freqresp(omegas)[0][0, 0]

    def x_cmd_fdesired_shaping(omegas):
        s = co.tf([1, 0], [1])
        return 0.5 * (1 / (1 + 0.02 * s)).freqresp(omegas)[0][0, 0]

    def x_cmd_fnoise_shaping(omegas):
        s = co.tf([1, 0], [1])
        return 0.5 * (1 / (1 + 0.02 * s)).freqresp(omegas)[0][0, 0]

    def x_cmd_x_env_shaping(omegas):
        s = co.tf([1, 0], [1])
        return 1.2 * (1 / (1 + 0.04 * s)).freqresp(omegas)[0][0, 0]

    f_error_fdesired_upper_bound = Ss.Qsyn.lambda_log_interpolate([[0.1, 1.05], [10, 1.05], [50, 1.05], [200, 1.05]])

    def f_error_fdesired_upper_bound(omegas):
        s = co.tf([1, 0], [1])
        return 1.1 * (s / (2 + s)).freqresp(omegas)[0][0, 0]

    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 1400,
        'Nsteps': 1000,
        'freqs': np.linspace(2, np.pi / Ts, 1000),
        'resp_delay': 3,  # number of delayed time step

        'shape-time-delay': 1,
        'shape-time': [
            ['step', (2, 0), rejection_force_resp, 1.0],
            ['step', (0, 1), rejection_desired_resp, 1.0],
        ],

        # 'shape-freq': [
        #     [(1, 0), co.tf([0], [1]), 2, 100]
        # ],

        'dc-gain': [
            # [(0, 0), 1],  # (fm, fd)
            # [(0, 1), 0],  # (fm, xe)
        ],

        'reg2': 1e-5,

        'constraint-freq': [
            # [(2, 0), f_error_fdesired_upper_bound, False],
            # [(3, 3), env_aug_shaping_function, False],
            [(1, 0), x_cmd_fdesired_shaping, False],
            [(1, 2), x_cmd_fnoise_shaping, False],
            [(1, 1), x_cmd_x_env_shaping, False],
            # [(4, 4), lambda omegas: np.ones_like(omegas) * 50],
        ],

        'constraint-nyquist-stability': [
            [(3, 3), (-0.5, 0), 1.57, (2, 60)],  # distance from (-1, 0): 0.5
            # [(3, 3), (0, 0), 0.643, (8, 35)],  # distance from (-1, 0): 0.6
            # [(3, 3), (0, 0), 0.775, (1, 34)],  # distance from (-1, 0): 0.7
            [(3, 3), (-0.5, 0), 1.57, (60, 220)]
        ],

        'additional-time-vars': [
            ['step', (0, 0)],
            ['step', (1, 0)],
            ['step', (1, 1)]
        ],
        'additional-freq-vars': [
            (4, 4), (0, 0), (1, 0), (3, 3), (1, 2), (2, 0), (1, 1)
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
                (2, 1, "bode_mag", (1, 1), omega_interested),
                (2, 1, "bode_mag", (1, 2), omega_interested),
                (2, 1, "bode_mag", x_cmd_fnoise_shaping),
                (2, 1, "bode_mag", x_cmd_fdesired_shaping),
                (2, 1, "bode_mag", x_cmd_x_env_shaping),

                # (2, 0, 'q')
                (2, 0, 'bode_mag', (2, 0)),
                (2, 0, 'bode_mag', f_error_fdesired_upper_bound)
            ]
        }
        Ss.Qsyn.Q_synthesis_analysis(
            data, analysis_dict, output_descriptions=mo.PlantV2.output_descriptions,
            input_descriptions=mo.PlantV2.input_descriptions)

    if input("Convert controller to state-space and save for later analysis?") == "y":
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            data['Qtaps'], data['Pyu'])
        np.savez("Jan09_controller_statespace_general_configuration.npz",
                 A=K_Qparam_ss.A, B=K_Qparam_ss.B,
                 C=K_Qparam_ss.C, D=K_Qparam_ss.D, dt=K_Qparam_ss.dt)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()

    if input("Print controller for execution") == "y":
        from importlib import reload
        reload(print_controllers)
        print_controllers.print_controller("Q_syn1_0", data, scale_output=1e-3)

if __name__ == '__main__':
    synthesize_controller_general_configuration()
    # synthesize_controller_scaled_1dof_configuration()
