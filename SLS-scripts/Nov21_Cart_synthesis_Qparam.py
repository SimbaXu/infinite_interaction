"""This script synthesizes a 1-dof controller for 1 direction of the
3-dimensional Cartersian admittance control problem.

"""
import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import Jan09_plant_pool as pool
import yaml


def print_controller(relative_file_path, Qtaps, zPyu_tf):
    """ Print a 3-dimensional DelayedFeedback controller.

    A delayed-feedback controller is defined as

         feedback(Q, z^-1 * zPyu, sign=-1)

    Every dimension shares the same feedback law.
    """
    num = [float(x) for x in zPyu_tf.num[0][0]]
    den = [float(x) for x in zPyu_tf.den[0][0]]
    num = num + [0] * (len(den) - len(num))
    taps = [float(x) for x in np.array(Qtaps).flatten()]
    ctrl_dict = {
        'xff_taps': taps,
        'xfb_b': num,
        'xfb_a': den,
        'yff_taps': taps,
        'yfb_b': num,
        'yfb_a': den,
        'zff_taps': taps,
        'zfb_b': num,
        'zfb_a': den,
    }
    output_dict = {
        'diagonal3D': {
            'filter_dev': ctrl_dict
        }
    }
    with open(relative_file_path, 'w') as f:
        yaml.dump(output_dict, f, default_flow_style=False)
    print("--> Wrote a controller to {:}".format(relative_file_path))


def synthesize_teaching_controller_general_configuration():
    """Synthesize a controller for teaching.

    The output controller is in the general configuration.
    """
    # constant
    Ts = 0.008  # sampling internval
    s = co.tf([1, 0], [1])

    k_env_nom = 40
    m_desired = 1.5
    b_desired = 8
    response_time = 0.4

    # desired response from w3 to z1
    plant_nominal = pool.PlantV3.plant(K_env=40, m_tip=0.07, K_env_aug=1000)

    # Specifications
    noise_atten_func = Ss.Qsyn.lambda_log_interpolate(
        [[0.1, 0.1], [10, 0.1], [20, 0.02], [25, 0.01], [30, 0.01], [50, 0.005], [200, 0.002]], preview=False)

    # fmeasure_xenv_desired_tf = co.c2d(co.tf([m_desired * k_env_nom, b_desired * k_env_nom, 0], [m_desired, b_desired, 0 + k_env_nom]), Ts)
    fmeasure_xenv_desired_tf = co.c2d(40 * response_time * s / (1 + response_time * s), Ts)
    xcmd_xenv_desired_tf = co.c2d(1 / (1 + response_time * s), Ts)

    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 1400,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),

        # different objective
        'shape-time-delay': 3,  # number of delayed time step
        'shape-time': [
            # ['step', (0, 1), fmeasure_xenv_desired_tf, 1]
            ['step', (1, 1), xcmd_xenv_desired_tf, 1]
        ],
        'shape-freq': [
            # [(0, 1), fmeasure_xenv_desired_tf, 'inf', 1.0]
        ],
        'constraint-freq': [
            # [(1, 1), noise_atten_func, False],
        ],

        'reg2': 1,

        # # passivity: Re[H_ij(e^jwT)] >= 0, for all w <= w_c
        # 'passivity': [
        #     [(0, 0), 0.001, 20]
        # ],

        # DC gain
        'dc-gain': [
            [(0, 1), 0]
        ],

        'constraint-nyquist-stability': [
            [(3, 3), (-0.5, 0), 0.57, (2, 40)],  # distance from (-1, 0): 0.5
            [(3, 3), (-0.5, 0), 1.57, (40, 220)]
        ],
        'additional-freq-vars': [
            (3, 3), (1, 1),
        ],
        'additional-time-vars': [
            ["step", (1, 1)],
            ["step", (0, 1)],
        ]
    }

    if input("Design controller?") == 'y':
        synthesis_result = Ss.Qsyn.Q_synthesis(plant_nominal, design_dict)

    if input("View immediate synthesis result?") == 'y':
        omega_interested = [1, 5, 10, 20, 40, 60, 80, 100, 200, 300]
        analysis_dict = {
            'row_col': (3, 2),
            'freqs': design_dict['freqs'],
            'recipe': [
                (0, 1, "step", (0, 1)),
                (0, 1, "step", fmeasure_xenv_desired_tf),
                (0, 0, "step", (1, 1)),
                (1, 0, "nyquist", (3, 3), omega_interested),

                (2, 1, "bode_mag", (1, 1), omega_interested),
                # (2, 1, "bode_mag", (1, 1), omega_interested),
                # (2, 1, "bode_mag", (1, 2), omega_interested),
            ]
        }
        from importlib import reload; reload(Ss)
        Ss.Qsyn.Q_synthesis_analysis(
            synthesis_result, analysis_dict, output_descriptions=pool.PlantV3.output_descriptions,
            input_descriptions=pool.PlantV3.input_descriptions)

    if input("Convert controller to state-space and save for later analysis?") == "y":
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            synthesis_result['Qtaps'], synthesis_result['Pyu'])
        np.savez("Nov21_synthesize_teaching_controller_general_configuration.npz",
                 A=K_Qparam_ss.A, B=K_Qparam_ss.B,
                 C=K_Qparam_ss.C, D=K_Qparam_ss.D, dt=K_Qparam_ss.dt)

    if input("Print controller for execution? y/[n]: ") == "y":
        DATA_DIR='~/catkin_ws/src/infinite_interaction/config/teaching_experiment'
        import Jan09_print_controllers as print_controllers
        from importlib import reload
        reload(print_controllers)
        print_controllers.print_controller("Q_syn_admittance_v0", synthesis_result, scale_output=1, DATA_DIR=DATA_DIR)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    synthesize_teaching_controller_general_configuration()
