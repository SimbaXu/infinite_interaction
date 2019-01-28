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


def synthesize_teaching_controller_general_configuration():
    """Synthesize a controller for teaching.

    The output controller is in the general configuration.
    """
    # constant
    Ts = 0.008  # sampling internval
    s = co.tf([1, 0], [1])

    k_env_nom = 20
    m_desired = 1.2
    b_desired = 8

    K_env_aug = 1.0

    # desired response from w3 to z1
    plant_nominal = pool.PlantV3.plant(K_env=20, m_tip=0.08, K_env_aug=K_env_aug)

    def xcmd_xenv_upper_bound(omegas):
        s = co.tf([1, 0], [1])
        return (3 / (1 + 0.03 * s) ** 3).freqresp(omegas)[0][0, 0]

    fmeasure_xenv_desired_tf = co.c2d(co.tf([m_desired * k_env_nom, b_desired * k_env_nom, 0], [m_desired, b_desired, 0 + k_env_nom]), Ts)

    f_env_aug_f_robust_tf = Ss.Qsyn.lambda_log_interpolate([
        [0.1, 50], [25, 50], [40, 0.2], [100, 0.05]], preview=False)

    def f_env_aug_f_robust_shaping_for_human(omegas):
        s = co.tf([1, 0], [1])
        # return np.ones_like(omegas)
        return 400 * (1 + 0.157 * 1j * omegas)

    def f_env_aug_f_robust_shaping_for_spring(omegas):
        s = co.tf([1, 0], [1])
        # return np.ones_like(omegas)
        return 800 * (1 + 0.05 * 1j * omegas)

    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 1400,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),

        # different objective
        'shape-time-delay': 6,  # number of delayed time step
        'shape-time': [
            ['step', (0, 1), fmeasure_xenv_desired_tf, 1]
            # ['step', (1, 1), xcmd_xenv_desired_tf, 1]
        ],

        'constraint-freq': [
            # [(1, 1), xcmd_xenv_upper_bound, False],
            ([3, 3, f_env_aug_f_robust_shaping_for_human], f_env_aug_f_robust_tf, False)
        ],

        'reg2': 1,

        # DC gain
        'dc-gain': [
            [(0, 1), 0]
        ],

        'constraint-nyquist-stability': [
            # shaping for robust stability against human
            [(3, 3), f_env_aug_f_robust_shaping_for_human, (-0.0, 0), 0.30, (3, 15)],
            [(3, 3), f_env_aug_f_robust_shaping_for_human, (-0.5, 0), 1.57, (15, 220)],

            # # shaping for robust stability against spring
            # [(3, 3), f_env_aug_f_robust_shaping_for_spring, (-0.0, 0), 0.20, (3, 25)],
            # [(3, 3), f_env_aug_f_robust_shaping_for_spring, (-0.5, 0), 1.57, (25, 220)],

        ],

        'additional-freq-vars': [
            (3, 3), (1, 1), (1, 3), (2, 3)
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
                (1, 0, "nyquist", (3, 3, f_env_aug_f_robust_shaping_for_human), omega_interested),
                (1, 1, "nyquist", (3, 3, f_env_aug_f_robust_shaping_for_spring), omega_interested),

                (2, 0, "bode_mag", (3, 3, f_env_aug_f_robust_shaping_for_human), omega_interested),
                (2, 0, "bode_mag", [f_env_aug_f_robust_tf]),
                
                # (2, 1, "bode_mag", (1, 1), omega_interested),
                # (2, 1, "bode_mag", xcmd_xenv_upper_bound),
                (2, 1, "bode_phs", (1, 3), omega_interested),
                (2, 1, "bode_phs", (2, 3), omega_interested),
                (2, 1, "bode_phs", (3, 3, f_env_aug_f_robust_shaping_for_human), omega_interested),
            ]
        }
        from importlib import reload; reload(Ss)
        Ss.Qsyn.Q_synthesis_analysis(
            synthesis_result, analysis_dict, output_descriptions=pool.PlantV3.output_descriptions,
            input_descriptions=pool.PlantV3.input_descriptions)

    if input("Print controller for execution? y/[n]: ") == "y":
        DATA_DIR='~/catkin_ws/src/infinite_interaction/config/teaching_experiment'
        import Jan09_print_controllers as print_controllers
        from importlib import reload
        reload(print_controllers)
        # profile_name = "Q_syn_admittance_v2"
        profile_name = "Q_syn_admittance_v3_Xaxis"
        print_controllers.print_controller(
            profile_name, synthesis_result, scale_output=1, DATA_DIR=DATA_DIR)

    if input("Convert controller to state-space and save for later analysis?") == "y":
        K_Qparam_ss = Ss.Qsyn.form_Q_feedback_controller_ss(
            synthesis_result['Qtaps'], synthesis_result['Pyu'])
        np.savez("Nov21_synthesize_teaching_controller_general_configuration.npz",
                 A=K_Qparam_ss.A, B=K_Qparam_ss.B,
                 C=K_Qparam_ss.C, D=K_Qparam_ss.D, dt=K_Qparam_ss.dt)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


if __name__ == '__main__':
    synthesize_teaching_controller_general_configuration()
