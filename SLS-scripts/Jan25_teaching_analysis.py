"""A program that analyzes multiple controllers for admittance
control.
"""
import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_plant_pool as mo
from collections import OrderedDict


if __name__ == '__main__':
    Ts = 0.008
    s = co.tf([1, 0], [1])
    plants = OrderedDict([
        ('nominal', mo.PlantV3.plant(K_env=30, m_tip=0.07, K_env_aug=200)),
        ('freespace', mo.PlantV3.plant(K_env=1e-3, m_tip=0.07)),
        ('stiffer', mo.PlantV3.plant(m_tip=0.07, K_env=100)),
        ('stiffest', mo.PlantV3.plant(m_tip=0.07, K_env=200)),
        ('spring', mo.PlantV3.plant(m_tip=0.07, K_env=400)),
    ])

    # controller descriptions:
    c0 = mo.Controllers.Admittance(m=2.5, b=12, k=0)
    # c0 = mo.Controllers.Admittance(m=6.0, b=18, k=0)
    # c0 = mo.Controllers.Qsyn(filename="Nov21_synthesize_teaching_controller_general_configuration.npz")

    # analysis mode, nothing special, just step output and Nyquist to check the condition of robust stability.
    omega_interested = [1, 5, 10, 20, 40, 80, 100, 200]
    analysis_dict = {
        'row_col': (2, 2),
        'freqs': np.logspace(-3, 2.56, 500),
        'recipe': [
            (0, 0, "step", (0, 1)),
            (0, 1, "step", (0, 0)),

            (1, 0, "nyquist", (3, 3), omega_interested),
            (1, 1, "bode_mag", (1, 1), omega_interested),
        ]
    }
    for plant_key in plants:
        Ss.analysis(plants[plant_key], c0, analysis_dict,
                    input_descriptions=mo.PlantV2.input_descriptions,
                    output_descriptions=mo.PlantV2.output_descriptions,
                    controller_name=plant_key, nb_sim_steps=500)
