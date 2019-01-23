""" A program for analyzing multiple controllers for force control.
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

    # plants definition:
    # - Link stiffness is identified to be 60 N/mm. To be clear, this stiffness level is
    #   not from the link only, it has effect from the foam table surface as well.
    # - The softest environment is created from those foam sheets. Stiffness level 3.8 N/mm
    # - The stiffest environment are the blocks. Their stiffnesses can be up to 200 N/mm. Because
    #   of this omega_add is assign -30 to check for stability.
    plants = OrderedDict([
        ('nominal', mo.PlantV2.plant(K_env=3.8, omega_add=-30, K_env_aug=40)),
        ('stiff', mo.PlantV2.plant(K_env=20, K_env_aug=10)),
        ('stiffer', mo.PlantV2.plant(K_env=40, K_env_aug=10)),
        ('stiffest', mo.PlantV2.plant(K_env=60, K_env_aug=10)),
        ('stiffest2', mo.PlantV2.plant(K_env=80, K_env_aug=10))
    ])

    # controller descriptions:
    # c0 = mo.Controllers.PI_v1(0, 12e-1)
    # c0 = mo.Controllers.PI_v1(0, 12e-5)
    c0 = mo.Controllers.Qsyn(filename="Jan09_controller_statespace_general_configuration.npz")

    # analysis mode, nothing special, just step output and Nyquist to check the condition of robust stability.
    omega_interested = [1, 5, 10, 20, 40, 80, 100, 200]
    analysis_dict = {
        'row_col': (3, 2),
        'freqs': np.logspace(-3, 2.56, 500),
        'recipe': [
            (0, 0, "step", (0, 0)),
            (0, 1, "step", (0, 1)),
            (1, 0, "nyquist", (3, 3), omega_interested),
            (1, 1, "bode_mag", (1, 0), omega_interested),
            (2, 1, "bode_mag", (0, 0), omega_interested),
            (2, 0, "bode_mag", (2, 0), omega_interested),
        ]
    }
    for plant_key in plants:
        Ss.analysis(plants[plant_key], c0, analysis_dict,
                    input_descriptions=mo.PlantV2.input_descriptions,
                    output_descriptions=mo.PlantV2.output_descriptions,
                    controller_name=plant_key, nb_sim_steps=500)
