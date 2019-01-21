""" A program for analyzing multiple controllers for force control.
"""
import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_plant_pool as mo


if __name__ == '__main__':
    Ts = 0.008
    s = co.tf([1, 0], [1])

    # with some stuff
    # Pz_design = mo.PlantV1.plant(K_env=500e3, k_tip=500, b_tip=30)
    Pz_design = mo.PlantV1.plant(K_env=500, K_env_aug=4000)
    c0 = mo.Controllers.PI_v1(0, 10e-3)
    # c0 = Controllers.Qsyn()

    analysis_dict = {
        'row_col': (2, 2),
        'freqs': np.logspace(-3, 2, 500),
        'recipe': [
            (0, 0, "step", (0, 0)),
            (0, 1, "step", (1, 0)),
            (1, 0, "step", (2, 1)),
            (1, 1, "nyquist", (3, 3)),
        ]
    }
    Ss.analysis(Pz_design, c0, analysis_dict,
                in_idxname=mo.PlantV1.input_descriptions,
                out_idxname=mo.PlantV1.output_descriptions)
