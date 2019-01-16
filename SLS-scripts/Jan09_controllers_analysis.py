import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import Jan09_model_derivation as mo


class Controllers:
    @staticmethod
    def PI_v1(Kp, Ki):
        s = co.tf([1, 0], [1])
        K_ctime_single = Kp + Ki / s
        K_dtime_single = co.c2d(K_ctime_single, 0.008)
        PI_v1 = - Ss.tf_blocks([[K_dtime_single, - K_dtime_single]])
        return PI_v1

    @staticmethod
    def Qsyn(filename="Jan09_K_ss_params.npz"):
        data = np.load(filename)
        K = co.StateSpace(data['A'], data['B'], data['C'], data['D'], data['dt'])
        return K


if __name__ == '__main__':
    Ts = 0.008
    s = co.tf([1, 0], [1])

    Pz_design = mo.PlantV1.plant()
    # c0 = Controllers.PI_v1(1e-3, 20e-3)
    c0 = Controllers.Qsyn()

    analysis_dict = {
        'row_col': (2, 2),
        'freqs': np.linspace(0.1, 400, 100),
        'recipe': [
            (0, 0, "step", (0, 0)),
            (0, 1, "step", (0, 1)),
            (1, 0, "step", (1, 1)),
            (1, 1, "bode_mag", (1, 1))
        ]
    }
    Ss.analysis(Pz_design, c0, analysis_dict,
                in_idxname=mo.PlantV1.w_idx_map,
                out_idxname=mo.PlantV1.z_idx_map)
