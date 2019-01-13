import SLSsyn as Ss
import control as co
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx
import Jan09_model_derivation as mo


Ts = 0.008

if __name__ == '__main__':
    Pz_design = mo.PlantV1.plant()
    step_responses(Pz_design, 0.5)  # preview

    design_dict = {
        'ny': 2,
        'nu': 1,
        'Ntaps': 800,
        'Nsteps': 1000,
        'freqs': np.linspace(1e-2, np.pi / Ts, 1000),
        'resp_delay': 2,  # number of delayed time step

        # DC gain
        'dc-gain': [
            [(0, 0), 1],  # (fm, fd)
            [(0, 1), 0],  # (fm, xe)
        ],

        'freq-constraints': [
            [(0, 2), lambda omegas: np.ones_like(omegas)],
            [(1, 1), lambda omegas: np.ones_like(omegas)],
        ]
    }
    if input("Design controller?") == 'y':
        K_Qparam, data = Ss.Qsyn.Q_synthesis(Pz_design, design_dict)

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()
