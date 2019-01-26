import control as co
import numpy as np
import os.path as path
import yaml
import SLSsyn as Ss
from Jan09_print_controllers import print_controller


def generate_and_print(profile_name, m, b, k, Ts=0.008):
    """Generate an admittance controller from the given setting and print.

    Note that the first input is reference force, which is always
    zero. The second input is force measured with sensor.

    """
    print("""
Printing admittance controller
------------------------------

  > name: {:}
  > m:={:.1f}, b:={:.1f}, k:={:.1f}
  > discretization time step: Ts:={:.3f}
""".format(profile_name, m, b, k, Ts))
    s = co.tf([1, 0], [1])
    admittance_model = 1.0 / (m * s ** 2 + b * s + k)
    admittance_model_discrete = co.c2d(admittance_model, Ts)
    controller = Ss.tf_blocks([[0, admittance_model_discrete]])
    print(controller)

    print_controller(profile_name, controller, scale_output=1,
                     DATA_DIR='~/catkin_ws/src/infinite_interaction/config/teaching_experiment')
    


if __name__ == '__main__':
    # parameters and profile

    # profile_name = 'admittance_v0'
    # m = 2.5
    # b = 12
    # k = 0

    profile_name = 'admittance_v1'
    m = 6
    b = 18
    k = 0

    generate_and_print(profile_name, m, b, k)


