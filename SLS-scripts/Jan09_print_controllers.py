import Jan09_plant_pool as pool
import control as co
import os.path as path
import yaml


def print_controller(controller_name, controller_or_Qdata, scale_output=1e-3):
    """Print a controller to yaml file.

    Controllers can be either a pair of discrete-time filters, or a
    dictionary of optimization result.
    """
    DATA_DIR = "~/catkin_ws/src/infinite_interaction/config/sliding_experiment"
    if type(controller_or_Qdata) == co.TransferFunction:
        assert controller_or_Qdata.inputs == 2
        assert controller_or_Qdata.outputs == 1

        num1 = [float(x) for x in controller_or_Qdata.num[0][0]]
        den1 = [float(x) for x in controller_or_Qdata.den[0][0]]
        num1 = num1 + [0] * (len(den1) - len(num1))

        num2 = [float(x) for x in controller_or_Qdata.num[0][1]]
        den2 = [float(x) for x in controller_or_Qdata.den[0][1]]
        num2 = num2 + [0] * (len(den2) - len(num2))

        controller_data = {
            'type': 'integral',
            'scale_output': scale_output,
            'filter00_b': num1,
            'filter00_a': den1,
            'filter01_b': num2,
            'filter01_a': den2,
        }

    elif type(controller_or_Qdata) == dict:  # Q synthesis dictionary
        Ntaps = int(controller_or_Qdata['Qtaps'].shape[0] / 2)
        taps1 = controller_or_Qdata['Qtaps'][:Ntaps].tolist()
        taps2 = controller_or_Qdata['Qtaps'][Ntaps:].tolist()
        zPyu = controller_or_Qdata['zPyu']

        num1 = [float(x) for x in zPyu.num[0][0]]
        den1 = [float(x) for x in zPyu.den[0][0]]
        num1 = num1 + [0] * (len(den1) - len(num1))

        num2 = [float(x) for x in zPyu.num[1][0]]
        den2 = [float(x) for x in zPyu.den[1][0]]
        num2 = num2 + [0] * (len(den2) - len(num2))

        controller_data = {
            'type': 'Q_synthesis',
            'scale_output': scale_output,
            'filter_Q_00_taps': taps1,
            'filter_Q_01_taps': taps2,
            'filter_zPyu_00_b': num1,
            'filter_zPyu_00_a': den1,
            'filter_zPyu_10_b': num2,
            'filter_zPyu_10_a': den2,
        }

    # save controller to disk
    full_file_path = path.join(path.expanduser(DATA_DIR), controller_name + ".yaml")
    output_dict = {
        "controllers": {controller_name: controller_data}
    }
    with open(full_file_path, "w+") as f:
        yaml.dump(output_dict, f, default_flow_style=False)


if __name__ == "__main__":

    # PI controller integral0_0:
    # This controller works with environments having stiffness around 3N/mm to 4 N/mm
    # Examples: foam sheets, and so on
    # Unstable with stiffer environment
    # controller_to_print = pool.Controllers.PI_v1(0, 12e-1)
    # controller_name = "integral0_0"

    # # PI controller integral1_0:
    # # This controller works with stiff environments having stiffness around 60 N/mm
    # controller_to_print = pool.Controllers.PI_v1(0, 1e-1)
    # controller_name = "integral1_0"

    # PI controller integral2_0:
    # This controller moves in free-space at the same speed as Qsyn1_0
    controller_to_print = pool.Controllers.PI_v1(0, 7.5e-1)
    controller_name = "integral2_0"

    # # PI controller integral3_0:
    # # This controller always produces zero motion
    # controller_to_print = pool.Controllers.PI_v1(0, 0e-1)
    # controller_name = "integral3_0"

    print_controller(controller_name, controller_to_print, scale_output=1e-3)
