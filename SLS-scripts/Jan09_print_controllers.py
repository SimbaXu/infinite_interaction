import Jan09_plant_pool as pool
import control as co
import os.path as path
import yaml


def print_controller(controller_name, controller_or_Qdata):
    """Print a controller to yaml file.

    Controllers can be either a pair of discrete-time filters, or a
    dictionary of optimization result.
    """
    DATA_DIR = "~/catkin_ws/src/infinite_interaction/config/hybrid_force_controllers"
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
            'filter00_b': num1,
            'filter00_a': den1,
            'filter01_b': num2,
            'filter01_a': den2,
        }

    elif type(controller_or_Qdata) == dict:  # Q synthesis dictionary
        Ntaps = controller_or_Qdata['Qtaps'].shape[0] / 2
        taps1 = controller_or_Qdata['Qtaps'][:Ntaps]
        taps2 = controller_or_Qdata['Qtaps'][Ntaps:]
        zPyu = controller_or_Qdata['zPyu']

        num1 = [float(x) for x in zPyu.num[0][0]]
        den1 = [float(x) for x in zPyu.den[0][0]]
        num1 = num1 + [0] * (len(den1) - len(num1))

        num2 = [float(x) for x in zPyu.num[1][0]]
        den2 = [float(x) for x in zPyu.den[1][0]]
        num2 = num2 + [0] * (len(den2) - len(num2))

        controller_data = {
            'type': 'Q_synthesis',
            'taps00': taps1,
            'taps01': taps2,
            'zPyu_filter00_b': num1,
            'zPyu_filter00_a': den1,
            'zPyu_filter10_b': num2,
            'zPyu_filter10_a': den2,
        }

    # save controller to disk
    full_file_path = path.join(path.expanduser(DATA_DIR), controller_name + ".yaml")
    output_dict = {
        "controllers": {controller_name: controller_data}
    }
    with open(full_file_path, "w+") as f:
        yaml.dump(output_dict, f, default_flow_style=False)


if __name__ == "__main__":
    # declare name and controller to print
    controller_to_print = pool.Controllers.PI_v1()
    controller_name = "integral0_0"

    print_controller(controller_name, controller_to_print)
