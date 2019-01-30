import rosbag
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path as path


def plot_teaching_experiment_figures():
    DATA_DIR = "~/Dropbox/ros_data/2018-2019-Convex-controller-synthesis/2019-experiment-nb3-teaching/bagfiles"

    bagfiles = {
        "Q": '2019-01-28-18-56-43-Qsyn.bag',
        "Ad_sameforce": '2019-01-28-19-17-25-Admittance-sameforce-good.bag',
        "Ad_stable": '2019-01-28-18-44-26-Admittance-stable.bag',
    }

    command_topic = "/debugger/cartesian_command"
    position_topic = "/debugger/cartesian_position_measurement"
    force_topic = "/debugger/cartesian_force_measurement"
    all_topic = [command_topic, position_topic, force_topic]

    data = []
    for controller_id in bagfiles:
        bag = rosbag.Bag(path.join(path.expanduser(DATA_DIR), bagfiles[controller_id]))
        for msg in bag.read_messages(all_topic):
            data.append([msg.topic, controller_id, msg.timestamp.to_sec()] + list(msg.message.data[:3]))
    df_master = pd.DataFrame(data=data, columns=['topic', 'controller', 'time', 'x', 'y', 'z'])

    # flag to get data
    fg_Q = df_master['controller'] == "Q"
    fg_Ad_sameforce = df_master['controller'] == "Ad_sameforce"
    fg_Ad_stable = df_master['controller'] == "Ad_stable"
    fg_force = df_master['topic'] == force_topic
    fg_position = df_master['topic'] == position_topic

    # plot stuff
    plt.plot(df_master[fg_Q & fg_force]['y'].values, label="Q")
    plt.plot(df_master[fg_Ad_sameforce & fg_force]['y'].values, label="Ad_sameforce")
    plt.plot(df_master[fg_Ad_stable & fg_force]['y'].values, label="Ad stable")
    plt.legend()
    plt.show()

    plt.plot(df_master[fg_Q & fg_position]['y'].values, label="Q")
    plt.plot(df_master[fg_Ad_sameforce & fg_position]['y'].values, label="Ad_sameforce")
    plt.plot(df_master[fg_Ad_stable & fg_position]['y'].values, label="Ad stable")
    plt.legend()
    plt.show()

    plt.plot(df_master[fg_Q & fg_position]['z'].values, label="Q")
    plt.plot(df_master[fg_Ad_sameforce & fg_position]['z'].values, label="Ad_sameforce")
    plt.plot(df_master[fg_Ad_stable & fg_position]['z'].values, label="Ad stable")
    plt.legend()
    plt.show()

    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()


def plot_force_control_experiment_figures():
    DATA_DIR = "~/Dropbox/ros_data/2018-2019-Convex-controller-synthesis/2019-experiment-nb2-force-control/bagfiles"

    bagfiles = {
        "Q": '2019-01-25-15-36-29-Qsyn.bag',
        "Q2": '2019-01-25-20-08-41-Qsyn-take2.bag',
        "PI_samespeed": '2019-01-25-15-44-08-integral-samespeed.bag',
        "PI_stiff": '2019-01-25-15-39-57-integral-stiff.bag',
    }

    command_topic = "/debugger/cartesian_command"
    position_topic = "/debugger/cartesian_position_measurement"
    force_topic = "/debugger/cartesian_force_measurement"
    all_topic = [command_topic, position_topic, force_topic]

    data = []
    for controller_id in bagfiles:
        bag = rosbag.Bag(path.join(path.expanduser(DATA_DIR), bagfiles[controller_id]))
        for msg in bag.read_messages(all_topic):
            data.append([msg.topic, controller_id, msg.timestamp.to_sec()] + list(msg.message.data[:3]))
    df_master = pd.DataFrame(data=data, columns=['topic', 'controller', 'time', 'x', 'y', 'z'])
    # flag to get data
    fg_Q = df_master['controller'] == "Q"
    fg_Q2 = df_master['controller'] == "Q2"
    fg_PI_samespeed = df_master['controller'] == "PI_samespeed"
    fg_PI_stiff = df_master['controller'] == "PI_stiff"
    fg_force = df_master['topic'] == force_topic
    fg_position = df_master['topic'] == position_topic

    # plot stuff
    plt.plot(df_master[fg_Q & fg_force]['z'].values, label="Q")
    plt.plot(df_master[fg_Q2 & fg_force]['z'].values, label="Q take 2")
    plt.plot(df_master[fg_PI_samespeed & fg_force]['z'].values, label="PI samespeed")
    plt.plot(df_master[fg_PI_stiff & fg_force]['z'].values, label="PI stiff")
    plt.legend()
    plt.show()

    
if __name__ == '__main__':
    plot_teaching_experiment_figures()


