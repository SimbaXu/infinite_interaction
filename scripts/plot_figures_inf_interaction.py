import rosbag
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
import os.path as path

# set plot style with seaborn
import seaborn as sns
sns.set('paper')
sns.set_style('dark')


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

    sns.set('paper', style='dark', font_scale=0.9)

    ###########################################################################
    #                                 Task 1 plot                                 #
    ###############################################################################
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3.5, 3))
    idx_Q = 7843
    idx_Ad_sameforce = 7840
    idx_Ad_stable = 7330
    time_array = np.arange(-50, -50 + 600) * 0.008
    
    axs[0].plot(time_array, df_master[fg_Q & fg_force]['y'].values[idx_Q: idx_Q + 600],
                label="CCSa")
    axs[0].plot(time_array,
        df_master[fg_Ad_sameforce & fg_force]['y'].values[idx_Ad_sameforce: idx_Ad_sameforce+600],
        label="CLa1")
    axs[0].plot(time_array,
        df_master[fg_Ad_stable & fg_force]['y'].values[idx_Ad_stable: idx_Ad_stable + 600],
        label="CLa2")
    axs[0].grid()
    axs[0].set_ylim([-2, 9])

    axs[1].plot(time_array, df_master[fg_Q & fg_position]['y'].values[idx_Q: idx_Q + 600], label="CCSa")
    axs[1].plot(time_array, df_master[fg_Ad_sameforce & fg_position]['y'].values[idx_Ad_sameforce: idx_Ad_sameforce+600], label="CLa1")
    axs[1].plot(time_array, df_master[fg_Ad_stable & fg_position]['y'].values[idx_Ad_stable: idx_Ad_stable + 600], label="CLa2")
    axs[1].legend()
    axs[1].grid()
    axs[1].set_ylim([-0.4, 0.4])
    plt.tight_layout()
    plt.savefig("2019_Jan_teaching_experiment.pdf")
    plt.show()

    ###########################################################################
    #                                 Task 2 plot                                 #
    ###############################################################################
    idx_Q = 11649
    idx_Ad_sameforce = 11845
    idx_Ad_stable = 12237
    idx_len = 450
    time_array = np.arange(idx_len) * 0.008

    thres_Q = 0.11
    thres_Ad_sameforce = 0.197
    thres_Ad_stable = 0.158

    fig, axs = plt.subplots(1, 1, figsize=(3.5, 1.75))
    plt.plot(time_array, df_master[fg_Q & fg_position]['z'].values[idx_Q: idx_Q + idx_len] + thres_Q, label="CCSa")
    plt.plot(time_array, (df_master[fg_Ad_sameforce & fg_position]['z'].values[idx_Ad_sameforce: idx_Ad_sameforce + idx_len] + thres_Ad_sameforce) * 0.8, label="CLa1")
    plt.plot(time_array, df_master[fg_Ad_stable & fg_position]['z'].values[idx_Ad_stable: idx_Ad_stable + idx_len] + thres_Ad_stable, label="CLa2")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("2019_Jan_teaching_experiment_touchdown.pdf")
    plt.show()


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
    fg_position = df_master['topic'] == command_topic

    idx_Q =   500 + 600 + 774 - 250
    idx_PI1 = 500 + 0 + 1054 - 250 + 35
    idx_PI2 = 500 + 3260 + 847 - 250 - 30
    idx_len = 50 * 125 + 300
    time_array = np.arange(idx_len) * 0.008

    vert_Q = -0.0201
    vert_PI1 = -0.0044
    vert_PI2 = -0.0078

    sns.set("paper", style="dark", font_scale=0.9)
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3.5, 3))
    # plot stuff
    axs[0].plot(time_array, df_master[fg_Q2 & fg_force]['z'].values[idx_Q: idx_Q + idx_len],
                label="CCSf", zorder=10.0)
    axs[0].plot(time_array, df_master[fg_PI_samespeed & fg_force]['z'].values[idx_PI1: idx_PI1 + idx_len],
                label="CLf1")
    axs[0].plot(time_array, df_master[fg_PI_stiff & fg_force]['z'].values[idx_PI2: idx_PI2 + idx_len],
                label="CLf2")

    axs[0].set_ylim(-2, 20)
    axs[0].plot([-2, 55], [5, 5], '--', c='C4')
    axs[0].set_yticks([0, 5, 10, 15, 20])
    axs[0].grid()

    axs[1].plot(time_array, df_master[fg_Q2 & fg_position]['z'].values[idx_Q: idx_Q + idx_len] - vert_Q,
                label="CCSf", zorder=10.0)
    axs[1].plot(time_array, df_master[fg_PI_samespeed & fg_position]['z'].values[idx_PI1: idx_PI1 + idx_len] - vert_PI1,
                label="CLf1")
    axs[1].plot(time_array, df_master[fg_PI_stiff & fg_position]['z'].values[idx_PI2: idx_PI2 + idx_len] - vert_PI2,
                label="CLf2")
    axs[1].legend()
    axs[1].grid()
    plt.tight_layout()
    plt.savefig("2019_Jan_hybrid_force_control.pdf")
    plt.show()


def make_superimposed_videos_teaching():
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

    ###########################################################################
    #                                 Task 1 plot                                 #
    ###############################################################################

    sns.set(context="paper", style='dark')
    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(3.2, 1.2))
    idx_Q = 7843
    idx_Ad_sameforce = 7840
    idx_Ad_stable = 7330
    idx_len = 1800
    time_array = np.arange(idx_len) * 0.008

    # video settings
    idx_wait = 240  # 2 seconds

    line_Q, = axs.plot(time_array,
                       df_master[fg_Q & fg_force]['y'].values[idx_Q: idx_Q + idx_len],
                       label="[Qa]")

    line_Ad1, = axs.plot(time_array,
        df_master[fg_Ad_sameforce & fg_force]['y'].values[idx_Ad_sameforce: idx_Ad_sameforce+idx_len],
        label="[Ad1]")
    line_Ad2, = axs.plot(time_array,
        df_master[fg_Ad_stable & fg_force]['y'].values[idx_Ad_stable: idx_Ad_stable + idx_len],
        label="[Ad2]")
    axs.grid()
    axs.set_ylim([-10, 10])
    axs.set_yticks([-10, -5, 0, 5, 10])
    axs.set_xticklabels([])
    axs.set_title("Force in the Y direction (N)")
    plt.tight_layout(pad=0.6, w_pad=0)

    # For animating the moving horizon
    line, = axs.plot([0], [0], lw=2, zorder=11)

    def animate(i):
        # Change line appearance
        if i < (idx_len + idx_wait):
            line_Q.set_linewidth(.5)
            line_Ad1.set_linewidth(0.5)
            line_Ad2.set_linewidth(2.0)

            line_Q.set_zorder(1)
            line_Ad1.set_zorder(1)
            line_Ad2.set_zorder(10)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([df_master[fg_Ad_stable & fg_force]['y'].values[idx_Ad_stable + time_idx] + 0.5,
                            df_master[fg_Ad_stable & fg_force]['y'].values[idx_Ad_stable + time_idx] - 0.5])

        elif i < 2 * (idx_len + idx_wait):
            line_Q.set_linewidth(0.5)
            line_Ad1.set_linewidth(2)
            line_Ad2.set_linewidth(0.5)

            line_Q.set_zorder(1)
            line_Ad1.set_zorder(10)
            line_Ad2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([df_master[fg_Ad_sameforce & fg_force]['y'].values[idx_Ad_sameforce + time_idx] + 0.5,
                            df_master[fg_Ad_sameforce & fg_force]['y'].values[idx_Ad_sameforce + time_idx] - 0.5])

        elif i < 3 * (idx_len + idx_wait):
            line_Q.set_linewidth(2.0)
            line_Ad1.set_linewidth(0.5)
            line_Ad2.set_linewidth(0.5)

            line_Q.set_zorder(10)
            line_Ad1.set_zorder(1)
            line_Ad2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([df_master[fg_Q & fg_force]['y'].values[idx_Q + time_idx] + 0.5,
                            df_master[fg_Q & fg_force]['y'].values[idx_Q + time_idx] - 0.5])

        return (line, line_Q, line_Ad1, line_Ad2)

    ani = animation.FuncAnimation(fig, animate, interval=64, blit=True,
                                  frames=range(0, 3 * (idx_len + idx_wait), 8))
    ani.save("clip_teaching_part1.mp4", codec="png", dpi=300,
             savefig_kwargs={'transparent': False, 'facecolor': (1, 1, 1, 0)})
    plt.show()

    ###########################################################################
    #                                 Task 2 plot                                 #
    ###############################################################################
    idx_Q = 11549
    idx_Ad_sameforce = 10945
    idx_Ad_stable = 12137
    idx_len = 1240
    time_array = np.arange(idx_len) * 0.008

    # video settings
    idx_wait = 240  # 2 seconds

    thres_Q = 0.11
    thres_Ad_sameforce = 0.197
    thres_Ad_stable = 0.158

    Q_data = df_master[fg_Q & fg_position]['z'].values[idx_Q: idx_Q + idx_len] + thres_Q
    Ad1_data = df_master[fg_Ad_sameforce & fg_position]['z'].values[idx_Ad_sameforce: idx_Ad_sameforce + idx_len] + thres_Ad_sameforce
    Ad2_data = df_master[fg_Ad_stable & fg_position]['z'].values[idx_Ad_stable: idx_Ad_stable + idx_len] + thres_Ad_stable

    fig, axs = plt.subplots(1, 1, figsize=(3.2, 1.2))
    line_Q, = plt.plot(time_array, Q_data, label="[Qa]")
    line_Ad1, = plt.plot(time_array, Ad1_data, label="[Ad1]")
    line_Ad2, = plt.plot(time_array, Ad2_data, label="[Ad2]")
    plt.grid()
    axs.set_xticklabels([])

    axs.set_title("Position in the Z direction (m)")
    plt.tight_layout(pad=0.6, w_pad=0)

    line, = axs.plot([0], [0], lw=2, zorder=11)

    def animate(i):
        # Change line appearance
        if i < (idx_len + idx_wait):
            line_Q.set_linewidth(.5)
            line_Ad1.set_linewidth(0.5)
            line_Ad2.set_linewidth(2.0)

            line_Q.set_zorder(1)
            line_Ad1.set_zorder(1)
            line_Ad2.set_zorder(10)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Ad2_data[time_idx] + 0.01,
                            Ad2_data[time_idx] - 0.01])

        elif i < 2 * (idx_len + idx_wait):
            line_Q.set_linewidth(0.5)
            line_Ad1.set_linewidth(2)
            line_Ad2.set_linewidth(0.5)

            line_Q.set_zorder(1)
            line_Ad1.set_zorder(10)
            line_Ad2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Ad1_data[time_idx] + 0.01,
                            Ad1_data[time_idx] - 0.01])

        elif i < 3 * (idx_len + idx_wait):
            line_Q.set_linewidth(2.0)
            line_Ad1.set_linewidth(0.5)
            line_Ad2.set_linewidth(0.5)

            line_Q.set_zorder(10)
            line_Ad1.set_zorder(1)
            line_Ad2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Q_data[time_idx] + 0.01,
                            Q_data[time_idx] - 0.01])

        return (line, line_Q, line_Ad1, line_Ad2)

    ani = animation.FuncAnimation(fig, animate, interval=64, blit=True,
                                  frames=range(0, 3 * (idx_len + idx_wait), 8))
    ani.save("clip_teaching_part2.mp4", codec="png", dpi=300,
             savefig_kwargs={'transparent': False, 'facecolor': (1, 1, 1, 0)})
    plt.show()


def make_superimposed_video_hybrid_force_control():
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
    fg_position = df_master['topic'] == command_topic

    idx_Q =   500 + 600 + 774 - 250 + 200
    idx_PI1 = 500 + 0 + 1054 - 250 + 35 + 200
    idx_PI2 = 500 + 3260 + 847 - 250 - 30 + 200
    idx_len = 30 * 125
    # idx_len = 300
    time_array = np.arange(idx_len) * 0.008

    idx_wait = 240  # 2 seconds
    
    Q_data = df_master[fg_Q2 & fg_force]['z'].values[idx_Q: idx_Q + idx_len]
    Pi1_data = df_master[fg_PI_samespeed & fg_force]['z'].values[idx_PI1: idx_PI1 + idx_len]
    Pi2_data = df_master[fg_PI_stiff & fg_force]['z'].values[idx_PI2: idx_PI2 + idx_len]

    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(3.2, 1.2))
    # plot stuff
    line_Q, = axs.plot(time_array, Q_data, label="[Qi]", zorder=10.0)
    line_Pi1, = axs.plot(time_array, Pi1_data, label="[Pi1]")
    line_Pi2, = axs.plot(time_array, Pi2_data, label="[Pi2]")

    axs.set_ylim(-2, 20)
    axs.set_yticks([0, 5, 10, 15, 20])

    axs.set_title("Force in the Z direction (N)")
    plt.tight_layout(pad=0.6, w_pad=0)

    line, = axs.plot([0], [0], lw=2, zorder=11)
    axs.grid()

    def animate(i):
        # Change line appearance
        if i < (idx_len + idx_wait):
            line_Q.set_linewidth(.5)
            line_Pi1.set_linewidth(0.5)
            line_Pi2.set_linewidth(2.0)

            line_Q.set_zorder(1)
            line_Pi1.set_zorder(1)
            line_Pi2.set_zorder(10)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Pi2_data[time_idx] + 0.5,
                            Pi2_data[time_idx] - 0.5])

        elif i < 2 * (idx_len + idx_wait):
            line_Q.set_linewidth(0.5)
            line_Pi1.set_linewidth(2)
            line_Pi2.set_linewidth(0.5)

            line_Q.set_zorder(1)
            line_Pi1.set_zorder(10)
            line_Pi2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Pi1_data[time_idx] + 0.5,
                            Pi1_data[time_idx] - 0.5])

        elif i < 3 * (idx_len + idx_wait):
            line_Q.set_linewidth(2.0)
            line_Pi1.set_linewidth(0.5)
            line_Pi2.set_linewidth(0.5)

            line_Q.set_zorder(10)
            line_Pi1.set_zorder(1)
            line_Pi2.set_zorder(1)

            time_idx = max(i % (idx_len + idx_wait) - idx_wait, 0)
            line.set_xdata([time_idx * 0.008, time_idx * 0.008])
            line.set_ydata([Q_data[time_idx] + 0.5,
                            Q_data[time_idx] - 0.5])

        return (line, line_Q, line_Pi1, line_Pi2)

    ani = animation.FuncAnimation(fig, animate, interval=64, blit=True,
                                  frames=range(0, 3 * (idx_len + idx_wait), 8))
    ani.save("animation_force_experiment.mp4", codec="png", dpi=300,
             savefig_kwargs={'transparent': False, 'facecolor': (1, 1, 1, 0)})
    plt.show()


if __name__ == '__main__':
    plot_teaching_experiment_figures()


