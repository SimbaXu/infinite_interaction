import rospy
import rosbag
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys


def analysis_freq(extracted_data, cmd):
    """
    """
    try:
        indices_str = cmd.split(" ")[2]
        tmin, tmax = map(lambda x: float(x), cmd.split(" ")[1].split(','))
    except Exception as e:
        print e
        return 0
    figure, axs = plt.subplots(2, 1, sharex=True)
    for cmd in indices_str.split(','):
        idx = int(cmd.split('/')[0])
        key = keys[idx]
        if extracted_data[key]['type'] == 'std_msgs/Float64':
            # get imin, imax
            imin = 0
            imax = len(extracted_data[key]['t']) - 1
            for i in range(len(extracted_data[key]['t'])):
                if extracted_data[key]['t'][i] > tmin and imin == 0:
                    imin = int(i)
                if extracted_data[key]['t'][i] > tmax and imax == len(extracted_data[key]['t']) - 1:
                    imax = int(i)
                    break
            frequency = (extracted_data[key]['t'][imax] - extracted_data[key]['t'][imin]) / (imax - imin)
            # fft
            N = imax - imin
            dft = np.fft.rfft(extracted_data[key]['data'][imin:imax])
            omegas = np.linspace(0, np.pi, dft.shape[0]) * frequency
            axs[0].scatter(omegas, np.absolute(dft), label=key)
            axs[1].scatter(omegas, np.angle(dft), label=key)

    axs[1].set_xlabel("frequency (rad/s)")
    axs[0].set_ylabel("magnitude")
    axs[1].set_ylabel("angle")
    axs[0].set_xlimits([0, 20])
    plt.show()


def analysis_view(extracted_data, cmd_string):
    """ View and analyze extracted data.
    """

    indices_str = cmd_string.split(" ")[1]
    figure, ax = plt.subplots(1, 1)
    for cmd in indices_str.split(','):
        if cmd[0] == 'd':
            idx = int(cmd[1:].split('/')[0])
            key = keys[idx]
            if extracted_data[key]['type'] == 'std_msgs/Float64':
                ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data']) / 0.008, '^-', label="diff/{:}".format(key))

            elif extracted_data[key]['type'] == 'sensor_msgs/JointState':
                if cmd[1:].split('/')[1] == '0':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 0]) / 0.008, 'x-', label='diff/{:}/j1'.format(key))
                if cmd[1:].split('/')[1] == '1':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 1]) / 0.008, 'x-', label='diff/{:}/j2'.format(key))
                if cmd[1:].split('/')[1] == '2':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 2]) / 0.008, 'x-', label='diff/{:}/j3'.format(key))
                if cmd[1:].split('/')[1] == '3':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 3]) / 0.008, 'x-', label='diff/{:}/j4'.format(key))
                if cmd[1:].split('/')[1] == '4':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 4]) / 0.008, 'x-', label='diff/{:}/j5'.format(key))
                if cmd[1:].split('/')[1] == '5':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 5]) / 0.008, 'x-', label='diff/{:}/j6'.format(key))

            elif extracted_data[key]['type'] == 'geometry_msgs/WrenchStamped':
                if cmd[1:].split('/')[1] == 'fx':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 0]) / 0.008, 'o-', label='{:}/fx'.format(key))
                if cmd[1:].split('/')[1] == 'fy':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 1]) / 0.008, 'o-', label='{:}/fy'.format(key))
                if cmd[1:].split('/')[1] == 'fz':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 2]) / 0.008, 'o-', label='{:}/fz'.format(key))
                if cmd[1:].split('/')[1] == 'taux':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 3]) / 0.008, 'o-', label='{:}/taux'.format(key))
                if cmd[1:].split('/')[1] == 'tauy':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 4]) / 0.008, 'o-', label='{:}/tauy'.format(key))
                if cmd[1:].split('/')[1] == 'tauz':
                    ax.plot(extracted_data[key]['t'][1:], np.diff(extracted_data[key]['data'][:, 5]) / 0.008, 'o-', label='{:}/tauz'.format(key))

        else:
            idx = int(cmd.split('/')[0])
            key = keys[idx]
            if extracted_data[key]['type'] == 'std_msgs/Float64':
                ax.plot(extracted_data[key]['t'], extracted_data[key]['data'], '^-', label=key)

            elif extracted_data[key]['type'] == 'sensor_msgs/JointState':
                jnt_idx = int(cmd[1:].split('/')[1])
                ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, jnt_idx], 'x-', label='{:}/j{:d}'.format(key, jnt_idx + 1))

            elif extracted_data[key]['type'] == 'geometry_msgs/WrenchStamped':
                if cmd.split('/')[1] == 'fx':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 0], 'o-', label='{:}/fx'.format(key))
                if cmd.split('/')[1] == 'fy':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 1], 'o-', label='{:}/fy'.format(key))
                if cmd.split('/')[1] == 'fz':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 2], 'o-', label='{:}/fz'.format(key))
                if cmd.split('/')[1] == 'taux':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 3], 'o-', label='{:}/taux'.format(key))
                if cmd.split('/')[1] == 'tauy':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 4], 'o-', label='{:}/tauy'.format(key))
                if cmd.split('/')[1] == 'tauz':
                    ax.plot(extracted_data[key]['t'], extracted_data[key]['data'][:, 5], 'o-', label='{:}/tauz'.format(key))

    plt.legend()
    plt.show()




if __name__ == '__main__':
    # parse command line
    # DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-18-13-37-07.bag'
    # DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-18-16-11-24-simple.bag'
    # DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-18-slow-move.bag'
    DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-19-freqresp.bag'
    # DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-19-freqresp-250Hz.bag'
    # DATA_PATH = '~/Dropbox/ros_data/2018-Sep-15-admittance-data/2018-09-19-freqresp-500Hz.bag'

    print("Input Parameters\n -- DATA_PATH: {:}".format(DATA_PATH))
    # load data
    print("Start extracting data")
    bag = rosbag.Bag(os.path.expanduser(DATA_PATH))
    extracted_data = {}
    for topic in bag.get_type_and_topic_info()[1].keys():
        topic_tuple = bag.get_type_and_topic_info()[1][topic]

        if topic_tuple.msg_type == 'sensor_msgs/JointState':
            print(" -- Found topic [{:}], type: [{:}]".format(topic, topic_tuple.msg_type))
            topic_t = []
            topic_data = []
            for _, msg, t in bag.read_messages(topic):
                topic_t.append(t.to_sec())
                topic_data.append(msg.position)
            topic_t = np.array(topic_t)
            topic_data = np.array(topic_data)
            extracted_data[topic] = {'t': topic_t,
                                     'data': topic_data,
                                     'type': topic_tuple.msg_type}

        if topic_tuple.msg_type == 'std_msgs/Float64':
            print(" -- Found topic [{:}], type: [{:}]".format(topic, topic_tuple.msg_type))
            topic_t = []
            topic_data = []
            for _, msg, t in bag.read_messages(topic):
                topic_t.append(t.to_sec())
                topic_data.append(msg.data)
            topic_t = np.array(topic_t)
            topic_data = np.array(topic_data)
            extracted_data[topic] = {'t': topic_t,
                                     'data': topic_data,
                                     'type': topic_tuple.msg_type}

        if topic_tuple.msg_type == 'geometry_msgs/WrenchStamped':
            print(" -- Found topic [{:}], type: [{:}]".format(topic, topic_tuple.msg_type))
            topic_t = []
            topic_data = []
            for _, msg, t in bag.read_messages(topic):
                topic_t.append(t.to_sec())
                topic_data.append(
                    (msg.wrench.force.x, msg.wrench.force.y, msg.wrench.force.z,
                     msg.wrench.torque.x, msg.wrench.torque.y, msg.wrench.torque.z)
                )
            topic_t = np.array(topic_t)
            topic_data = np.array(topic_data)
            import ipdb; ipdb.set_trace()
            extracted_data[topic] = {'t': topic_t,
                                     'data': topic_data,
                                     'type': topic_tuple.msg_type}

    # zeroing
    t0 = extracted_data.values()[0]['t'][0]
    for key in extracted_data.keys():
        extracted_data[key]['t'] -= t0

    intro_msg = "\n\nSignals indices:\n"
    keys = extracted_data.keys()
    for i, key in enumerate(keys):
        intro_msg += " -- Idx: {:d}, topic: {:}\n".format(i, key)
    print intro_msg

    # processing
    help_msgs = "Available commands:\n -- VIEW i,j,k   : plot signal i,j,k"
    print(help_msgs)
    while True:
        cmd = raw_input("cmd: ")
        if cmd in ["help", "HELP", 'h']:
            print help_msgs

        elif cmd.split(" ")[0] in ["v", "VIEW", 'view']:
            analysis_view(extracted_data, cmd)

        elif cmd.split(" ")[0] in ["f", "freqresp"]:
            analysis_freq(extracted_data, cmd)
