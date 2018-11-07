import rospy
import rosbag
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
import cvxpy as cvx
import scipy.signal as signal
import scipy.io


def analysis_freq(extracted_data, cmd, keys):
    """Spectral analysis of given data. 

    Signals are selected using the cmd string.

    Plot the magnitude/angle of selected signals at different
    frequencies. Magnitude is simply the absolute value of the complex
    valued frequency response.

    Args:
        extracted_data (dict): Data.
        cmd (str): String contains information on what signals to extract,
                   as well as the time interval to be used.
        keys(list): Sorted dictionary keys. Use in signal selection.
    """
    try:
        indices_str = cmd.split(" ")[2]
        ts = map(lambda x: float(x), cmd.split(" ")[1].split(','))
        if len(ts) == 2:
            tmin, tmax = ts
        else:
            tmin = ts[0]
            tmax = tmin + 12.56
            # this value is chosen so that the increment in frequency
            # response is exactly 0.5 rad/s
    except Exception as e:
        print e
        return 0
    figure, axs = plt.subplots(2, 1, sharex=True)
    axs[0].set_yscale('log')
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
            # sampling period
            Ts = (extracted_data[key]['t'][imax] - extracted_data[key]['t'][imin]) / (imax - imin)

            # fft, normalized so that the dft is exactly the amplitude
            # of sinuisoidal signals
            N = imax - imin
            dft = np.fft.rfft(extracted_data[key]['data'][imin:imax]) / N * 2
            omegas = np.linspace(0, np.pi, dft.shape[0]) / Ts
            axs[0].scatter(omegas, np.absolute(dft), label="{:}(N={:d})".format(key, N))
            axs[1].scatter(omegas, np.angle(dft), label="{:}(N={:d})".format(key, N))

        elif extracted_data[key]['type'] == 'sensor_msgs/JointState':
            jnt_idx = int(cmd[1:].split('/')[1])
            # get imin, imax
            imin = 0
            imax = len(extracted_data[key]['t']) - 1
            for i in range(len(extracted_data[key]['t'])):
                if extracted_data[key]['t'][i] > tmin and imin == 0:
                    imin = int(i)
                if extracted_data[key]['t'][i] > tmax and imax == len(extracted_data[key]['t']) - 1:
                    imax = int(i)
                    break
            # sampling period
            Ts = (extracted_data[key]['t'][imax] - extracted_data[key]['t'][imin]) / (imax - imin)

            # fft
            N = imax - imin
            dft = np.fft.rfft(extracted_data[key]['data'][imin:imax, jnt_idx]) / N * 2
            omegas = np.linspace(0, np.pi, dft.shape[0]) / Ts
            axs[0].scatter(omegas, np.absolute(dft), label="{:}/{:}(N={:d})".format(key, jnt_idx, N))
            axs[1].scatter(omegas, np.angle(dft), label="{:}/{:}(N={:d})".format(key, jnt_idx, N))

        elif extracted_data[key]['type'] == 'geometry_msgs/WrenchStamped':
            try:
                data_idx = int(cmd[1:].split('/')[1])
            except ValueError as e:
                if cmd[1:].split('/')[1] == "fx":
                    data_idx = 0
                if cmd[1:].split('/')[1] == "fy":
                    data_idx = 1
                if cmd[1:].split('/')[1] == "fz":
                    data_idx = 2
                if cmd[1:].split('/')[1] == "taux":
                    data_idx = 3
                if cmd[1:].split('/')[1] == "taux":
                    data_idx = 4
                if cmd[1:].split('/')[1] == "taux":
                    data_idx = 5
                
            # get imin, imax
            imin = 0
            imax = len(extracted_data[key]['t']) - 1
            for i in range(len(extracted_data[key]['t'])):
                if extracted_data[key]['t'][i] > tmin and imin == 0:
                    imin = int(i)
                if extracted_data[key]['t'][i] > tmax and imax == len(extracted_data[key]['t']) - 1:
                    imax = int(i)
                    break
            # sampling period
            Ts = (extracted_data[key]['t'][imax] - extracted_data[key]['t'][imin]) / (imax - imin)

            # fft
            N = imax - imin
            dft = np.fft.rfft(extracted_data[key]['data'][imin:imax, data_idx]) / N * 2
            omegas = np.linspace(0, np.pi, dft.shape[0]) / Ts
            axs[0].scatter(omegas, np.absolute(dft), label="{:}/{:}(N={:d})".format(key, data_idx, N))
            axs[1].scatter(omegas, np.angle(dft), label="{:}/{:}(N={:d})".format(key, data_idx, N))


    axs[1].set_xlabel("frequency (rad/s)")
    axs[0].set_ylabel("magnitude")
    axs[1].set_ylabel("angle (rad)")
    axs[0].set_xlim(0, 50)
    axs[0].legend()
    axs[1].legend()
    axs[0].grid()
    axs[1].grid()
    plt.show()


def analysis_view(extracted_data, cmd_string, keys):
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


def extract_data_from_bag(DATA_PATH):
    """ Extract data from a bag file from `DATA_PATH`.

    Args
        DATA_PATH (str): Path to data file.

    Returns
        (dict): Extracted data
    """
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
            extracted_data[topic] = {'t': topic_t,
                                     'data': topic_data,
                                     'type': topic_tuple.msg_type}
    return extracted_data


def analysis_z(extracted_data, cmd, keys, gam=0):
    """Analyze data and identify model in discrete time.

    Note: This script analyse interaction between two scalar signal
    only. MIMO is possible, but that is for a future project.

    cmd = "z tmin,tmax order idx1,idx2"

    Params:
        z: Indicator that analysis_z is to be ran.
        tmin, tmax: Two ends of the analysis interval.
        order: Desired Order.
        idx1: Index of the input signal.
        idx2: Index of the output signal.

    """
    tmin, tmax = map(float, cmd.split(" ")[1].split(","))
    Norder = int(cmd.split(" ")[2])
    idx1, idx2 = cmd.split(" ")[3].split(",")
    print("-- Analysis interval: [{:.3f}, {:.3f}]\n"
          "-- Desired Norder: {:d}\n"
          "-- Input/Output indices: {:}, {:}".format(
              tmin, tmax, Norder, idx1, idx2))
    data1 = get_data(extracted_data, idx1, keys)
    data2 = get_data(extracted_data, idx2, keys)

    Nstart1 = np.argmin(np.abs(data1['t'] - tmin))
    Nend1 = np.argmin(np.abs(data1['t'] - tmax))

    Nstart2 = np.argmin(np.abs(data2['t'] - data1['t'][Nstart1]))
    Nend2 = Nstart2 + (Nend1 - Nstart1)

    # notation: Let a=[a0..a[T-1]] and b=[b[0]...b[T-1]] be the
    # coefficients. By definition we have:
    # a[0] x[T] + a[1] x[T-1] + ... + a[T-1] x[0] = b[0] y[T] + b[1] y[T-1] + ... + b[T-1] y[0] 
    # rearrange the equations to obtain:
    # X a = Y b,
    # where X and Y are not the vectors of xs and ys, but a matrix
    # formed by rearranging the vectors appropriately. Note that T is
    # the degree (Norder + 1) of the relevant polynomial coefficients.

    # x[0] is substracted from X as do y[0] is substracted from Y.

    xoffset = 1.0
    X = []
    Y = []
    for i in range(Norder - 1, Nend1 - Nstart1):
        X.append(data1['data'][Nstart1 + i: Nstart1 + i + Norder][::-1] - xoffset)
        Y.append(data2['data'][Nstart2 + i: Nstart2 + i + Norder][::-1])
    X = np.array(X)
    Y = np.array(Y)

    N = Nend1 - Nstart1

    # optimization
    a = cvx.Variable(Norder)
    b = cvx.Variable(Norder)
    obj = cvx.sum_squares(X * a - Y * b) / N
    reg = gam * cvx.norm1(a) + gam * cvx.norm1(b)
    constraints = [b[0] == 1]
    prob = cvx.Problem(cvx.Minimize(obj + reg), constraints)
    prob.solve(solver='MOSEK')

    print("obj: {:f}\nreg: {:f}\na={:}\nb={:}".format(
        obj.value, reg.value, a.value, b.value))

    # validation by running the filter obtained on data
    xin = data1['data'][Nstart1: Nend1] - xoffset
    yact = data2['data'][Nstart2: Nend2]
    bval = np.array(b.value).reshape(-1)
    aval = np.array(a.value).reshape(-1)
    ypred = signal.lfilter(bval, aval, xin)
    plt.plot(- 160 * xin, label='xin')
    plt.plot(yact, label='yactual')
    # plt.plot(ypred, '--', label='ypredict')
    plt.legend()
    plt.show()

    scipy.io.savemat("11_1_J3_human.mat", {'x': xin, 'y': yact})

def get_data(extracted_data, idx_str, keys):
    """ Extract a dictionary using string command.
    """
    if len(idx_str.split(",")) == 1:
        # single index command; e.g. "1" or "2"
        key = keys[int(idx_str)]
        return extracted_data[key]
    elif len(idx_str.split(",")) == 2:
        # double index command; e.g. "1,2" or "2,0"
        key = keys[int(idx_str.split(',')[0])]
        sub_idx = int(idx_str.split(',')[1])
        return {'t': extracted_data[key]['t'],
                'data': extracted_data[key]['data'][:,sub_idx],
                'type': extracted_data[key]['type']}
    else:
        raise NotImplementedError("Unable to parse {:}".format(idx_str))


if __name__ == '__main__':
    # parse command line
    try:
        DATA_PATH = sys.argv[1]
    except:
        print("No bag file given as argument! Exit now!")
        exit(0)

    print("Input Parameters\n -- DATA_PATH: {:}".format(DATA_PATH))

    # load data
    extracted_data = extract_data_from_bag(DATA_PATH)

    # sort dictionary w.r.t keys
    keys = extracted_data.keys()
    keys.sort()

    # zeroing time
    t0 = extracted_data.values()[0]['t'][0]
    for key in keys:
        extracted_data[key]['t'] -= t0

    intro_msg = "\n\nSignals indices:\n"
    for i, key in enumerate(keys):
        intro_msg += " -- Idx: {:2d}, topic: {:}\n".format(i, key)
    print intro_msg

    # processing
    help_msgs = ("Available commands:"
                 "\n -- [h or help]                                "
                 "\n    > to show this help"
                 "\n\n -- [view or v] i,j,k                          "
                 "\n    > plot signal i,j,k"
                 "\n\n -- [freqresp or f] [tmin,tmax or tmin] i,j,k  "
                 "\n    > plot the frequency response of signals i,j,k between tmin and tmax."
                 "\n    > or tmin and tmin + 12.56.")

    # cmd = "z 60,80 3 2,10"
    # analysis_z(extracted_data, cmd, keys)
    # import IPython
    # if IPython.get_ipython() is None:
    #     IPython.embed()

    print(help_msgs)
    while True:
        cmd = raw_input("-\: ")
        if cmd in ["help", "HELP", 'h']:
            print help_msgs

        elif cmd.split(" ")[0] in ["v", "VIEW", 'view']:
            analysis_view(extracted_data, cmd, keys)

        elif cmd.split(" ")[0] in ["f", "freqresp"]:
            analysis_freq(extracted_data, cmd, keys)

        elif cmd.split(" ")[0] in ["z"]:
            pass
