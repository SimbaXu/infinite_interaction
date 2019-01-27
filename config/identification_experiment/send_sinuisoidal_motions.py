import rospy
import numpy as np
from std_msgs.msg import Float64MultiArray


if __name__ == '__main__':
    rospy.init_node("command_center")
    # timing
    Ts = 0.008
    rate = rospy.Rate(125)

    # parameter
    x_radius = 0.01
    y_radius = 0
    z_radius = 0

    freqs = 10.0  # rad/s

    setpoints = np.zeros(3)
    setpoints_publisher = rospy.Publisher(
        "/identification_node/setpoints", Float64MultiArray, queue_size=1)

    index = 0
    time_now = rospy.get_time()
    while(not rospy.is_shutdown()):

        setpoints[0] = x_radius * np.sin(freqs * rospy.get_time() - time_now)
        setpoints[1] = y_radius * np.sin(freqs * rospy.get_time() - time_now)
        setpoints[2] = z_radius * np.sin(freqs * rospy.get_time() - time_now)
        print(setpoints)

        setpoints_msg = Float64MultiArray()
        setpoints_msg.data = setpoints

        setpoints_publisher.publish(setpoints_msg)

        index += 1
        rate.sleep()
