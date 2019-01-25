import rospy
import numpy as np
from std_msgs.msg import Float64MultiArray


if __name__ == '__main__':
    rospy.init_node("command_center")
    commands = rospy.get_param("/hybrid_force_controller/commands")
    desired_force = commands['desired_force']

    # timing
    Ts = 0.008
    rate = rospy.Rate(125)

    # handle waypoints
    waypoints = np.array(commands['waypoints'])
    duration = waypoints[-1, 0]
    t_array = np.arange(0, duration, Ts)
    x_array = np.interp(t_array, waypoints[:, 0], waypoints[:, 1])
    y_array = np.interp(t_array, waypoints[:, 0], waypoints[:, 2])

    setpoints = np.zeros(3)
    setpoints_publisher = rospy.Publisher(
        "/hybrid_force_controller/setpoints", Float64MultiArray, queue_size=1)

    index = 0
    while(not rospy.is_shutdown() and index < t_array.shape[0]):
        setpoints[0] = x_array[index]
        setpoints[1] = y_array[index]
        setpoints[2] = desired_force

        setpoints_msg = Float64MultiArray()
        setpoints_msg.data = setpoints

        setpoints_publisher.publish(setpoints_msg)

        index += 1
        rate.sleep()

