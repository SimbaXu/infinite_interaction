import control as co
import numpy as np
import rospy
from std_msgs.msg import Float64MultiArray


class Identificator(object):
    """ Listen to identification data and fitting.
    """
    def __init__(self, base_ns):
        self.base_ns = base_ns
        self.data_subscriber = rospy.Subscriber(
            base_ns + "/identificator/raw_data", Float64MultiArray, self.callback)
        self.data = np.zeros((50, 2))  # (force, displacement) data
        self.date_len = self.data.shape[0]
        self.data_current_index = 0

    def callback(self, msg):
        """

        :type msg: Float64MultiArray
        """
        self.data[self.data_current_index] = msg.data
        self.data_current_index = (self.data_current_index + 1) % self.data_len


if __name__ == '__main__':
    rospy.init_node("logic_node")
    identificator = Identificator("/hybrid_force_controller")
    rospy.spin()

