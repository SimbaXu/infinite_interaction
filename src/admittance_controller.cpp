#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>

class FTSensor {
public:
    double fx, fy, fz;
    FTSensor(): fx(0), fy(0), fz(0) {};
    void ft_callback(const geometry_msgs::WrenchStampedConstPtr& msg);

};

void FTSensor::ft_callback(const geometry_msgs::WrenchStampedConstPtr &msg) {
  fx = msg->wrench.force.x;
  fy = msg->wrench.force.y;
  fz = msg->wrench.force.z;

}

void ft_callback(const geometry_msgs::WrenchStampedConstPtr& msg){
  ROS_INFO_STREAM(msg->wrench.force.x);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "sim_proc");
  ros::NodeHandle n("~");
  FTSensor ft_sensor;
  ros::Subscriber ft_sub = n.subscribe("/netft/ft_sensor/raw", 3, &FTSensor::ft_callback, &ft_sensor);
  auto t = ros::Duration(0.1);
  while (!ros::isShuttingDown()){
    ros::spinOnce();
    ROS_INFO_STREAM(ft_sensor.fx << ", " << ft_sensor.fy);
    t.sleep();
  }
  return 0;
}


