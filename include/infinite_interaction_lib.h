//
// Created by hung on 8/9/18.
//

#ifndef PROJECT_INFINITE_INTERACTION_LIB_H
#define PROJECT_INFINITE_INTERACTION_LIB_H

#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float64.h>
#include <sensor_msgs/JointState.h>
#include <vector>

class FTSensorHandler {
public:
    double fx, fy, fz;
    FTSensorHandler(): fx(0), fy(0), fz(0) {};
    void signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg);

};

class JointPositionHandler{
public:
    std::vector<double> joint_position;
    JointPositionHandler();
    void signal_callback(const sensor_msgs::JointStateConstPtr &msg);

};

JointPositionHandler::JointPositionHandler() {
    joint_position.resize(6);
}

void JointPositionHandler::signal_callback(const sensor_msgs::JointStateConstPtr &msg) {
    for (int i = 0; i < 6; ++i) {
        joint_position[i] = msg->position[i];
    }
}

class JointPositionController {
public:
    /*! Constructor for a joint position controller class.
     */
    explicit JointPositionController(std::string name_space, ros::NodeHandle& nh);
    void set_joint_positions(std::vector<double>& jnt_positions);

private:
    std::string _name_space;
    std::vector<ros::Publisher> _jnt_pubs;
    ros::NodeHandle _nh;

};

JointPositionController::JointPositionController(std::string name_space, ros::NodeHandle &nh) : _name_space(name_space), _nh(nh) {
    for (int i=0; i < 6; i++){
        std::string jnt_topic = "/" + _name_space + "/j" + std::to_string(i + 1) + "/command";
        ros::Publisher jnt_pub = _nh.advertise<std_msgs::Float64>(jnt_topic, 5);
        _jnt_pubs.push_back(jnt_pub);
    }
}

void JointPositionController::set_joint_positions(std::vector<double>& jnt_positions) {
    for(int i=0; i < 6; i++){
        std_msgs::Float64 msg;
        msg.data = jnt_positions[i];
        _jnt_pubs[i].publish(msg);
    }
}

void FTSensorHandler::signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg) {
    fx = msg->wrench.force.x;
    fy = msg->wrench.force.y;
    fz = msg->wrench.force.z;

}
#endif //PROJECT_INFINITE_INTERACTION_LIB_H
