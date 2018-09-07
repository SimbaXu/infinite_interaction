#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Float64.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <infinite_interaction_lib.h>
// openrave
#include <openrave-core.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

using std::cout;


void SetViewer(OpenRAVE::EnvironmentBasePtr penv, const std::string& viewername)
{
    OpenRAVE::ViewerBasePtr viewer = RaveCreateViewer(penv,viewername);
    BOOST_ASSERT(!!viewer);

    // attach it to the environment:
    penv->Add(viewer);

    // finally call the viewer's infinite loop (this is why a separate thread is needed)
    bool showgui = true;
    viewer->main(showgui);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "admittance_node");
    ros::NodeHandle node_handle("~");
    // parameter
    std::string name_space = "denso";
    std::string ft_sensor_topic = "/netft/ft_sensor/raw";
    std::string joint_state_topic = "/denso/joint_states";
    std::string scenefilename = "robots/denso_handle.robot.xml";
    std::string viewername = "qtosg";
    std::string robotname = "denso_handle";
    // FTSensor
    FTSensorHandler ft_handler;
    ros::Subscriber ft_sub = node_handle.subscribe(ft_sensor_topic, 3, &FTSensorHandler::signal_callback, &ft_handler);
    // position handler and command
    JointPositionHandler position_handler;
    ros::Subscriber position_handler_sub = node_handle.subscribe(joint_state_topic, 3, &JointPositionHandler::signal_callback, &position_handler);
    // command
    JointPositionController position_ctrl(name_space, node_handle);
    std::vector<double> position_cmd = {0, 0, 0, 0, 0, 0};
    // openrave
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr penv = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Debug);
    boost::thread thviewer(boost::bind(SetViewer,penv,viewername));
    penv->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot;
    robot = penv->GetRobot(robotname);
    ros::Rate rate(125);
    while(!ros::isShuttingDown()){
        // retrieve the latest wrench measurement, offseted and filtered

        // project the wrench to the joint space
        robot->SetActiveDOFValues(position_handler.joint_position);

        // run the wrench through the joint filter
        rate.sleep();
        ROS_DEBUG_STREAM("In control loop " << position_handler.joint_position[0]);
        ROS_DEBUG_STREAM(ft_handler.fx << ft_handler.fy);
    }

  return 0;
}


