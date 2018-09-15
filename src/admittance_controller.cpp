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
    ros::init(argc, argv, "admittance_control");
    ros::NodeHandle node_handle("~");
    // parameter
    std::string name_space = "denso";
    std::string ft_sensor_topic = "/netft/ft_sensor/raw";
    std::string joint_state_topic = "/denso/joint_states";
    std::string scenefilename = "robots/denso_handle.robot.xml";
    std::string viewername = "qtosg";
    std::string robotname = "denso_handle";
    std::vector<double> wrench_ref;
    node_handle.getParam("ft_sensor_ref", wrench_ref);
    std::vector<double> position_ref;
    node_handle.getParam("joint_pos_ref", position_ref);
    if (wrench_ref.size() != 6){
        ROS_ERROR_STREAM("Reference wrench from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    if (position_ref.size() != 6){
        ROS_ERROR_STREAM("Reference position from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    // FTsensor handler
    FTSensorHandler ft_handler(wrench_ref);
    ros::Subscriber ft_sub = node_handle.subscribe(ft_sensor_topic, 3, &FTSensorHandler::signal_callback, &ft_handler);
    // position handler
    JointPositionHandler position_handler;
    ros::Subscriber position_handler_sub = node_handle.subscribe(joint_state_topic, 3, &JointPositionHandler::signal_callback, &position_handler);
    ros::Duration(0.5).sleep(); ros::spinOnce(); // update current joint position
    if (!position_handler.received_msg()){
        ROS_FATAL("Have not received messages to update initial joint position! \n-- Terminating!");
        exit(0);
    }
    // position controller and command
    JointPositionController position_act(name_space, node_handle);
    dVector initial_jnt_position = position_handler.get_latest_jnt_position();
    ROS_INFO_STREAM("Initial position: " << initial_jnt_position[0] << ", " << initial_jnt_position[1] << ", "
                            << initial_jnt_position[2] << ", " << initial_jnt_position[3] << ", "
                            << initial_jnt_position[4] << ", " << initial_jnt_position[5]);
    // position command
    std::vector<DiscreteTimeFilter> admittance_filters;
    for (int i = 0; i < 6; ++i) {
        dVector cof_a, cof_b;
        node_handle.getParam("j" + std::to_string(i + 1) + "/tf_f2q_a", cof_a);
        node_handle.getParam("j" + std::to_string(i + 1) + "/tf_f2q_b", cof_b);
        double jnt_initial_pos = position_handler.get_latest_jnt_position()[i] - position_ref[i];
        DiscreteTimeFilter filter_(cof_a, cof_b, jnt_initial_pos);
        admittance_filters.push_back(filter_);
    }
    // openrave
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr penv = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
    boost::thread thviewer(boost::bind(SetViewer,penv,viewername));
    penv->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot;
    robot = penv->GetRobot(robotname);
    auto manip = robot->GetActiveManipulator();
    ros::Rate rate(125);
    // temp variable
    std::vector<double> position_cmd = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench = {0, 0, 0, 0, 0, 0};
    std::vector<double> tau_prj(6), tau_prj1(6), tau_prj2(6), force(3), torque(3);
    std::vector<double> jacobian, jacobian_rot;
    // controller
    while(!ros::isShuttingDown()){
        auto tstart = ros::Time::now();
        // call all callbacks
        ros::spinOnce();
        // retrieve the latest wrench measurement, offseted and filtered
        ft_handler.get_latest_wrench(force, torque);

        // project to joint torque space
        robot->SetActiveDOFValues(position_handler.get_latest_jnt_position());
        manip->CalculateJacobian(jacobian);
        manip->CalculateAngularVelocityJacobian(jacobian_rot);

        matrix_mult(jacobian, force, tau_prj1);
        matrix_mult(jacobian_rot, torque, tau_prj2);
        matrix_add(tau_prj1, tau_prj2, tau_prj);

        // run through joint position filter
        for (int i = 0; i < 6; ++i) {
            double y_ = admittance_filters[i].compute(tau_prj[i]);
            position_cmd[i] = position_ref[i] + y_;
            if (i==0){
                ROS_DEBUG_STREAM("j0: y_=" << y_);
            }
        }
        // send command
        position_act.set_joint_positions(position_cmd);
        // record time required
        auto tend = ros::Time::now();
        ros::Duration tdur = tend - tstart;
        // report, clean up then sleep
        ROS_DEBUG_STREAM("force: " << force[0] << ", " << force[1] << ", " << force[2]);
        ROS_DEBUG_STREAM("torque: " << torque[0] << ", " << torque[1] << ", " << torque[2]);
        ROS_DEBUG_STREAM("tau: " << tau_prj[0] << ","
                                 << tau_prj[1] << ","
                                 << tau_prj[2] << ","
                                 << tau_prj[3] << ","
                                 << tau_prj[4] << ","
                                 << tau_prj[5]);
        ROS_DEBUG_STREAM("comp time: " << tdur.toSec() * 1000 << "ms");
        rate.sleep();
    }

  return 0;
}


