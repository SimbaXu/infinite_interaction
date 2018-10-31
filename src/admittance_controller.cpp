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

// TODO: Change this node name into something more descriptive like
// joint admittance controller


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
    ros::init(argc, argv, "joint_admittance_controller");
    ros::NodeHandle node_handle("~");
    // general parameters
    std::string name_space = "denso";
    std::string ft_sensor_topic = "/netft/raw";
    std::string joint_state_topic = "/denso/joint_states";
    std::string scenefilename = "robots/denso_handle.robot.xml";
    std::string viewername = "qtosg";
    std::string robotname = "denso_handle";
    std::string ftmanipname = "FTsensor";

    // FT sensor data acquisition setup via FTsensor handler
    // - define wrench filtering and
    // - offseting
    std::vector<double> wrench_ref, FT_filter_b, FT_filter_a;
    node_handle.getParam("ft_sensor_ref", wrench_ref);
    if (wrench_ref.size() != 6){
        ROS_ERROR_STREAM("Reference wrench from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    node_handle.getParam("ft_sensor_filters/b", FT_filter_b);
    node_handle.getParam("ft_sensor_filters/a", FT_filter_a);
    FTSensorHandler ft_handler(wrench_ref, FT_filter_b, FT_filter_a);
    ros::Subscriber ft_sub = node_handle.subscribe(ft_sensor_topic, 3, &FTSensorHandler::signal_callback, &ft_handler);
    bool FT_debug;
    node_handle.param<bool>("debug", FT_debug, false);
    if (FT_debug){
        ft_handler.set_debug(node_handle);
    }
    // Robot Joint Position data acquisition (from ros control) via JointPositionHandler
    std::vector<double> position_ref;
    node_handle.getParam("joint_pos_ref", position_ref);
    if (position_ref.size() != 6){
        ROS_ERROR_STREAM("Reference position from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    JointPositionHandler position_handler;
    ros::Subscriber position_handler_sub = node_handle.subscribe(joint_state_topic, 3, &JointPositionHandler::signal_callback, &position_handler);
    ros::Duration(0.5).sleep(); ros::spinOnce(); // wait for a few second before updating the current joint position
    if (!position_handler.received_msg()){
        ROS_FATAL("Have not received messages to update initial joint position! \n-- Terminating!");
        exit(0);
    }

    // Robot Joint Position commanding interface
    JointPositionController position_act(name_space, node_handle);
    dVector initial_jnt_position = position_handler.get_latest_jnt_position();
    ROS_INFO_STREAM("Initial position: " << initial_jnt_position[0] << ", " << initial_jnt_position[1] << ", "
		    << initial_jnt_position[2] << ", " << initial_jnt_position[3] << ", "
		    << initial_jnt_position[4] << ", " << initial_jnt_position[5]);

    // Individual Joint Controllers
    std::vector<DiscreteTimeFilter> admittance_filters;
    for (int i = 0; i < 6; ++i) {
        dVector cof_a, cof_b;
        node_handle.getParam("j" + std::to_string(i + 1) + "/tf_f2q_a", cof_a);
        node_handle.getParam("j" + std::to_string(i + 1) + "/tf_f2q_b", cof_b);
        double jnt_initial_pos = position_handler.get_latest_jnt_position()[i] - position_ref[i];
        DiscreteTimeFilter filter_(cof_a, cof_b, jnt_initial_pos);
        admittance_filters.push_back(filter_);
    }

    // Publisher for torque inputs to the controller (for logging/debugging purposes)
    ExternalTorquePublisher torque_pub(name_space, node_handle);

    // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr penv = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
    boost::thread thviewer(boost::bind(SetViewer,penv,viewername));
    penv->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot;
    robot = penv->GetRobot(robotname);
    robot->SetActiveManipulator(ftmanipname);
    auto manip = robot->GetActiveManipulator();
    ros::Rate rate(125);

    // temp variable
    std::vector<double> position_cmd = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench = {0, 0, 0, 0, 0, 0};
    std::vector<double> tau_prj(6), tau_prj1(6), tau_prj2(6), force(3), torque(3);
    std::vector<double> jacobian, jacobian_rot;
    OpenRAVE::Transform T_wee;
    int step_idx = 0;

    // Control loop
    while(!ros::isShuttingDown()){
        auto tstart = ros::Time::now();

        // call all callbacks
        ros::spinOnce();

        // retrieve the latest wrench measurement, offseted and filtered
        ft_handler.get_latest_wrench(force, torque);
        OpenRAVE::RaveVector<double> rave_force(force[0], force[1], force[2]);
        OpenRAVE::RaveVector<double> rave_torque(torque[0], torque[1], torque[2]);

        // transform offset wrench to joint torque
        robot->SetActiveDOFValues(position_handler.get_latest_jnt_position());
        manip->CalculateJacobian(jacobian);
        manip->CalculateAngularVelocityJacobian(jacobian_rot);
        T_wee = manip->GetEndEffectorTransform();
        auto rave_force_ = T_wee.rotate(rave_force);
        auto rave_torque_ = T_wee.rotate(rave_torque);
        force[0] = rave_force_.x; force[1] = rave_force_.y; force[2] = rave_force_.z;
        torque[0] = rave_torque_.x; torque[1] = rave_torque_.y; torque[2] = rave_torque_.z;

        auto jacobian_T = mat_transpose(jacobian, 6);
        auto jacobian_rot_T = mat_transpose(jacobian_rot, 6);
        matrix_mult(jacobian_T, force, tau_prj1);
        matrix_mult(jacobian_rot_T, torque, tau_prj2);
        matrix_add(tau_prj1, tau_prj2, tau_prj);

        // Compute joint actuation commands from projected joint torque
        for (int i = 0; i < 6; ++i) {
            double y_ = admittance_filters[i].compute(tau_prj[i]);
            position_cmd[i] = position_ref[i] + y_;
        }

        // Send joint position command
        position_act.set_joint_positions(position_cmd);
        torque_pub.publish_joint_torques(tau_prj);

        // record time required
        auto tend = ros::Time::now();
        ros::Duration tdur = tend - tstart;

        // report, clean up then sleep
        ROS_DEBUG_STREAM_THROTTLE(1, "force: " << force[0] << ", " << force[1] << ", " << force[2]);
        ROS_DEBUG_STREAM_THROTTLE(1, "torque: " << torque[0] << ", " << torque[1] << ", " << torque[2]);
        ROS_DEBUG_STREAM_THROTTLE(1, "comp time: " << tdur.toSec() * 1000 << "ms");
        // log data
//        ROS_DEBUG_STREAM("tau: " << step_idx << ","
//                                 << tau_prj[0] << "," << tau_prj[1] << "," << tau_prj[2] << ","
//                                 << tau_prj[2] << "," << tau_prj[4] << "," << tau_prj[5]);
//        ROS_DEBUG_STREAM("cmd: " << step_idx << ","
//                                 << position_cmd[0] << "," << position_cmd[1] << "," << position_cmd[2] << ","
//                                 << position_cmd[2] << "," << position_cmd[4] << "," << position_cmd[5]);
        step_idx = (step_idx + 1) % 2147483640;
        rate.sleep();
    }

  return 0;
}


