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
    ros::init(argc, argv, "sin_motion_generator");
    ros::NodeHandle node_handle("~");
    // parameter
    std::string name_space = "denso";
    std::string ft_sensor_topic = "/netft/raw";
    std::string joint_state_topic = "/denso/joint_states";
    std::string scenefilename = "robots/denso_handle.robot.xml";
    std::string viewername = "qtosg";
    std::string robotname = "denso_handle";
    std::string ftmanipname = "FTsensor";
    std::vector<double> wrench_ref;
    node_handle.getParam("ft_sensor_ref", wrench_ref);
    std::vector<double> position_ref;
    node_handle.getParam("joint_pos_ref", position_ref);
    double accel_duration, reach_duration;
    node_handle.getParam("accel_duration", accel_duration);
    if (!node_handle.getParam("reach_duration", reach_duration)){
        ROS_ERROR_STREAM("No reach duration found! Use default value (5 seconds)");
        reach_duration = 5;
    };
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
    // torque inputs
    ExternalTorquePublisher torque_pub(name_space, node_handle);

    // openrave
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr penv = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);

    penv->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot;
    robot = penv->GetRobot(robotname);
    robot->SetActiveManipulator(ftmanipname);
    auto manip = robot->GetActiveManipulator();
    int cmd_rate = 125;
    ros::Rate rate(cmd_rate);
    // temp variable
    std::vector<double> position_cmd = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench = {0, 0, 0, 0, 0, 0};
    std::vector<double> tau_prj(6), tau_prj1(6), tau_prj2(6), force(3), torque(3);
    std::vector<double> jacobian, jacobian_rot;
    OpenRAVE::Transform T_wee;
    // controller
    bool done = false;
    while(!ros::isShuttingDown() and !done){
        // ask for motion parameters
        double t_duration, omega, amplitude;
        int jnt_idx;
        cout << "Enter [joint index], [duration], [angular velocity] and amplitude. Separated by space." << std::endl;
        std::cin >> jnt_idx >> t_duration >> omega >> amplitude;
        // move the robot to the reference position in 5 seconds
        ros::spinOnce();
        dVector position_init = position_handler.get_latest_jnt_position();
        ROS_INFO_STREAM("Starting experiment. Moving to reference position!. Will take " << reach_duration << " seconds.");
        for (int j = 0; j < (reach_duration * cmd_rate + 1); ++j) {
            double alpha = (double)(j) / (reach_duration * cmd_rate);
            for (int i = 0; i < 6; ++i) {
                position_cmd[i] = position_init[i] * (1 - alpha) + position_ref[i] * alpha;
            }
            position_act.set_joint_positions(position_cmd);
            rate.sleep();
        }

        // execute sinuisoidal motion
        double t_current = 0;
        double amplitude_scaled = 0;
        auto tstart = ros::Time::now();
        ros::Duration(1).sleep();
        ROS_INFO_STREAM("Start generating sinuisoidal motion.");
        while (!ros::isShuttingDown() and t_current < t_duration){
            auto tstart_loop = ros::Time::now();
            // call all callbacks
            ros::spinOnce();
            // retrieve the latest wrench measurement, offseted and filtered
            ft_handler.get_latest_wrench(force, torque);
            OpenRAVE::RaveVector<double> rave_force(force[0], force[1], force[2]);
            OpenRAVE::RaveVector<double> rave_torque(torque[0], torque[1], torque[2]);

            // project to joint torque space
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

            // sinuisoidal trajectory. in the initial 2 seconds the amplitude linearly increases.
            if (t_current > accel_duration){
                amplitude_scaled = amplitude;
            }
            else {
                amplitude_scaled = t_current / accel_duration * amplitude;
            }
            position_cmd[jnt_idx] = position_ref[jnt_idx] + amplitude_scaled * sin(omega * t_current);

            // send command
            position_act.set_joint_positions(position_cmd);
            torque_pub.publish_joint_torques(tau_prj);

            // record time required
            auto tend = ros::Time::now();
            ros::Duration tdur = tend - tstart_loop;
            // report, clean up then sleep
            ROS_DEBUG_STREAM_THROTTLE(1, "force: " << force[0] << ", " << force[1] << ", " << force[2]);
            ROS_DEBUG_STREAM_THROTTLE(1, "torque: " << torque[0] << ", " << torque[1] << ", " << torque[2]);
            ROS_DEBUG_STREAM_THROTTLE(1, "comp time: " << tdur.toSec() * 1000 << "ms");
            t_current = (tend - tstart).toSec();
            rate.sleep();
        }

        std::cout << "Continue? Any letter to continue, n to stop.";
        std::string cont;
        std::cin >> cont;
        if (cont == "n"){
            std::cout << "Quitting!" << std::endl;
            break;
        }
    }

    return 0;
}

