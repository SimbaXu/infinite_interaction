#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Float64.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <include/infinite_interaction/infinite_interaction_lib.h>
// openrave
#include <openrave-core.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>


using std::cout;


// Generate a vector of reference command
// reference motion:
// mo_type := 1,  x(t) = min(M, dM * t) * sin(omega * t)
// mo_type := 2,  x(t+1) = M * sin(phase[t]);  phase[t + 1] = phase[t] + omega / t_duration * t
// motion type 1 generates motion with a single frequency
// motion type 2 is used to generate smooth motion that have varying
// frequencies.
dVector generate_reference_command(int mo_type, double total_duration, double omega, double amplitude){
    dVector jnt_pos_vec(125 * 3);
    // any motion generated using this function has a mandatory 3 seconds rests
    std::fill(jnt_pos_vec.begin(), jnt_pos_vec.end(), 0);
    if (mo_type == 1) {
        double amplitude_scaled, t_current, accel_duration, jnt_pos_cur=0;
        accel_duration = 2;
        // sinuisoidal trajectory. in the initial 2 seconds the amplitude linearly increases.
        t_current = 0;
        while(t_current < total_duration){
            if (t_current > accel_duration){
                amplitude_scaled = amplitude;
            }
            else {
                amplitude_scaled = t_current / accel_duration * amplitude;
            }
            jnt_pos_cur = amplitude_scaled * sin(omega * t_current);
            jnt_pos_vec.push_back(jnt_pos_cur);
            t_current += 0.008;
        }
    }
    else if (mo_type == 2){
        double t_current = 0, jnt_pos_cur, phase=0, omega_rate = omega / (0.8 * total_duration);
        while(t_current < total_duration){
            if (t_current < 0.8 * total_duration){
                phase += omega_rate * t_current * 0.008;
            }
            else {
                // swift slow down
                phase += (omega - 4 * omega_rate * (t_current - 0.8 * total_duration)) * 0.008;
            }
            jnt_pos_cur = amplitude * sin(phase);
            jnt_pos_vec.push_back(jnt_pos_cur);
            t_current += 0.008;
        }
    }
    return jnt_pos_vec;
}


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
    // get parameters
    std::vector<double> wrench_ref;
    if (not node_handle.getParam("ft_sensor_ref", wrench_ref)){
        ROS_ERROR_STREAM("Unable to get reference FT wrench from param sever! Have you loaded experiment_params.yaml? \n -- Exitting!");
        exit(0);
    }
    std::vector<double> jnt_pos_ref;
    node_handle.getParam("joint_pos_ref", jnt_pos_ref);
    double accel_duration, reach_duration;
    node_handle.getParam("accel_duration", accel_duration);
    if (!node_handle.getParam("reach_duration", reach_duration)){
        ROS_ERROR_STREAM("No [reach_duration] found on parameter server! Use default value (5 seconds)");
        reach_duration = 5;
    };
    int cmd_rate;
    if (!node_handle.getParam("cmd_rate", cmd_rate)){
        ROS_ERROR_STREAM("No [cmd_rate] found on parameter server! Use default value (125 Hz)");
        cmd_rate = 125;
    };

    if (jnt_pos_ref.size() != 6){
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
    dVector initial_jnt_position = position_handler.get_latest_jnt();
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
    ros::Rate rate(cmd_rate);
    // temp variable
    std::vector<double> jnt_pos_cmd = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench = {0, 0, 0, 0, 0, 0};
    std::vector<double> tau_prj(6), tau_prj1(6), tau_prj2(6), force(3), torque(3);
    std::vector<double> jacobian, jacobian_rot;
    OpenRAVE::Transform T_wee;
    // controller
    bool done = false;
    double t_duration, omega, amplitude;
    int jnt_idx, mo_type;
    while(!ros::isShuttingDown() and !done){
        // ask for motion parameters
        cout << "Enter [motion_type: 0 or 1] [joint_index: starting from 0], [duration], [angular_velocity] and [amplitude]. Separated by space."
                << std::endl << "e.g. [1 2 20 24 0.1]";
        std::cin >> mo_type >> jnt_idx >> t_duration >> omega >> amplitude;


        // generate commanding signal
        dVector jnt_pos_i = generate_reference_command(mo_type, t_duration, omega, amplitude);

        // move the robot to the reference position in [reach_duration] second
        ros::spinOnce();
        dVector position_init = position_handler.get_latest_jnt();
        ROS_INFO_STREAM("Starting experiment. Moving to reference position!. Will take " << reach_duration << " seconds.");
        for (int j = 0; j < (reach_duration * cmd_rate + 1); ++j) {
            double alpha = (double)(j) / (reach_duration * cmd_rate);
            for (int i = 0; i < 6; ++i) {
                jnt_pos_cmd[i] = position_init[i] * (1 - alpha) + jnt_pos_ref[i] * alpha;
            }
            position_act.send_jnt_command(jnt_pos_cmd);
            rate.sleep();
        }

        // execute sinuisoidal motion
        int n = 0;
        auto tstart = ros::Time::now();
        ros::Duration(1).sleep();
        ROS_INFO_STREAM("Start generating sinuisoidal motion.");
        while (!ros::isShuttingDown() and n < jnt_pos_i.size()){
            auto tstart_loop = ros::Time::now();
            // call all callbacks
            ros::spinOnce();
            // retrieve the latest wrench measurement, offseted and filtered
            ft_handler.get_latest_wrench(force, torque);
            OpenRAVE::RaveVector<double> rave_force(force[0], force[1], force[2]);
            OpenRAVE::RaveVector<double> rave_torque(torque[0], torque[1], torque[2]);

            // project to joint torque space
            robot->SetActiveDOFValues(position_handler.get_latest_jnt());
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

            jnt_pos_cmd[jnt_idx] = jnt_pos_ref[jnt_idx] + jnt_pos_i[n++];

            // send command
            position_act.send_jnt_command(jnt_pos_cmd);
            torque_pub.publish_joint_torques(tau_prj);

            // record time required
            auto tend = ros::Time::now();
            ros::Duration tdur = tend - tstart_loop;
            // report, clean up then sleep
            ROS_DEBUG_STREAM_THROTTLE(1, "force: " << force[0] << ", " << force[1] << ", " << force[2]);
            ROS_DEBUG_STREAM_THROTTLE(1, "torque: " << torque[0] << ", " << torque[1] << ", " << torque[2]);
            ROS_DEBUG_STREAM_THROTTLE(1, "comp time: " << tdur.toSec() * 1000 << "ms");
            rate.sleep();
        }

        std::cout << "Continue? Any letter to continue, n to stop.";
        std::string cont;
        std::cin >> cont;
        if (cont == "n"){
            std::cout << "Quitting!" << std::endl;
            ros::shutdown();
        }
    }

    return 0;
}

