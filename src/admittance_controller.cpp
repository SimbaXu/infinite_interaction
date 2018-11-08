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

void print_dVector(dVector x, std::string name){
    std::cout << name << ": ";
    for(int i=0; i < x.size(); ++i){
        std::cout << x[i] << ", ";
    }
    std::cout << std::endl;
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "joint_admittance_controller");
    ros::NodeHandle node_handle("~");

    // general parameters
    std::string name_space = "denso";
    std::string ft_topic = "/netft/raw";
    std::string jnt_state_topic = "/denso/joint_states";
    std::string viewername = "qtosg";
    std::string scenefilename, robot_name, manip_name, ft_name;
    if (!node_handle.getParam("/scene_path", scenefilename)){
        ROS_ERROR_STREAM("Unable to find [/scene_path]. Have you loaded all parameters?");
        ros::shutdown();
    }
    if (!node_handle.getParam("/robot_name", robot_name)){
        ROS_ERROR_STREAM("Unable to find [/robot_name]. Have you loaded all parameters?");
        ros::shutdown();
    }
    if (!node_handle.getParam("/manip_frame", manip_name)){
        ROS_ERROR_STREAM("Unable to find [/manip_name]. Have you loaded all parameters?");
        ros::shutdown();
    }
     if (!node_handle.getParam("/ftsensor_frame", ft_name)){
        ROS_ERROR_STREAM("Unable to find [/ftsensor_name]. Have you loaded all parameters?");
        ros::shutdown();
    }
    ros::Rate rate(125);
    bool viewer;
    node_handle.param("/viewer", viewer, true);

    // FT sensor data acquisition setup via FTsensor handler
    // - define wrench filtering and
    // - offseting
    std::vector<double> wrench_ref, FT_filter_b, FT_filter_a;
    node_handle.getParam("ft_sensor_ref", wrench_ref);
    if (wrench_ref.size() != 6){
        ROS_ERROR_STREAM("Reference wrench from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    std::string ft_filter_path;
    node_handle.getParam("ft_sensor_filter", ft_filter_path);
    node_handle.getParam(ft_filter_path + "/b", FT_filter_b);
    node_handle.getParam(ft_filter_path + "/a", FT_filter_a);
    FTSensorHandler ft_handler(wrench_ref, FT_filter_b, FT_filter_a);
    ros::Subscriber ft_subscriber = node_handle.subscribe(ft_topic, 3, &FTSensorHandler::signal_callback, &ft_handler);
    bool FT_debug;
    node_handle.param<bool>("debug", FT_debug, false);
    if (FT_debug){
        ft_handler.set_debug(node_handle);
    }

    // Robot Joint Position data acquisition (from ros control) via JointPositionHandler
    std::vector<double> jnt_pos_ref;
    node_handle.getParam("joint_pos_ref", jnt_pos_ref);
    if (jnt_pos_ref.size() != 6){
        ROS_ERROR_STREAM("Reference position from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
        exit(0);
    }
    JointPositionHandler jnt_pos_handler;
    ros::Subscriber jnt_pos_subscriber = node_handle.subscribe(jnt_state_topic, 3, &JointPositionHandler::signal_callback, &jnt_pos_handler);
    ros::Duration(0.5).sleep(); ros::spinOnce(); // wait for a few second before updating the current joint position
    if (!jnt_pos_handler.received_msg()){
        ROS_FATAL("Have not received messages to update initial joint position! \n-- Terminating!");
        exit(0);
    }

    // Robot Joint Position commanding interface
    JointPositionController jnt_pos_act(name_space, node_handle);
    dVector jnt_pos_init = jnt_pos_handler.get_latest_jnt_position();
    ROS_INFO_STREAM("Initial position: " << jnt_pos_init[0] << ", " << jnt_pos_init[1] << ", "
		    << jnt_pos_init[2] << ", " << jnt_pos_init[3] << ", "
		    << jnt_pos_init[4] << ", " << jnt_pos_init[5]);

    // Publisher for torque inputs to the controller (for logging/debugging purposes)
    ExternalTorquePublisher torque_pub(name_space, node_handle);

    // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr env_ptr = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
    if (viewer) boost::thread thviewer(boost::bind(SetViewer,env_ptr,viewername));  // create viewer
    env_ptr->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot_ptr;
    robot_ptr = env_ptr->GetRobot(robot_name);
    robot_ptr->SetActiveManipulator(ft_name);
    robot_ptr->SetActiveDOFs(std::vector<int> {0, 1, 2, 3, 4, 5});
    auto manip_ptr = robot_ptr->GetActiveManipulator();

    // Interaction controller selection: depending on the loaded parameter, difference "blocks" will
    // be initialized. This is where most of the work should be done: creating custom classes for the
    // task at hand and initializing them. The control loop is actually trivial.
    std::string controller_id, controller_type;
    if (!node_handle.getParam("/active_controller", controller_id)){
        ROS_ERROR_STREAM("No active controller found. Check if parameters have been loaded to parameter sever!");
        ros::shutdown();
    }
    if (!node_handle.getParam("/" + controller_id + "/type", controller_type)){
        ROS_ERROR_STREAM("Unable to decide controller type. Check if parameters have been loaded to parameter sever!");
        ros::shutdown();
    }

    std::shared_ptr<LTI> force_map_ptr;
    std::shared_ptr<LTI> controller_ptr;
    std::shared_ptr<LTI> position_map_ptr;
    // Initialize the controllers for a given task
    if (controller_type == "joint_admittance"){
        force_map_ptr = std::make_shared<InfInteraction::JointTorqueFromWrenchProjector>(robot_ptr, ft_name);
        // Init Joint Controllers, then throw 'em in a Controller Collection
        std::vector<std::shared_ptr<LTI > > jnt_controllers;
        for (int i = 0; i < 6; ++i) {
            std::string jnt_filter_path, jnt_filter_type;
            node_handle.getParam("j" + std::to_string(i + 1) + "/filter", jnt_filter_path);
            jnt_filter_type = jnt_filter_path.substr(1, 8);
            if (jnt_filter_type == "iir_siso"){
                dVector cof_a, cof_b;
                if (not node_handle.getParam(jnt_filter_path + "/b", cof_b)){
                    ROS_ERROR_STREAM("Unable to find filter: " << jnt_filter_path << ". Shutting ROS down");
                    ros::shutdown();
                }
                node_handle.getParam(jnt_filter_path + "/a", cof_a);
                jnt_controllers.push_back(std::make_shared<DiscreteTimeFilter>(cof_b, cof_a, jnt_pos_init[i] - jnt_pos_ref[i]));
            }
            else if (jnt_filter_type == "fir_siso") {
                int T, nu, ny;
                dVector L, MB2;
                if (not node_handle.getParam(jnt_filter_path + "/T", T)){
                    ROS_ERROR_STREAM("Unable to find filter: " << jnt_filter_path << ". Shutting ROS down");
                    ros::shutdown();
                }
                node_handle.getParam(jnt_filter_path + "/nu", nu);
                node_handle.getParam(jnt_filter_path + "/ny", ny);
                node_handle.getParam(jnt_filter_path + "/L", L);
                node_handle.getParam(jnt_filter_path + "/MB2", MB2);
                jnt_controllers.push_back(std::make_shared<FIRsrfb>(T, ny, nu, L, MB2, dVector{jnt_pos_init[i] - jnt_pos_ref[i]}));
            }
            else {
                throw std::invalid_argument("Unknown filter kind");
            }
        }
        controller_ptr = std::make_shared<InfInteraction::ControllerCollection>(jnt_controllers);
        // inverse kinematic: simply add the reference joint position to each position output to obtain the joint command.
        position_map_ptr = std::make_shared<InfInteraction::SimpleOffset>(jnt_pos_ref);
    }

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

        // transform offset wrench to joint torque
        OpenRAVE::RaveVector<double> rave_force(force[0], force[1], force[2]);
        OpenRAVE::RaveVector<double> rave_torque(torque[0], torque[1], torque[2]);
        robot_ptr->SetActiveDOFValues(jnt_pos_handler.get_latest_jnt_position());
        manip_ptr->CalculateJacobian(jacobian);
        manip_ptr->CalculateAngularVelocityJacobian(jacobian_rot);
        T_wee = manip_ptr->GetEndEffectorTransform();
        auto rave_force_ = T_wee.rotate(rave_force);
        auto rave_torque_ = T_wee.rotate(rave_torque);
        force[0] = rave_force_.x; force[1] = rave_force_.y; force[2] = rave_force_.z;
        torque[0] = rave_torque_.x; torque[1] = rave_torque_.y; torque[2] = rave_torque_.z;

        auto jacobian_T = mat_transpose(jacobian, 6);
        auto jacobian_rot_T = mat_transpose(jacobian_rot, 6);
        matrix_mult(jacobian_T, force, tau_prj1);
        matrix_mult(jacobian_rot_T, torque, tau_prj2);
        matrix_add(tau_prj1, tau_prj2, tau_prj);

        // new pipeline
        ft_handler.get_latest_wrench(wrench);
        force_map_ptr->set_state(jnt_pos_handler.get_latest_jnt_position());
        auto tau0 = force_map_ptr->compute(wrench);
        auto new_un = controller_ptr->compute(tau0);
        position_cmd = position_map_ptr->compute(new_un);

        // Send joint position command
        jnt_pos_act.set_joint_positions(position_cmd);
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


