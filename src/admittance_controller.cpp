#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Float64.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <infinite_interaction/infinite_interaction_lib.h>
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
    ros::init(argc, argv, "admittance_controller");
    ros::NodeHandle node_handle("~");

    // load parameters from ROS parameter server
    std::string robot_ns = "denso";
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
    bool viewer, debug;
    node_handle.param("viewer", viewer, true);
    node_handle.param("debug", debug, false);

    // TODO: incorporate these parameters into a launch file
    std::string robot_control_method = "direct";
    std::string robot_ip_addr = "192.168.0.21";

    /* Hardware handles
     *
     */
    // Robot controllers
    std::shared_ptr<HWHandle::AbstractRobotController> robot_handle;
    if (robot_control_method == "ros_control"){
        robot_handle = std::make_shared<HWHandle::JointPositionController> (robot_ns, node_handle);
    }
    else if (robot_control_method == "direct"){
       robot_handle = std::make_shared<HWHandle::RC8HWController> (robot_ip_addr, 1);
    }
    DoubleVector jnt_init = robot_handle->get_latest_jnt();

    // FT sensor data acquisition setup via FTsensor handler
    FTSensorHandle ft_handle(node_handle, ft_topic);
    if (debug){
        ft_handle.set_debug(node_handle);
    }

    // debugger which publishes data to topics
    std::string debug_ns = "/debugger";
    InfInteraction::TopicDebugger debugger(debug_ns, node_handle);
    int wId, yId, uId, qId, qcmdId;
    wId = debugger.register_multiarray("w_n"); // meeasured wrench
    yId = debugger.register_multiarray("y_n"); // force_ output
    uId = debugger.register_multiarray("u_n"); // position command
    qId = debugger.register_multiarray("q_n"); // joint command
    qcmdId = debugger.register_multiarray("qcmd_n"); // joint command

    // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr env_ptr = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
    if (viewer) boost::thread thviewer(boost::bind(SetViewer, env_ptr, viewername));  // create viewer
    env_ptr->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot_ptr;
    robot_ptr = env_ptr->GetRobot(robot_name);
    auto manip_ptr = robot_ptr->GetActiveManipulator();
    // set active dofs to [0,..,5]
    robot_ptr->SetActiveDOFs(std::vector<int> {0, 1, 2, 3, 4, 5});

    /* Control loop components construct.
     *
     * Interaction controller selection: depending on the loaded
     * parameter, difference "blocks" will be initialized. This is
     * where most of the work should be done: creating custom classes
     * for the task at hand and initializing them. The control loop is
     * very simple.
     */
    std::string controller_id, controller_type;
    if (!node_handle.getParam("/active_controller", controller_id)){
        ROS_ERROR_STREAM("No active controller found. Check if parameters have been loaded to parameter sever!");
        ros::shutdown();
    }
    if (!node_handle.getParam("/" + controller_id + "/type", controller_type)){
        ROS_ERROR_STREAM("Unable to decide controller type. Check if parameters have been loaded to parameter sever!");
        ros::shutdown();
    }
    // Initialize controlling blocks for the given task
    std::shared_ptr<SignalBlock> force_map_ptr;
    std::shared_ptr<SignalBlock> controller_ptr;
    std::shared_ptr<SignalBlock> position_map_ptr;
    if (controller_type == "joint_admittance"){
        // set reference wrench
        std::vector<double> wrench_offset;
        std::string ft_filter_path;
        if (!node_handle.getParam("/" + controller_id + "/ft_sensor_ref", wrench_offset)){
            ROS_WARN_STREAM("Reference wrench not found! Setting zero reference");
            wrench_offset.resize(6);
            std::fill(wrench_offset.begin(), wrench_offset.end(), 0);
        }
        if (wrench_offset.size() != 6){
            ROS_ERROR_STREAM("Reference wrench from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
	    ros::shutdown();
        }
        ft_handle.set_wrench_offset(wrench_offset);

        // reference position: the zero joint position. Required for this type of controller.
        std::vector<double> jnt_ref;
        node_handle.getParam("/" + controller_id + "/joint_pos_ref", jnt_ref);
        if (jnt_ref.size() != 6){
            ROS_ERROR_STREAM("Reference position from param sever is invalid! Have you loaded the parameters? \n -- Exitting!");
            ros::shutdown();
        }
        // mapping from external wrench to joint torque
        force_map_ptr = std::make_shared<InfInteraction::JointTorqueFromWrenchProjector>(robot_ptr, ft_name);
        // Init Joint Controllers, then throw 'em in a Controller Collection
        std::vector<std::shared_ptr<SignalBlock > > jnt_controllers;
        for (int i = 0; i < 6; ++i) {
            std::string jnt_filter_path, jnt_filter_type;
            node_handle.getParam("j" + std::to_string(i + 1) + "/filter", jnt_filter_path);
            jnt_filter_type = jnt_filter_path.substr(1, 8);
            if (jnt_filter_type == "iir_siso"){
                DoubleVector cof_a, cof_b;
                if (not node_handle.getParam(jnt_filter_path + "/b", cof_b)){
                    ROS_ERROR_STREAM("Unable to find filter: " << jnt_filter_path << ". Shutting ROS down");
                    ros::shutdown();
                }
                node_handle.getParam(jnt_filter_path + "/a", cof_a);
                jnt_controllers.push_back(std::make_shared<DiscreteTimeFilter>(cof_b, cof_a, jnt_init[i] - jnt_ref[i]));
            }
            else if (jnt_filter_type == "fir_siso") {
                int T, nu, ny;
                DoubleVector L, MB2;
                if (not node_handle.getParam(jnt_filter_path + "/T", T)){
                    ROS_ERROR_STREAM("Unable to find filter: " << jnt_filter_path << ". Shutting ROS down");
                    ros::shutdown();
                }
                node_handle.getParam(jnt_filter_path + "/nu", nu);
                node_handle.getParam(jnt_filter_path + "/ny", ny);
                node_handle.getParam(jnt_filter_path + "/L", L);
                node_handle.getParam(jnt_filter_path + "/MB2", MB2);
                jnt_controllers.push_back(std::make_shared<FIRsrfb>(T, ny, nu, L, MB2, DoubleVector{jnt_init[i] - jnt_ref[i]}));
            }
            else {
                throw std::invalid_argument("Unknown filter kind");
            }
        }
        controller_ptr = std::make_shared<InfInteraction::ControllerCollection>(jnt_controllers);
        // inverse kinematic: simply add the reference joint position to each position output to obtain the joint command.
        position_map_ptr = std::make_shared<InfInteraction::SimpleOffset>(jnt_ref);
    }
    else if (controller_type == "cartesian_3D_admittance"){
        // get current wrench reading and set it as the offset
        ros::Duration(0.3).sleep(); ros::spinOnce();  // make sure to receive at least a wrench reading before continue
        DoubleVector wrench_offset;
        ft_handle.get_latest_wrench(wrench_offset);
        ft_handle.set_wrench_offset(wrench_offset);
        // from external wrench to joint torque
        robot_ptr->SetActiveDOFValues(jnt_init); // set the robot's initial configuraiton before initializing force_ projection block
        force_map_ptr = std::make_shared<InfInteraction::Wrench2CartForceProjector>(robot_ptr, ft_name);
        // Cartesian admittance controller is a fir_mino type. The initial position (Cartesian, 3D) is zero.
        // The initial pose is kept fixed.
        std::string filter_path;
        int T, nu, ny;
        DoubleVector L, MB2;
        if (!node_handle.getParam("/" + controller_id + "/filter", filter_path)){
            ROS_FATAL_STREAM("Does not find any filter in controller id: " << controller_id << ". Shutting down!");
            ros::shutdown();
        }
        if (not node_handle.getParam(filter_path + "/T", T)){
            ROS_ERROR_STREAM("Unable to find filter: " << filter_path << ". Shutting down");
            ros::shutdown();
        }
        node_handle.getParam(filter_path + "/nu", nu);
        node_handle.getParam(filter_path + "/ny", ny);
        node_handle.getParam(filter_path + "/L", L);
        node_handle.getParam(filter_path + "/MB2", MB2);
        controller_ptr = std::make_shared<FIRsrfb>(T, ny, nu, L, MB2);
        // Inverse Kinematics: Converts Cartesian signals to joint coordinates
        double gam, gam2;  // weights for regulation of dq, and of (q_current - q_init)
        node_handle.param<double>("/" + controller_id + "/gam", gam, 0.01);
        node_handle.param<double>("/" + controller_id + "/gam2", gam2, 0.001);
        position_map_ptr = std::make_shared<InfInteraction::CartPositionTracker>(robot_ptr, manip_name, jnt_init, gam, gam2);
    }
     else if (controller_type == "cartesian_3D_admittance_Qparam"){
        // get current wrench reading and set it as the offset
        ros::Duration(0.3).sleep(); ros::spinOnce();  // make sure to receive at least a wrench reading before continue
        DoubleVector wrench_offset;
        ft_handle.get_latest_wrench(wrench_offset);
        ft_handle.set_wrench_offset(wrench_offset);
        // from external wrench_measure to joint torque
        robot_ptr->SetActiveDOFValues(jnt_init); // set the robot's initial configuraiton before initializing force_ projection block
        force_map_ptr = std::make_shared<InfInteraction::Wrench2CartForceProjector>(robot_ptr, ft_name);
        // The initial pose is kept fixed.
        std::string filter_path;
        if (!node_handle.getParam("/" + controller_id + "/filter", filter_path)){
            ROS_FATAL_STREAM("Unable to find a filter in controller: " << controller_id << ". Shutting down!");
            ros::shutdown();
        }
        DoubleVector fb_b, fb_a, ff_taps;
        if (not node_handle.getParam(filter_path + "/xff_taps", ff_taps)){
            ROS_ERROR_STREAM("Unable to find filter: " << filter_path << ". Shutting down");
            ros::shutdown();
        }
        node_handle.getParam(filter_path + "/xff_taps", ff_taps);
        node_handle.getParam(filter_path + "/xfb_b", fb_b);
        node_handle.getParam(filter_path + "/xfb_a", fb_a);
        auto xff = std::make_shared<DiscreteTimeFilter> (ff_taps);
        auto xfb = std::make_shared<DiscreteTimeFilter> (fb_b, fb_a);
        auto xfilter = std::make_shared<InfInteraction::DelayFeedback>(xff, xfb, -1);

        node_handle.getParam(filter_path + "/yff_taps", ff_taps);
        node_handle.getParam(filter_path + "/yfb_b", fb_b);
        node_handle.getParam(filter_path + "/yfb_a", fb_a);
        auto yff = std::make_shared<DiscreteTimeFilter> (ff_taps);
        auto yfb = std::make_shared<DiscreteTimeFilter> (fb_b, fb_a);
        auto yfilter = std::make_shared<InfInteraction::DelayFeedback>(yff, yfb, -1);

        node_handle.getParam(filter_path + "/zff_taps", ff_taps);
        node_handle.getParam(filter_path + "/zfb_b", fb_b);
        node_handle.getParam(filter_path + "/zfb_a", fb_a);
        auto zff = std::make_shared<DiscreteTimeFilter> (ff_taps);
        auto zfb = std::make_shared<DiscreteTimeFilter> (fb_b, fb_a);
        auto zfilter = std::make_shared<InfInteraction::DelayFeedback>(zff, zfb, -1);

        std::vector<std::shared_ptr<SignalBlock > > ctrls;
        ctrls.push_back(xfilter);
        ctrls.push_back(yfilter);
        ctrls.push_back(zfilter);

        controller_ptr = std::make_shared<InfInteraction::ControllerCollection>(ctrls);

        // Inverse Kinematics: Converts Cartesian signals to joint coordinates
        double gam, gam2;  // weights for regulation of dq, and of (q_current - q_init)
        node_handle.param<double>("/" + controller_id + "/gam", gam, 0.01);
        node_handle.param<double>("/" + controller_id + "/gam2", gam2, 0.001);
        position_map_ptr = std::make_shared<InfInteraction::CartPositionTracker>(robot_ptr, manip_name, jnt_init, gam, gam2);
    }
    else {
        ROS_ERROR_STREAM("controller type not recognized! Quitting.");
        ros::shutdown();
    }

    /* Real-Time (RT) control loop.
     *
     * The program enters an infinite control loop. In this, it takes
     * wrench measurements, transform it to joint positions through an
     * intricate set of operations before sending the values to the
     * robot.
     */
    // pre-defined loop variables
    std::vector<double> joint_cmd = {0, 0, 0, 0, 0, 0}, joint_measure(6);
    std::vector<double> wrench_measure = {0, 0, 0, 0, 0, 0}, y_n, u_n;

    // rt setups
    timespec slp_dline_spec, sent_spec, wake_spec;
    if (RTUtils::set_policy_fifo() == 0) {
      ROS_INFO("sch policy set to fifo");
    }
    else{
      ROS_INFO("sch policy remained default");
    }

    const int FOUR_MS = 4000000;
    int diff_nsec; // time difference between (i) the deadline to wake up and (ii) the time when the message is received.
    clock_gettime(CLOCK_MONOTONIC, &slp_dline_spec);
    int cy_idx = 0;  // cycle index
    // Control loop
    while(!ros::isShuttingDown()){
        // call all callbacks
        ros::spinOnce();

        // collect current wrench and joint position
        ft_handle.get_latest_wrench(wrench_measure);
        joint_measure = robot_handle->get_latest_jnt();

        // wrench transform: convert measured wrench to control output (Cartesian force_)
        force_map_ptr->set_state(joint_measure);
        y_n = force_map_ptr->compute(wrench_measure);

        // control law: run control output through controller
        u_n = controller_ptr->compute(y_n);

        // inverse kinematics: map control displacement to actual joint command
        position_map_ptr->set_state(joint_measure);
        joint_cmd = position_map_ptr->compute(u_n);

        // stage 1 finished; sleep
        RTUtils::increment_timespec(slp_dline_spec, FOUR_MS);
        clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &slp_dline_spec, NULL); // sleep until slp_dline_spec

        // send joint position command
        clock_gettime(CLOCK_MONOTONIC, &wake_spec);
        robot_handle-> send_jnt_command(joint_cmd);
        clock_gettime(CLOCK_MONOTONIC, &sent_spec);

        // publish debug information
        if (debug) {
            debugger.publish_multiarray(wId, wrench_measure);
            debugger.publish_multiarray(yId, y_n);
            debugger.publish_multiarray(uId, u_n);
            debugger.publish_multiarray(qId, joint_measure);
            debugger.publish_multiarray(qcmdId, joint_cmd);

            // publish time discrepancy
	    if (cy_idx % 125 == 0){
	      RTUtils::diff_timespec(diff_nsec, slp_dline_spec, wake_spec);
	      ROS_DEBUG_STREAM("wake - deadline-to-send: " << diff_nsec << " nsec (this value should be very small, ideally less than 1e4)");
	      RTUtils::diff_timespec(diff_nsec, slp_dline_spec, sent_spec);
	      ROS_DEBUG_STREAM("sent - deadline-to-send: " << diff_nsec << " nsec (this value should be small, ideally less than 1e5)");
	    }
        }
        cy_idx = (cy_idx + 1) % 2147483640;

        // stage 2 finsished, sleep
        RTUtils::increment_timespec(slp_dline_spec, FOUR_MS);
        clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &slp_dline_spec, NULL);
    }

    return 0;
}


