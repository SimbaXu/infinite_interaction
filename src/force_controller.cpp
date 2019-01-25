#include <ros/ros.h>
#include <infinite_interaction/infinite_interaction_lib.h>

// openrave
#include <openrave-core.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

const double MM = 1e-3;
const double Ts = 0.008;
inline double ABS(double x) {return (x > 0) ? x: -x;}

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


/*! A combined class computes controls from sensory data.
 *
 * ControllerManager can
 * - manages multiple Z-axis force controllers,
 * - load them from ROS parameter server,
 * - switch between them,
 * - know how to execute the feedback structure of the Q controller, and
 * - scale output signal accordingly.
 *
 * NOTE: Currently this class is only used to handle one controller.
 */
class SISOForceControllerManager {
    // A vector of controllers, each is a vector of pointers to filters. If there are two filters per vector,
    // it is a PI. If there are four filters per vector, it is a Q-synthesis controller. Sequence: Q0, Q1, zPyu0, zPyu1
    std::vector<std::vector<std::shared_ptr<DiscreteTimeFilter> > > controller_ptrs_;
    int active_controller_idx_;  // Index of the active controller
    int nb_filters_; // Nb. of filters an individual controller have. This must be the same for all controllers in the pool.
    float scale_output_;  // Scale output with this value before returning.
    std::string controller_type_str_; // Type of the controller: Q_synthesis or PI
    int controller_type_;  // Q_synthesis: 0, PI: 1

    // internal states
    double beta0_=0, beta1_=0;

public:
    /*! Initialize from a ros parameter path, load controllers from the server.
     *
     * User should store each individual force controller in a yaml file, and load it
     * to ROS parameter server before this initialization. A controller `\controllers\awesome_controller`
     * should be defined like so:
     *
     * controllers:
     *   awesome_controller:
     *     filter_Q_00_a: [0, 1, 2]
     *     filter_Q_00_b: [0, 1, 2]
     *     filter_Q_01_a: [0, 1, 2]
     *     filter_Q_01_b: [0, 1, 2]
     *
     *     filter_zPyu_00_a: [0, 1, 2]
     *     filter_zPyu_00_b: [0, 1, 2]
     *     filter_zPyu_01_a: [0, 1, 2]
     *     filter_zPyu_01_b: [0, 1, 2]
     *
     * or if it is a PI type:
     *
     * controllers:
     *   awesome_controller:
     *     filter00_a: [0, 1, 2]
     *     filter00_b: [0, 1, 2]
     *     filter01_a: [0, 1, 2]
     *     filter01_b: [0, 1, 2]
     *
     * @param controller_path Path to the hybrid controller profile.
     */
    SISOForceControllerManager(std::string controller_path, ros::NodeHandle nh): active_controller_idx_(0), scale_output_(1) {
        // controller type
        if (!nh.getParam(controller_path + "/type", controller_type_str_))
        {
            ROS_ERROR_STREAM("Unable to load controller " << controller_path << ". Are the controllers loaded?");
        }
        else {
            ROS_INFO_STREAM("Loading a controller of [" << controller_type_str_ << "] type.");
        }
        if (controller_type_str_ == "Q_synthesis") controller_type_ = 0;
        else if (controller_type_str_ == "PI") controller_type_ = 1;
        else {
            ROS_ERROR_STREAM("Unknown controller type: " << controller_type_str_ << " .Shutting down.");
            ros::shutdown();
        }

        // load individual controllers
        std::vector<std::string> all_controller_paths;  // paths to controllers
        if (!nh.getParam(controller_path + "/controllers", all_controller_paths)){
            ROS_ERROR_STREAM("Unable to load individual controllers of " << controller_path << ". Are the controllers loaded?");
        }
        else {
            ROS_INFO_STREAM("Preparing to load " << all_controller_paths.size() << " individual controller!");
        }

        std::vector<double> b, a, taps;
        std::vector<std::shared_ptr<DiscreteTimeFilter> > filter_ptrs;  // a vector of pointers to filter, representing an individual force controller
        bool ret;
        for (int i = 0; i < all_controller_paths.size(); ++i) {
            // obtain scale factor for controller output
            ret = nh.getParam(all_controller_paths[i] + "/scale_output", scale_output_);
            if (!ret) ROS_ERROR_STREAM("Unable to load scale_output for individual controller: " << i);
            else{
                ROS_INFO_STREAM("Loading an individual controller of type: [" << controller_type_str_ <<"]. Output is scaled by: " << scale_output_);
            }

            // if is a PI type controller
            if (controller_type_str_ == "PI"){
                ret = true;
                ret = ret && nh.getParam(all_controller_paths[i] + "/filter00_a", a);
                ret = ret && nh.getParam(all_controller_paths[i] + "/filter00_b", b);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter00 = std::make_shared<DiscreteTimeFilter>(b, a);

                ret = ret && nh.getParam(all_controller_paths[i] + "/filter01_a", a);
                ret = ret && nh.getParam(all_controller_paths[i] + "/filter01_b", b);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter01 = std::make_shared<DiscreteTimeFilter>(b, a);

                filter_ptrs.clear();
                filter_ptrs.push_back(filter00);
                filter_ptrs.push_back(filter01);
                controller_ptrs_.push_back(filter_ptrs);

            }

            if (controller_type_str_ == "Q_synthesis"){
	        // sequence to store filters of a Q_synthesis controllers:
	        // Q filters first, row-major orders
	        // zPyu filters, row-major orders
                ret = true;

                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_Q_00_taps", taps);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter_Q_00 = std::make_shared<DiscreteTimeFilter>(taps);

                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_Q_01_taps", taps);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter_Q_01 = std::make_shared<DiscreteTimeFilter>(taps);

                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_zPyu_00_a", a);
                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_zPyu_00_b", b);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter_zPzu_00 = std::make_shared<DiscreteTimeFilter>(b, a);

                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_zPyu_10_a", a);
                ret = ret && nh.getParam(all_controller_paths[i] + "/filter_zPyu_10_b", b);
                if (!ret) ROS_ERROR_STREAM("Unable to load individual controller: " << i);
                auto filter_zPzu_10 = std::make_shared<DiscreteTimeFilter>(b, a);

                filter_ptrs.clear();
                filter_ptrs.push_back(filter_Q_00);
                filter_ptrs.push_back(filter_Q_01);
                filter_ptrs.push_back(filter_zPzu_00);
                filter_ptrs.push_back(filter_zPzu_10);
                controller_ptrs_.push_back(filter_ptrs);
            }
        }
    };

    void compute(const std::vector<double> & inputs, double & output){
        if (controller_type_ == 0)
            compute_Q(inputs, output);
        else if (controller_type_ == 1)
            compute_PI(inputs, output);
        else
        {
            ROS_ERROR_STREAM("Unknown controller type: " << controller_type_ << ". Shutting down!");
            ros::shutdown();
        }
    };

    /*! Compute the output of a 1-by-2 Q-based controller
     *
     * Implement the algorithm that I wrote in the paper. Equation (TODO: what is that?).
     *
     * @param inputs
     * @param output
     */
    void compute_Q(const std::vector<double> & inputs, double & output){
        output = (controller_ptrs_[active_controller_idx_][0]->compute(inputs[0] - beta0_)
                + controller_ptrs_[active_controller_idx_][1]->compute(inputs[1] - beta1_));
        beta0_ = controller_ptrs_[active_controller_idx_][2]->compute(output);
        beta1_ = controller_ptrs_[active_controller_idx_][3]->compute(output);
        output = output * scale_output_;
    }

    /*! Compute the output of a 1-by-2 linear controller
     *
     * Implement a basic 1-by-2 linear controller.
     *    output = P0 {input 0} + P1 {input 1}
     *
     * @param inputs
     * @param output
     */
    void compute_PI(const std::vector<double> & inputs, double & output){
        output = (controller_ptrs_[active_controller_idx_][0]->compute(inputs[0])
                  + controller_ptrs_[active_controller_idx_][1]->compute(inputs[1]));
        output = output * scale_output_;
    }
    /*! Switch active controller to one with given index
     *
     * @param controller_idx
     * @return true if success, false otherwise
     */
    bool switch_controller(int controller_idx){
        if (controller_idx < 0 or controller_idx > controller_ptrs_.size())
            // immediately return failure
            return false;
        else if (controller_idx == active_controller_idx_)
            // do nothing
            return true;
        // copy state
        for (int i = 0; i < nb_filters_; ++i) {
            if (not controller_ptrs_[controller_idx][i]->copy_state_from(controller_ptrs_[active_controller_idx_][i])) return false;
        }

        // switch controller
        active_controller_idx_ = controller_idx;
        return true;
    }
};


/*! Hybrid force controller with Q-synthesis/PI control law.
 *
 */
class HybridForceController {
    ros::NodeHandle nh_;
    std::string robot_ns_ = "denso";  // this namespace is used only when using the Denso's ros_control component
    bool debug_mode_;  // turn on class-wide debug

    std::string ft_topic_name_ = "/netft/raw";
    std::string viewer_name_ = "qtosg";
    bool viewer_on_;
    std::string scene_file_name_;  // openrave xml scenario file
    std::string robot_name_; // robot name in load scenario
    std::string manip_frame_name_;
    std::string ft_frame_name_;

    std::string robot_control_method_ = "direct";
    std::string robot_ip_addr_ = "192.168.0.21";

    std::shared_ptr<HWHandle::AbstractRobotController> robot_hw_ptr_;
    std::shared_ptr<FTSensorHandle> ft_hw_ptr_;
    std::vector<double> joint_init_;

    // debug-related members
    std::string debug_ns_ = "/debugger";
    std::shared_ptr<InfInteraction::TopicDebugger> debugger_ptr_;
    int debug_w_id_, debug_y_id_, debug_u_id_, debug_q_id_, debug_q_cmd_id_;

    // identificator
    std::string identificator_ns_ = "identificator";
    std::shared_ptr<InfInteraction::TopicDebugger> identificator_publisher_ptr_;
    std::vector<double> identify_data_;
    int identificator_pub_id_;

    // openrave
    OpenRAVE::EnvironmentBasePtr env_ptr_;
    OpenRAVE::RobotBasePtr robot_ptr_;
    OpenRAVE::RobotBase::ManipulatorPtr manip_ptr_;

    // core signal blocks
    std::shared_ptr<SignalBlock> wrench2force_map_ptr_;

    std::string force_controller_id_;
    std::shared_ptr<SISOForceControllerManager> Zaxis_force_controller_ptr_;
    double search_velocity_mm_sec_; /* Velocity to search for horizontal surface */
    double search_force_threshold_; /* Threshold for surface detection */
    double safety_force_threshold_;
    double initial_desired_force_;
    double distance_stopband_;  // (in meter) control = control if abs(control) > deadzone, otherwise 0
    double force_deadzone_; // if f_measure in [f_desired - f_deadzone, f_desired + f_deadzone], f_error = 0
    double joint_command_deadzone_;

    std::shared_ptr<SignalBlock> position2joint_map_ptr_;

    // variables used in RT loop
    std::vector<double> joint_cmd_;
    std::vector<double> joint_measure_;
    std::vector<double> wrench_measure_;  // wrench measurement
    std::vector<double> force_measure_;  // Cartesian force measurement
    std::vector<double> Zaxis_force_controller_inputs_;  // input vector to force controller
    double Zaxis_force_controller_output_ = 0;
    std::vector<double> cartesian_cmd_;

    int control_state_id_; // state index

    std::vector<double> setpoints_;  /* variables used for setpoints*/
    // NOTE on setpoint: In the second state of force control, this variable is freqently looked up
    // in the main control loop. In fact, exactly once per loop. The details of which index is used
    // for what signal is leave below.
    ros::Subscriber setpoints_subscriber_;

public:
    /*! Initialize a HybridForceController instance.
     *
     * This class provides functionality of a Cartesian Hybrid Force controller.
     *
     * User sends external commands to a HybridForceController via the ROS topic setup for
     * setpoints. See setpoint_subscriber_setup for more details on the topic exposed.
     *
     * User must load all ROS parameters before initializing this class.
     *
     * @param nh ROS Node Handle
     */
    HybridForceController(ros::NodeHandle nh): nh_(nh), setpoints_(10, 0), cartesian_cmd_(3, 0), joint_cmd_(6, 0), joint_measure_(6, 0),
    force_measure_(3, 0), wrench_measure_(6, 0), Zaxis_force_controller_inputs_(2, 0),
    safety_force_threshold_(20), identify_data_(2)
    {
        ROS_INFO("-- Initializing Hybrid Force Controller");

        // load parameters
        try_load_param("scene_path", scene_file_name_);
        try_load_param("robot_name", robot_name_);
        try_load_param("manip_frame", manip_frame_name_);
        try_load_param("ftsensor_frame", ft_frame_name_);
        try_load_param("viewer", viewer_on_, false);
        try_load_param("debug", debug_mode_, false);
        try_load_param("robot_control_method", robot_control_method_, std::string("ros_control"));

        // force control-related stuffs
        try_load_param("active_controller", force_controller_id_);
        try_load_param("search_velocity_mm_sec", search_velocity_mm_sec_, (double) 1);
        try_load_param("search_force_threshold", search_force_threshold_, (double) 2.0);
        try_load_param("safety_force_threshold", safety_force_threshold_, (double) 20);
        try_load_param("initial_desired_force", initial_desired_force_, (double) 2);
        try_load_param("positional_deadzone", distance_stopband_, (double) 0);
        try_load_param("force_deadzone", force_deadzone_, (double) 0);
        try_load_param("joint_command_deadzone", joint_command_deadzone_, (double) 0.0);

        // load force controller first, so that assertion error is thrown without affecting other
        Zaxis_force_controller_ptr_ = std::make_shared<SISOForceControllerManager>(force_controller_id_, nh_);

        // load hardware handles
        if (robot_control_method_ == "ros_control"){
            robot_hw_ptr_ = std::make_shared<HWHandle::JointPositionController> (robot_ns_, nh_);
        }
        else if (robot_control_method_ == "direct"){
            robot_hw_ptr_ = std::make_shared<HWHandle::RC8HWController> (robot_ip_addr_);
        }
        joint_init_ = robot_hw_ptr_->get_latest_jnt();

        ft_hw_ptr_ = std::make_shared<FTSensorHandle>(nh_, ft_topic_name_);
        if (debug_mode_)ft_hw_ptr_->set_debug(nh_);

        // communication
        debugger_setup();
        setpoint_subscriber_setup();
        identificator_publisher_ptr_ = std::make_shared<InfInteraction::TopicDebugger>(identificator_ns_, nh_);
        identificator_pub_id_ = identificator_publisher_ptr_->register_multiarray("raw_data");

        // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
        OpenRAVE::RaveInitialize(true); // start openrave core
        env_ptr_ = OpenRAVE::RaveCreateEnvironment(); // create the main environment
        if (viewer_on_) boost::thread thviewer(boost::bind(SetViewer, env_ptr_, viewer_name_));  // create viewer
        env_ptr_->Load(scene_file_name_); // load the scene
        robot_ptr_ = env_ptr_->GetRobot(robot_name_);
        manip_ptr_ = robot_ptr_->GetActiveManipulator();
        // set active dofs to [0,..,5]
        robot_ptr_->SetActiveDOFs(std::vector<int> {0, 1, 2, 3, 4, 5});

        // wrench2force map
        ros::Duration(0.3).sleep(); ros::spinOnce();  // make sure to receive at least a wrench reading before continue
        if (not ft_hw_ptr_->received_signal()){
            ROS_FATAL("Have not received any wrench signal. Shutting down.");
            ros::shutdown();
        }
        DoubleVector wrench_offset;
        ft_hw_ptr_->get_latest_wrench(wrench_offset);
        ft_hw_ptr_->set_wrench_offset(wrench_offset);

        // set the robot's initial configuration before initializing force projection block
        robot_ptr_->SetActiveDOFValues(joint_init_);
        wrench2force_map_ptr_ = std::make_shared<InfInteraction::Wrench2CartForceProjector>(robot_ptr_, ft_frame_name_);

        // inverse kinematics map
        // Inverse Kinematics: Converts Cartesian signals to joint coordinates
        double gam, gam2;  // weights for regulation of dq, and of (q_current - q_init)

        try_load_param("gam", gam, 0.01);
        try_load_param("gam2", gam2, 0.01);
        position2joint_map_ptr_ = std::make_shared<InfInteraction::CartPositionTracker>(robot_ptr_, manip_frame_name_, joint_init_, gam, gam2);
    };

    void setpoint_subscriber_setup(){
        setpoints_subscriber_ = nh_.subscribe("setpoints", 3, &HybridForceController::setpoint_callback_, this);
    }

    void debugger_setup(){
        // debugger publishing data to topics
        debugger_ptr_ = std::make_shared<InfInteraction::TopicDebugger>(debug_ns_, nh_);
        debug_w_id_ = debugger_ptr_->register_multiarray("wrench_measurement"); // meeasured wrench
        debug_y_id_ = debugger_ptr_->register_multiarray("cartesian_force_measurement"); // force_ output
        debug_u_id_ = debugger_ptr_->register_multiarray("cartesian_command"); // position command
        debug_q_id_ = debugger_ptr_->register_multiarray("joint_measure"); // joint command
        debug_q_cmd_id_ = debugger_ptr_->register_multiarray("joint_command"); // joint command
    }

    void debugger_publish_data(){
        debugger_ptr_->publish_multiarray(debug_w_id_, wrench_measure_);
        debugger_ptr_->publish_multiarray(debug_y_id_, force_measure_);
        debugger_ptr_->publish_multiarray(debug_u_id_, cartesian_cmd_);
        debugger_ptr_->publish_multiarray(debug_q_id_, joint_measure_);
        debugger_ptr_->publish_multiarray(debug_q_cmd_id_, joint_cmd_);
    }

    bool start(){
        char start_decision;
        std::cout << "Start hybrid force control? y/[n]: ";
        std::cin >> start_decision;
        if(start_decision != 'y'){
            ROS_INFO_STREAM("Trial arborts by user");
            return false;
        }

        timespec slp_dline_spec, sent_spec, wake_spec; // timespecs for time-keeping and time reporting
        // fifo policy for better (RT) responsiveness
        if (RTUtils::set_policy_fifo() == 0) {
            ROS_INFO("Scheduling policy set to fifo");
        }
        else{
            ROS_WARN("Set policy fail. Scheduling policy remained default");
        }
        const int FOUR_MS = 4000000;
        int diff_nsec; // time difference between (i) the deadline to wake up and (ii) the time when the message is received.
        int cy_idx = 0;  // cycle index
        double surface_height = 0;  // estimated surface height at the time of touch down. this value is used during the main control loop.
        double current_Zaxis_command;  // use to store current command position
        std::vector<double> new_joint_command(6);
        setpoints_[2] = initial_desired_force_; // initial force threshold set to the same as search_force_threshold_

        control_state_id_ = 1; // init control state to 1, touching down.
	control_state_id_ = 2;
	ROS_ERROR_STREAM("Initial state is set to 2 for testing. Change it back.");

	// two states loop designs
        //  [state 1] -> [states 2]
        //
        // Descriptpion:
        // - In [state 1], the robot moves downward until the vertical force component is greater
        //   than a prescribed threshold. When this occurs, the system enters [state 2];
        // - In [state 2], the robot moves according to instruction given by the setpoints. Particularly
        //      index 0: translation along the X-axis,
        //      index 1: translation along the Y-axis,
        //      index 2: desired vertical force,
        //      index 3: controller number. Must be integer.
        // - Whenever errors occurs, change to [state 3], then issues shutdown.
        clock_gettime(CLOCK_MONOTONIC, &slp_dline_spec);
        while(!ros::isShuttingDown()){
            ros::spinOnce();  // allow all callback in queue to run

            ft_hw_ptr_->get_latest_wrench(wrench_measure_);
            robot_hw_ptr_->get_latest_jnt(joint_measure_);

            // transform wrench measurement to Cartesian force
            wrench2force_map_ptr_->set_state(joint_measure_);
            force_measure_ = wrench2force_map_ptr_->compute(wrench_measure_);
            // safety measure
            if (ABS(force_measure_[0]) > safety_force_threshold_ or ABS(force_measure_[1]) > safety_force_threshold_ or ABS(force_measure_[2]) > safety_force_threshold_){
                ROS_ERROR_STREAM("Measure force too high. Shutdown to prevent damaging the robot!");
                control_state_id_ = 3;
            }

            if (control_state_id_ == 1){
                ROS_INFO_STREAM_THROTTLE(1, "Hybrid force controller moving down!");
                if (force_measure_[2] > search_force_threshold_) {
                    ROS_INFO_STREAM("Vertical component of measured force: "<< force_measure_[2] << ". Touch down successfuly.");
                    control_state_id_ = 2;
                    surface_height = cartesian_cmd_[2];
                }
                else {
                    cartesian_cmd_[2] += - search_velocity_mm_sec_ * MM * Ts;  // rate = 1e-2 mm/sec
                }
            }

            else if (control_state_id_ == 2){
                ROS_INFO_STREAM_ONCE("Hybrid force controller switches to active force control, listening to setpoints!");

                // switch controller
                if (!Zaxis_force_controller_ptr_->switch_controller(static_cast<int>(std::round(setpoints_[3])))){
                    ROS_ERROR_STREAM("Unable to switch controller. Shutting down.");
                    ros::shutdown();
                }

//                // compute Z-axsi force signal
//                Zaxis_force_controller_inputs_[1] = force_measure_[2];  // Z-axis force measurement
//                // if the force_measure_Z is not in [f_desired - deadzone, f_desired + deadzone], nothing changes
//                if (force_measure_[2] > force_deadzone_ + setpoints_[2]){
//                    Zaxis_force_controller_inputs_[0] = setpoints_[2] + force_deadzone_;  // Z-axis force measurement
//                }
//                else if (force_measure_[2] < setpoints_[2] - force_deadzone_){
//                    Zaxis_force_controller_inputs_[0] = setpoints_[2] - force_deadzone_;  // Z-axis force measurement
//                }
//                // otherwise, set the measured force to be the same as the desired value
//                else{
//                    Zaxis_force_controller_inputs_[0] = force_measure_[2];  // Z-axis force measurement
//                }

                // compute Z-axsi force signal
                Zaxis_force_controller_inputs_[0] = setpoints_[2];  // desired force
                // if the force_measure_Z is not in [f_desired - deadzone, f_desired + deadzone], nothing changes
                if (ABS(force_measure_[2] - setpoints_[2]) > force_deadzone_){
                    Zaxis_force_controller_inputs_[1] = force_measure_[2];  // Z-axis force measurement
                }
                    // otherwise, set the measured force to be the same as the desired value
                else{
                    Zaxis_force_controller_inputs_[1] = setpoints_[2];  // Z-axis force measurement
                }
                current_Zaxis_command = Zaxis_force_controller_output_;


                Zaxis_force_controller_ptr_->compute(Zaxis_force_controller_inputs_, Zaxis_force_controller_output_);

                if (ABS(Zaxis_force_controller_output_ - current_Zaxis_command) < distance_stopband_)
                {
                    ROS_DEBUG_STREAM_THROTTLE(1, "Commanding distance changes very little, setting to unchanged!");
                    Zaxis_force_controller_output_ = current_Zaxis_command;
                }
                else{
                    ROS_DEBUG_STREAM_THROTTLE(1, "Commanding distance changes: " << Zaxis_force_controller_output_ - current_Zaxis_command);
                }

                // compute cartesian command
                cartesian_cmd_[0] = setpoints_[0];
                cartesian_cmd_[1] = setpoints_[1];
                cartesian_cmd_[2] = surface_height + Zaxis_force_controller_output_;
            }
            else {
                ROS_FATAL_STREAM("Encounter illegal/non-action control_state_id_: " << control_state_id_ << ". Shutting down.");
                ros::shutdown();
            }
            //publish data for identification
            // a 2-vector which contains the measured vertical force and the measured cartesian command
            identify_data_[0] = force_measure_[2];
            identify_data_[1] = cartesian_cmd_[2];
            identificator_publisher_ptr_->publish_multiarray(identificator_pub_id_, identify_data_);

            // compute command to send
            position2joint_map_ptr_->set_state(joint_measure_);
            new_joint_command = position2joint_map_ptr_->compute(cartesian_cmd_);
            for(int i=0; i < 6; i++){
                if (ABS(joint_cmd_[i] - new_joint_command[i]) > joint_command_deadzone_){
                    joint_cmd_[i] = new_joint_command[i];
                }
            }

            // stage 1 finished; sleep
            RTUtils::increment_timespec(slp_dline_spec, FOUR_MS);
            clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &slp_dline_spec, NULL); // sleep until slp_dline_spec

            // send joint position command
            clock_gettime(CLOCK_MONOTONIC, &wake_spec);
            robot_hw_ptr_-> send_jnt_command(joint_cmd_);
            robot_ptr_->SetActiveDOFValues(joint_cmd_); // update in openrave viewer
            clock_gettime(CLOCK_MONOTONIC, &sent_spec);

            // publish debug information
            if (debug_mode_) {
                debugger_publish_data();

                // publish time discrepancy
                if (cy_idx % 125 == 0){
                    RTUtils::diff_timespec(diff_nsec, slp_dline_spec, wake_spec);
                    ROS_DEBUG_STREAM("wake - deadline-to-send: " << diff_nsec << " nsec (this value should be very small, ideally less than 1e4)");
                    RTUtils::diff_timespec(diff_nsec, slp_dline_spec, sent_spec);
                    ROS_DEBUG_STREAM("sent - deadline-to-send: " << diff_nsec << " nsec (this value should be small, ideally less than 1e5)");
                    ROS_DEBUG_STREAM("current state: " << control_state_id_);
                }
            }
            cy_idx = (cy_idx + 1) % 2147483640;

            // stage 2 finsished, sleep
            RTUtils::increment_timespec(slp_dline_spec, FOUR_MS);
            clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &slp_dline_spec, NULL);
        }

        // state 2: track force
        return true;
    }

    /*! Load parameter, shutdown if fail.
     *
     */
    template <typename T> void try_load_param(std::string param_path, T & param){
        if (!nh_.getParam(param_path, param)){
            ROS_ERROR_STREAM("Unable to find ["<< param_path <<"] on param server. Shutting down!");
            ros::shutdown();
        }
    }

    /*! Load parameter, if fail set default value.
     *
     * @tparam T
     * @param param_path
     * @param param
     * @param param_default
     */
    template <typename T> void try_load_param(std::string param_path, T & param, T param_default){
        if (!nh_.param<T>(param_path, param, param_default)){
            ROS_WARN_STREAM("Unable to find ["<< param_path <<"] on param server. Setting default value: " << param_default);
        }
        else {
            ROS_INFO_STREAM("Param loaded from [" << param_path << "]. Value: " << param);
        }
    }

private:
    void setpoint_callback_(const std_msgs::Float64MultiArray &multiarray_msg){
        ROS_DEBUG_STREAM_THROTTLE(1, "in setpoint callback");
        if (control_state_id_ == 1){
            ROS_ERROR_STREAM("Attempt to change setpoint during state 1. Error occurs.");
        }
        setpoints_[0] = multiarray_msg.data[0];
        setpoints_[1] = multiarray_msg.data[1];
        setpoints_[2] = multiarray_msg.data[2];
    }
};


int main(int argc, char **argv){
    ros::init(argc, argv, "hybrid_force_controller");
    ros::NodeHandle nh("~");
    HybridForceController controller(nh);
    if (!ros::isShuttingDown()){
        controller.start();
    }
    return 0;
}

