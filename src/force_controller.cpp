#include <ros/ros.h>
#include <infinite_interaction/infinite_interaction_lib.h>

// openrave
#include <openrave-core.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

const double MM = 1e-3;
const double Ts = 0.008;


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


class ControllerManager {
    // A vector of controllers, each is a vector of pointers to filters. If there are two filters per vector,
    // it is a PI. If there are four filters per vector, it is a Q-synthesis controller. Sequence: Q0, Q1, zPyu0, zPyu1
    std::vector<std::vector<std::shared_ptr<DiscreteTimeFilter> > > controller_ptrs_pool_;
    int active_controller_idx_ = 0;
    int nb_filters_;
    int controller_type_;  // 0: Q-synthesis, 1: PI

    // internal states
    double state0_=0, state1_=0;
public:
    /*! Initialize from a ROS parameter path
     *
     * @param controller_path
     */
    ControllerManager(std::string controller_path){};

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

    void load_controller(const std::vector<std::shared_ptr<DiscreteTimeFilter > > & new_controller){
        // a controller is a vector of filter.
        controller_ptrs_pool_.push_back(new_controller);
    }

    int get_nb_controller() {return controller_ptrs_pool_.size();}

    void compute_Q(const std::vector<double> & inputs, double & output){
        output = (controller_ptrs_pool_[active_controller_idx_][0]->compute(inputs[0] - state0_)
                + controller_ptrs_pool_[active_controller_idx_][1]->compute(inputs[1] - state1_));
        state0_ = controller_ptrs_pool_[active_controller_idx_][2]->compute(output);
        state1_ = controller_ptrs_pool_[active_controller_idx_][3]->compute(output);
    }
    void compute_PI(const std::vector<double> & inputs, double & output){
        output = (controller_ptrs_pool_[active_controller_idx_][0]->compute(inputs[0])
                  + controller_ptrs_pool_[active_controller_idx_][1]->compute(inputs[1]));
    }
    /*! Switch active controller to one with given index
     *
     * @param controller_idx
     * @return true if success, false otherwise
     */
    bool switch_controller(int controller_idx){
        if (controller_idx < 0 or controller_idx > controller_ptrs_pool_.size())
            // immediately return failure
            return false;
        else if (controller_idx == active_controller_idx_)
            // do nothing
            return true;
        // copy state
        for (int i = 0; i < nb_filters_; ++i) {
            if (not controller_ptrs_pool_[controller_idx][i]->copy_state_from(controller_ptrs_pool_[active_controller_idx_][i])) return false;
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
    // name space
    std::string robot_ns_ = "denso";
    bool debug_mode_;

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

    // openrave
    OpenRAVE::EnvironmentBasePtr env_ptr_;
    OpenRAVE::RobotBasePtr robot_ptr_;
    OpenRAVE::RobotBase::ManipulatorPtr manip_ptr_;

    // core signal blocks
    std::shared_ptr<SignalBlock> wrench2force_map_ptr_;

    std::string force_controller_id_;
    std::shared_ptr<ControllerManager> force_controller_ptr_;
    double search_velocity_mm_sec_; /* Velocity to search for horizontal surface */
    double search_force_threshold_; /* Threshold for surface detection */

    std::shared_ptr<SignalBlock> position2joint_map_ptr_;

    // variables used in RT loop
    std::vector<double> joint_cmd_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> joint_measure_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench_measure_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> force_measure_ = {0, 0, 0};
    double force_vert_measure_, force_vert_desired_ = -1;
    std::vector<double> force_controller_inputs_ = {0, 0};
    double force_controller_output_ = 0;
    std::vector<double> cartesian_cmd_ = {0, 0, 0};

public:
    /*! Initialize Hybrid Force controller
     *
     * Relevant parameters must be loaded prior to initialization.
     *
     * @param nh
     */
    HybridForceController(ros::NodeHandle nh): nh_(nh){
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

        debugger_setup();

        // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
        OpenRAVE::RaveInitialize(true); // start openrave core
        env_ptr_ = OpenRAVE::RaveCreateEnvironment(); // create the main environment
        OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
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
            return;
        }
        dVector wrench_offset;
        ft_hw_ptr_->get_latest_wrench(wrench_offset);
        ft_hw_ptr_->set_wrench_offset(wrench_offset);

        // set the robot's initial configuration before initializing force projection block
        robot_ptr_->SetActiveDOFValues(joint_init_);
        wrench2force_map_ptr_ = std::make_shared<InfInteraction::Wrench2CartForceProjector>(robot_ptr_, ft_frame_name_);

        // force scheme
        force_controller_ptr_ = std::make_shared<ControllerManager>("/" + force_controller_id_);

        // inverse kinematics map
        // Inverse Kinematics: Converts Cartesian signals to joint coordinates
        double gam, gam2;  // weights for regulation of dq, and of (q_current - q_init)

        try_load_param("gam", gam, 0.01);
        try_load_param("gam2", gam2, 0.01);
        position2joint_map_ptr_ = std::make_shared<InfInteraction::CartPositionTracker>(robot_ptr_, manip_frame_name_, joint_init_, gam, gam2);
    };

    void debugger_setup(){
        // debugger publishing data to topics
        debugger_ptr_ = std::make_shared<InfInteraction::TopicDebugger>(debug_ns_, nh_);
        debug_w_id_ = debugger_ptr_->register_multiarray("w_n"); // meeasured wrench
        debug_y_id_ = debugger_ptr_->register_multiarray("y_n"); // force_ output
        debug_u_id_ = debugger_ptr_->register_multiarray("u_n"); // position command
        debug_q_id_ = debugger_ptr_->register_multiarray("q_n"); // joint command
        debug_q_cmd_id_ = debugger_ptr_->register_multiarray("qcmd_n"); // joint command
    }

    void debugger_log(){
        debugger_ptr_->publish_multiarray(debug_w_id_, wrench_measure_);
        debugger_ptr_->publish_multiarray(debug_y_id_, force_measure_);
        debugger_ptr_->publish_multiarray(debug_u_id_, cartesian_cmd_);
        debugger_ptr_->publish_multiarray(debug_q_id_, joint_measure_);
        debugger_ptr_->publish_multiarray(debug_q_cmd_id_, joint_cmd_);
    }

    bool main(){
        char start_dec;
        std::cout << "Start hybrid force control? y/[n]: ";
        std::cin >> start_dec;
        if(start_dec != 'y'){
            ROS_INFO_STREAM("Trial arborts by user");
            return false;
        }
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

        // state 1: touch down
        // The robot moves slowly downward until it touches the board.
        while(!ros::isShuttingDown()){
            ros::spinOnce();

            ft_hw_ptr_->get_latest_wrench(wrench_measure_);
            robot_hw_ptr_->get_latest_jnt(joint_measure_);

            // wrench transform
            wrench2force_map_ptr_->set_state(joint_measure_);
            force_measure_ = wrench2force_map_ptr_->compute(wrench_measure_);

            if (force_measure_[2] > search_force_threshold_) {
                ROS_INFO_STREAM("Vertical component of measured force: "<< force_measure_[2] << ". Touch down successfuly.");
                break;
            }

            cartesian_cmd_[0] = 0;
            cartesian_cmd_[1] = 0;
            cartesian_cmd_[2] += - search_velocity_mm_sec_ * MM * Ts;  // rate = 1e-2 mm/sec

            // compute command to send
            position2joint_map_ptr_->set_state(joint_measure_);
            joint_cmd_ = position2joint_map_ptr_->compute(cartesian_cmd_);

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
                debugger_log();

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

        // state 2: move back
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
};


int main(int argc, char **argv){
    ros::init(argc, argv, "hybrid_force_controller");
    ros::NodeHandle nh("~");
    HybridForceController controller(nh);
    return controller.main();
}

