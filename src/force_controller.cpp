#include <ros/ros.h>
#include <infinite_interaction/infinite_interaction_lib.h>

// openrave
#include <openrave-core.h>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>


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


/*! Hybrid force controller with Q-synthesis/PI control law.
 *
 */
class ForceController {
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

    std::string active_force_controller_id_, active_force_controller_type_;

    // core signal blocks
    std::shared_ptr<SignalBlock> wrench2force_map_ptr_;
    std::shared_ptr<SignalBlock> active_force_controller_ptr_;
    std::shared_ptr<SignalBlock> position2joint_map_ptr_;

    // variables used in RT loop
    std::vector<double> joint_cmd_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> joint_measure_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> wrench_measure_ = {0, 0, 0, 0, 0, 0};
    std::vector<double> force_controller_inputs_ = {0, 0};
    double force_controller_output_ = 0;
    std::vector<double> cartesian_cmd_ = {0, 0, 0};

public:
    ForceController(ros::NodeHandle nh): nh_(nh){
        ROS_INFO("-- Initializing controller node");

        // load parameters
        try_load_param("/scene_path", scene_file_name_);
        try_load_param("/robot_name", robot_name_);
        try_load_param("/manip_frame", manip_frame_name_);
        try_load_param("/ftsensor_frame", ft_frame_name_);
        try_load_param("/viewer", viewer_on_, false);
        try_load_param("/debug", debug_mode_, false);

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

        // debugger publishing data to topics
        debugger_ptr_ = std::make_shared<InfInteraction::TopicDebugger>(debug_ns_, nh_);
        debug_w_id_ = debugger_ptr_->register_multiarray("w_n"); // meeasured wrench
        debug_y_id_ = debugger_ptr_->register_multiarray("y_n"); // force_ output
        debug_u_id_ = debugger_ptr_->register_multiarray("u_n"); // position command
        debug_q_id_ = debugger_ptr_->register_multiarray("q_n"); // joint command
        debug_q_cmd_id_ = debugger_ptr_->register_multiarray("qcmd_n"); // joint command

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
        dVector wrench_offset;
        ft_hw_ptr_->get_latest_wrench(wrench_offset);
        ft_hw_ptr_->set_wrench_offset(wrench_offset);
        robot_ptr_->SetActiveDOFValues(joint_init_); // set the robot's initial configuraiton before initializing force_ projection block
        wrench2force_map_ptr_ = std::make_shared<InfInteraction::Wrench2CartForceProjector>(robot_ptr_, ft_frame_name_);

    };

    bool main(){
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
        if (!nh_.getParam(param_path, param)){
            ROS_WARN_STREAM("Unable to find ["<< param_path <<"] on param server. Setting default value.");
            param = param_default;
        }
    }
};


int main(int argc, char **argv){
    ros::init(argc, argv, "admittance_controller");
    ros::NodeHandle nh("~");
    ForceController controller(nh);
    return controller.main();
}

