//
// Created by Hung Pham on 8/9/18.
//

#ifndef PROJECT_INFINITE_INTERACTION_LIB_H
#define PROJECT_INFINITE_INTERACTION_LIB_H

#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Header.h>
#include <sensor_msgs/JointState.h>
#include <vector>
#include <algorithm>
#include <openrave-core.h>

// time header for RT utility functions
#include <time.h>

// RC8
#include "denso_control/rc8_controller_interface.h"

typedef std::vector<double> dVector;
typedef std_msgs::Float64MultiArray MultiArrayMsg;


const int NANOSECS_IN_SEC = 1000000000;

inline double MIN(double x, double y)
{
    return x > y ? y : x;
}


inline double MAX(double x, double y)
{
    return x > y ? x : y;
}

// A common terminology:
// x: input signal
// y: output signal




// A generic Base class for a discrete-time processing "block" that receives at each time
// step a vector input and produces a vector output. In general, most controllers, kinematic
// maps, force_ maps implemented in this package are child classes of this class.
class SignalBlock {
public:
    // NOTE: virtual method allows derived class to call its implementation,
    // instead of the base class implementation.
    // compute an output, which is a vector, given an input which is also a vector
    virtual dVector compute(const dVector & y_n) = 0;
    // set the internal state of the block to x_n. This should be done before calling compute.
    // Mathematically, this object compute u_n as a function of y_n and x_n.
    virtual void set_state(const dVector & x_n) = 0;
    // virtual destructor allows Derived Classes to destroy its internal
    // states. This is always recommend when polymorphism is needed.
    // https://stackoverflow.com/questions/461203/when-to-use-virtual-destructors
    virtual int get_input_size() {return 1;};
    virtual int get_output_size() {return 1;};
    virtual int get_state_size() {return 1;};
    virtual ~SignalBlock() {};
};

namespace InfInteraction{

// Project measured wrench to each joint torque.
class JointTorqueFromWrenchProjector: public SignalBlock {
    OpenRAVE::RobotBasePtr robot_ptr_;
    OpenRAVE::RobotBase::ManipulatorPtr ft_sensor_ptr_;
    dVector jacobian, jacobian_rot, jacobian_T, jacobian_rot_T,  // Translational and rotational Jacobians, transposed
            force, torque,  // input force_, torque
            tau1, tau2, tau, // projected torque
            joint_current_;  // current joint position
    OpenRAVE::Transform T_wee_;
    OpenRAVE::RaveVector<double> rave_force_, rave_torque_, temp_vec_;
public:
    JointTorqueFromWrenchProjector(OpenRAVE::RobotBasePtr robot_ptr, std::string ft_sensor_frame);
    dVector compute(const dVector & wrench_measure);
    void set_state(const dVector & joint_measure);
    virtual int get_input_size() {return 6;};
    virtual int get_output_size() {return 6;};
    virtual int get_state_size() {return 6;};
};


// TODO: Current the projector assume that the robot's end-effector pose remains at the initial pose.
// This assumtption is not true for a more general class of problem. (But true for the Cartersian admittance
// task I am working on.

/*! Project measured wrench to Cartesian workspace
 *
 * NOTE: The robot needs to be initialized to the initial configuration, before initializing
 * this function. As
 */
class Wrench2CartForceProjector: public SignalBlock {
    OpenRAVE::Transform T_wee_;
    OpenRAVE::geometry::RaveVector<double> force_, force_rotated_;
public:
    Wrench2CartForceProjector(OpenRAVE::RobotBasePtr robot_ptr, std::string ft_frame_name);
    dVector compute(const dVector & wrench_measure);
    void set_state(const dVector & dummy_var);
};


/*! Form diagonal matrix of given controllers.
 *
 * Suppose (c1, c2, .., cN) is given as input. This controller expect an input with length N
 * and produces an output with length N at each call to compute. One has respectively:
 *
 *   y[i] = ci.compute{x[i]}
 */
class ControllerCollection: public SignalBlock {
    std::vector<std::shared_ptr<SignalBlock > > controllers;
    unsigned int nb_controllers;
public:
    ControllerCollection(std::vector<std::shared_ptr<SignalBlock > > controllers_);
    dVector compute(const dVector & x_n);
    void set_state(const dVector & x_n);
};

/*! A feedback interconnection with delayed feedback route.
 *
 * --- o -- forward block -------->
 *     ^                     | z^-1
 *     \--- feedback block --/
 *
 */
class DelayFeedback: public SignalBlock {
    std::shared_ptr<SignalBlock> _fwd_block, _fb_block;
    dVector y_last;
    unsigned int _input_size;
    int _sign;
public:
    DelayFeedback(std::shared_ptr<SignalBlock> fwd_block, std::shared_ptr<SignalBlock> fb_block, int sign=-1);
    dVector compute(const dVector & x_n);
    void set_state(const dVector & x_n) {};
};

// A Block that simply offset the input by an initially given vector.
class SimpleOffset: public SignalBlock {
    dVector offset;
public:
    SimpleOffset(const dVector & offset_);
    // Always return a vector that has the same size as the input vector. Each element in the new vector is
    // offset by the corresponding element in the internal offset vector.
    dVector compute(const dVector & y_n);
    void set_state(const dVector & x_n);
};


/*! Cartesian Position Tracker
 *
 * Compute Joint Position Command to track a given position pos[n] while keeping the orientation fixed at the
 * initial value.
 *
 * NOTE: At each iteration, users should call set_state
 * and update jnt_pos_current before calling compute to obtain another joint position.
 */
 class CartPositionTracker: public SignalBlock {
     dVector _jnt_pos_current,  /* Current robot joint position */
             _jnt_pos_init,     /* Initial robot joint position */
             jnt_pos_save_; /* Temporary vector to store data */
     OpenRAVE::RaveVector<OpenRAVE::dReal>
             _quat_init, /*Initial orientation*/
             _pos_init /*Initial position*/ ;
     OpenRAVE::RobotBasePtr robot_ptr;
     OpenRAVE::RobotBase::ManipulatorPtr manip_ptr;
     double _gam2,  /*Weight to limit (q_cur + dq - q_init)*/
             _gam;  /*Weight to limit dq*/
 public:
     /*!
      *
      * @param robot_ptr_
      * @param manip_frame
      * @param jnt_pos_init
      * @param gam
      * @param gam2
      */
     CartPositionTracker(OpenRAVE::RobotBasePtr robot_ptr_, std::string manip_frame, dVector jnt_pos_init, double gam=0.0, double gam2=0.0);

     /*! \brief Compute joint values that track the given the Cartesian position pos_n.
      *
      * This function solves a QP with Jacobians computed using OpenRAVE. A single
      * optimization problem is employed.
      *
      * \param pos_n Cartesian position of the end-effector w.r.t to the initial Cartesian position.
      * */
     dVector compute(const dVector & pos_n);
     void set_state(const dVector & jnt_pos_n);
 };


 /* A convenient class for publishing messages for debugging or logging with simple usage.
  * Register the topic by name, then publish messages from STL Vector<double>.
  */
 class TopicDebugger {
     std::vector<std::shared_ptr<ros::Publisher > > pub_vecs; /* Publishers vector */
     std::string debug_ns;  /*Base namespace*/
     ros::NodeHandle nh;

 public:
     TopicDebugger(std::string debug_ns, ros::NodeHandle & node_handle);
     /*! Register new topic
      *
      * @param topic_name
      * @return topic index.
      */
     int register_multiarray(std::string topic_name);
     /*! Publish message on selected topic
      *
      * @param topic_id
      * @param data
      */
     void publish_multiarray(int topic_id, const dVector & data);
 };
}


/*! Discrete-time MIMO FIR series+feedback (srfb) controller.
 *
 * The FIRsrfb filter receives a stream of vector-valued input (y[0],
 * y[1], ... ,y[T-1]) and produces a stream of vector-valued output
 * (u[0], u[1], ... ,u[T-1])
 *
 * NOTE: This controller has the form of a series+feedback with two FIR blocks:
 *  y[n]                             u[n]
 *     ---[L] ----o--------------------->
 *               -|              |
 *                \-----[MB2]----/
 *
 * The input and output streams are stored respectively in two
 * circular input and output banks `ybank` and `ubank`. Both are
 * attributes of the class. New data, which are u[n] and y[n], are
 * appended to the banks reversely. Consider `ybank`:
 *
 * y[n][0], y[n][1], y[n][2], ..., y[1][0], y[1][1], y[1][2], y[0][0], y[0][1], y[0][2]
 *
 * The input bank `ybank` has length ny * T; the output bank has
 * length nu * T, which are sufficient for T periods of inputs and
 * outputs.
 *
 * Filter coefficients vector L has shape (L, ny, nu) and filter
 * coefficients vector MB2 has shape (L, nu, nu). Both are stored
 * flatten in row-order.
 *
 * At the beginning of each call to `compute`:
 *   1. an input y[n] is fed to the controller;
 *   2. yidx point to y[n-1][0];
 *   3. uidx point to u[n-1][0].
 * Before the call to `compute` terminates:
 *   1. an output u[n] is to be returned;
 *   2. yidx point to y[n][0];
 *   3. uidx point to u[n][0].
 *
 * Input u[n] is computed using the following computations:
 *
 *   alpha = L[0] * y[n] + ... + L[T - 1] * y[n - T + 1]
 *   beta = M[1] * u[n - 1] + ... + M[T - 1] * u[n - T + 1]
 *   u[n] = alpha - beta
 * 
 * It can be shown that the above computations are equivalent to the
 * following equality:
 *
 *    L[0] * y[n] + ... + L[T - 1] * y[n - T + 1]
 *  = u[n] + M[1] * u[n - 1] + ... + M[T - 1] * u[n - T + 1]
 *
 * This form can easily be adapted to implement decoupled (diagonal)
 * control rules. As an example, consider the scalar control rule
 *
 *    u[n] / y[n] = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2)
 * 
 * This rule can be implemented with the coefficients: L:=[b0, b1, b2], MB2:=:[0, a1, a2]
 * 
 * Suppose there are two dofs, and the input and output are both
 * 2-vectors, one has instead the following coefficients:
 * 
 * L:=[b0 * I2, b1 * I2, b2 * I2] and MB2:=[0*I2, a1*I2, a2*I2]
 */
class FIRsrfb: public SignalBlock {
    unsigned int T, ny, nu, /*FIR banks shape*/
            yidx = 0, uidx = 0, ybank_sz, ubank_sz;  /*Current index of y[n][0] and u[n][0] respectively.*/
    dVector L, MB2,   /*FIR coefficients in row-order*/
            ybank, ubank,  /*Storage banks*/
            alpha, beta, un;  /*temp*/

    // initialize filter internal variables (input and output banks, etc)
    // This function should only be used in a constructor.
    void init_filter(dVector uinit);

public:
    // initialize the filter with zero initial output
    FIRsrfb(unsigned int T_, unsigned int ny_, unsigned int nu_, dVector L_, dVector MB2_);
    // initialize the filter with initial outputs being uinit for all time steps
    FIRsrfb(unsigned int T_, unsigned int ny_, unsigned int nu_, dVector L_, dVector MB2_, dVector uinit);
    dVector compute(const dVector & y_n);  /*See class docstring*/
    void set_state(const dVector & x_n);
};


/*! Discrete-time SISO filter.
 *
 * An implementation of a basic Discrete-time SISO filter.  Let x[n]
 * be the input and y[n] be the output, this filter implements the
 * following CCDE numerically:
 *
 *     b0 x[n] + b1 x[n-1] + ... b_M x[n-M] = a0 y[n] + a1 y[n-1] + ... + a_N y[n - N]
 * 
 * The z-transform of this filter has the standard form:
 *
 *     b0 + b1 z^-1 + ... + bM z^-M   Y
 *     ---------------------------- = -
 *     a0 + a1 z^-1 + ... + aN z^-N   X
 *
 * NOTE: This filter is initialized in a signal-processing oriented
 * way. For instance, suppose ([1], [2, 3]) are given, the following
 * ccde is realized:
 * 
 *          1 * x[n] = 2 * y[n] + 3 * y[n-1]
 *  
 * On another hand, in a control-theory oriented way, such as
 * python-control, the following ccde is realized:
 * 
 *          1 * x[n-1] = 2 * y[n] + 3 * y[n-1]
 * 
 * NOTE: a FIR filter is basically a ccde with zeroth order
 * denominator. This class can also be used to initialized such
 * configuration.
 */
class DiscreteTimeFilter: public SignalBlock {
    /*! Order of the filter. Is basically the order of the denominator. */
    int N, M;
    int n;       /* Current memory index. */
    int mem_sz;  /* Both memory banks, input and output banks, have `mem_sz` numbers. */
    /*! Store the last (order)-th inputs and outputs */
    dVector x_mem, y_mem;
    dVector b, a; // Coefficients of the filter. If is a FIR, a is unit.
public:
    /*! Initialize a standard ccde.
     *
     * \param b_in coefficients of the numerator
     * \param a_in coefficients of the denominator
     * \param y_initial_val initial condition of the output. NOTE: initial condition 
                            of input is default to zero.
     */
    DiscreteTimeFilter(dVector b_in, dVector a_in, double y_initial_val=0);
    /*! Initialize a FIR filter.
     *
     * \param taps Taps of the filter.
     */
    DiscreteTimeFilter(dVector taps);
    /*! Compute the output in the next time step.
     *
     *
     * @param x_n
     * @return y_n
     */
    double compute(double x_n);
    dVector compute(const dVector & x_n);
    void set_state(const dVector & x_n);
    /*! Copy the internal state of the given filter.
     *
     * @param other_filter
     * @return true if sucess, false otherwise.
     */
    bool copy_state_from(const std::shared_ptr<DiscreteTimeFilter> & other_filter);
};



/*! FT signal stream handler.
 *
 * A handler applies basic operations including collecting data and
 * filter them using DiscreteTimeFilter objects. Note that there are 6
 * individual filters, but all have the same coefficients.
 */
class FTSensorHandle {
    double fx, fy, fz, tx, ty, tz;
    std::vector<double> wrench_offset;
    std::vector<DiscreteTimeFilter> lp_filters;
    dVector b, a; /*!filter coefficients*/
    bool debug = false;
    ros::Publisher wrench_pub_debug;  // publisher for debugging
    ros::Subscriber ft_subscriber;
public:

    /*! Basic initialization. All internal data members are initialized to the default values.
     *
     * No ft_sensor topic is attached.
     */
    FTSensorHandle();

    /*! Basic initialization.
     *
     * All internal data members are initialized to the default values. This class listens to signal
     * from topic `ft_topic` and update its internal data members (fx, fy, fz, tx, ty, tz) respectively.
     */
    FTSensorHandle(ros::NodeHandle & nh, std::string ft_topic);
    explicit FTSensorHandle(const dVector &wrench_offset_input);; // with non-zero offset values
    FTSensorHandle(const dVector &wrench_offset_input, dVector b_in, dVector a_in); // with offset and low-pass filter
    void signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg);
    void get_latest_wrench(dVector &force, dVector &torque);
    void get_latest_wrench(dVector &wrench);
    void set_debug(ros::NodeHandle &nh);
    void set_wrench_offset(dVector wrench_offset_);
    void log_latest_wrench(const std_msgs::Header &);
    bool received_signal();
};

namespace HWHandle {

    class AbstractRobotController {
    public:
        virtual void send_jnt_command(std::vector<double> &jnt_cmds) = 0;

        virtual void get_latest_jnt(std::vector<double> &jnt_positions) = 0;

        virtual std::vector<double> get_latest_jnt() = 0;

        virtual ~AbstractRobotController() {};
    };

    /*! Joint Position Robot controller.
     *
     *  This class sends joint position commands through the topics
     *  exposed by ros-control, and receives the current joint positions
     *  by listening to joint state message. Retrieving and sending are
     *  accessible via member functions.
     */
    class JointPositionController : public AbstractRobotController {
        std::string _name_space;
        std::vector<ros::Publisher> _jnt_pubs;
        ros::NodeHandle _nh;
        std::vector<double> _joint_position;
        ros::Subscriber _jnt_pos_subscriber;
    public:
        /*! Constructor for a joint position controller class.
         */
        explicit JointPositionController(std::string name_space, ros::NodeHandle &nh);

        ~JointPositionController() {};

        void send_jnt_command(std::vector<double> &jnt_cmds);

        void get_latest_jnt(std::vector<double> &jnt_positions);

        std::vector<double> get_latest_jnt();

        void signal_callback(const sensor_msgs::JointStateConstPtr &msg);

        bool received_msg();
    };

    class RC8HWController : public AbstractRobotController {
        std::vector<double> _jnt;
        std::shared_ptr<denso_control::RC8ControllerInterface> _rc8_controller_ptr;
        bool _ret;
//        unsigned int slave_mode = SlaveMode::J0; // most smooth, but must be ran on RT machine
        unsigned int slave_mode_ = SlaveMode::J1; // less smooth
    public:
        /* Connect to the RC8 controller, start motor.
         *
         */
        explicit RC8HWController(std::string ip_addr);

        ~RC8HWController();

        void send_jnt_command(std::vector<double> &jnt_cmds);

        void get_latest_jnt(std::vector<double> &jnt_positions);

        std::vector<double> get_latest_jnt();
    };
}


class JointPositionHandler{
    std::vector<double> joint_position;
public:
    JointPositionHandler();
    void signal_callback(const sensor_msgs::JointStateConstPtr &msg);
    dVector get_latest_jnt();
    bool received_msg();
};




namespace RTUtils{
    /*! Increment the given timespec.
     *
     * @param tspec a timespace to increment
     * @param inc_nsec increment internal in nanoseconds
     */
    void increment_timespec(timespec& tspec, const int& inc_nsec);

    /*! Compute the difference between two timespecs.
     *
     * Result = timespec2 - timespec1
     *
     * @param diff_nsec difference in nanoseconds.
     * @param tspec1
     * @param tspec2
     */
    void diff_timespec(int& diff_nsec, timespec& tspec1, timespec& tspec2);

    /*! Set the current process' scheduling policy to fifo
     *
     * @return 0 if success, 1 otherwise.
     */
    int set_policy_fifo();
}

[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
void matrix_mult(std::vector<double> &A, std::vector<double> &x, std::vector<double> &y);
[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
void matrix_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
std::vector<double> mat_transpose(const std::vector<double>& M, int ncol);

#endif //PROJECT_INFINITE_INTERACTION_LIB_H
