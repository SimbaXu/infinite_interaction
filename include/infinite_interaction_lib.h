//
// Created by Hung Pham on 8/9/18.
//

#ifndef PROJECT_INFINITE_INTERACTION_LIB_H
#define PROJECT_INFINITE_INTERACTION_LIB_H

#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Header.h>
#include <sensor_msgs/JointState.h>
#include <vector>
#include <algorithm>
#include <openrave-core.h>

typedef std::vector<double> dVector;



// TODO: Rename this class to something more appropriate
// A generic Base class for a discrete-time processing "block" that receives at each time
// step a vector input and produces a vector output. In general, most controllers, kinematic
// maps, force maps implemented in this package are child classes of this class.
class LTI {
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
    virtual ~LTI() {};
};


namespace InfInteraction{

// Project measured wrench to each joint torque.
class JointTorqueFromWrenchProjector: public LTI {
    OpenRAVE::RobotBasePtr robot_ptr;
    OpenRAVE::RobotBase::ManipulatorPtr ft_sensor_ptr;
    dVector jacobian, jacobian_rot, jacobian_T, jacobian_rot_T,  // Translational and rotational Jacobians, transposed
            force, torque,  // input force, torque
            tau1, tau2, tau, // projected torque
            jnt_pos_current;  // current joint position
    OpenRAVE::Transform T_wee;
    OpenRAVE::RaveVector<double> rave_force, rave_torque, temp_vec;
public:
    JointTorqueFromWrenchProjector(OpenRAVE::RobotBasePtr robot_ptr_, std::string ft_sensor_frame);
    dVector compute(const dVector & u_n);
    void set_state(const dVector & x_n);
};


// TODO: Current the projector assume that the robot's end-effector pose remains at the initial pose.
// This assumtption is not true for a more general class of problem. (But true for the Cartersian admittance
// task I am working on.

/*! Project measured external acting wrench to Cartesian workspace
 *
 */
class Wrench2CartForceProjector: public LTI {
    OpenRAVE::Transform T_wee;
public:
    Wrench2CartForceProjector(OpenRAVE::RobotBasePtr robot_ptr_, std::string ft_sensor_frame);
    dVector compute(const dVector & wrench);
    void set_state(const dVector & x_n);
};


/*A collection of controller. Used to implement Joint Admittance controller, where each
 * joint is controlled individually.
 */
class ControllerCollection: public LTI {
    std::vector<std::shared_ptr<LTI > > controllers;
public:
    ControllerCollection(std::vector<std::shared_ptr<LTI > > controllers_);
    dVector compute(const dVector & y_n);
    void set_state(const dVector & x_n);
};

// A Block that simply offset the input by an initially given vector.
class SimpleOffset: public LTI {
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
 * initial value quat_init. At each iteration, the robot's position is set to jnt_pos_current. Therefore, for
 * applications in which the robot moves, which basically are every applications, users should call set_state
 * and update jnt_pos_current before calling compute.
 */
 class CartPositionTracker: public LTI {
     dVector jnt_pos_current  /* Current robot joint position */ ;
     OpenRAVE::RaveVector<OpenRAVE::dReal>
             quat_init, /*Initial orientation*/
             pos_init /*Initial position*/ ;
     OpenRAVE::RobotBasePtr robot_ptr;
     OpenRAVE::RobotBase::ManipulatorPtr manip_ptr;
 public:
     CartPositionTracker(OpenRAVE::RobotBasePtr robot_ptr_, std::string manip_frame, dVector jnt_pos_init);

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
class FIRsrfb: public LTI {
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
 * The filter is constructed by two arrays: the numeration and the denominator.
 * In the following, x[n] is the input while y[n] is the output.
 */
class DiscreteTimeFilter: public LTI {
    /*! Order of the filter. Is basically the order of the denominator. */
    int order;
    /*! Current index. */
    int n, mem_sz;
    /*! The last (order)-th input and output is stored. All are default to some initial condition. */
    dVector x_mem, y_mem;
    dVector b, a;
public:
    /*! Initialize the filter.
        *
        * \param b_in coefficients of the numerator
        * \param a_in coefficients of the denominator
        * \param y_initial_val initial condition of the output. initial condition of input is default to zero.
        */
    DiscreteTimeFilter(dVector b_in, dVector a_in, double y_initial_val=0);
    /*! Compute the output in the next time step.
     *
     * ccde: a0 x[n] + a1 x[n-1] + ... a_o x[n-o]  = b0 y[n] + b1 y[n-1] + ... + b_o y[n - o]
     *
     * @param x_n
     * @return y_n
     */
    double compute(double x_n);
    dVector compute(const dVector & x_n);
    void set_state(const dVector & x_n);
};



/*! FT signal stream handler.
 *
 * A handler applies basic operations including collecting data and
 * filter them using DiscreteTimeFilter objects. Note that there are 6
 * individual filters, but all have the same coefficients.
 */
class FTSensorHandler {
    double fx, fy, fz, tx, ty, tz;
    std::vector<double> wrench_offset;
    std::vector<DiscreteTimeFilter> lp_filters;
    dVector b, a; /*!filter coefficients*/
    bool debug = false;
    ros::Publisher wrench_pub_debug;  // publisher for debugging
public:
    FTSensorHandler();;
    explicit FTSensorHandler(const dVector &wrench_offset_input);; // with non-zero offset values
    FTSensorHandler(const dVector &wrench_offset_input, dVector b_in, dVector a_in); // with offset and low-pass filter
    void signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg);
    void get_latest_wrench(dVector &force, dVector &torque);
    void get_latest_wrench(dVector &wrench);
    void set_debug(ros::NodeHandle &nh);
    void set_wrench_offset(dVector wrench_offset_);
    void log_latest_wrench(const std_msgs::Header &);
};
class JointPositionHandler{
    std::vector<double> joint_position;
public:
    JointPositionHandler();
    void signal_callback(const sensor_msgs::JointStateConstPtr &msg);
    dVector get_latest_jnt_position();
    bool received_msg();
};



/*! A controller used for controlling the robot joint position.
 *
 * */
class JointPositionController {
public:
    /*! Constructor for a joint position controller class.
     */
    explicit JointPositionController(std::string name_space, ros::NodeHandle& nh);
    void set_joint_positions(std::vector<double>& jnt_positions);

private:
    std::string _name_space;
    std::vector<ros::Publisher> _jnt_pubs;
    ros::NodeHandle _nh;

};

class ExternalTorquePublisher {
public:
    explicit ExternalTorquePublisher(std::string name_space, ros::NodeHandle& nh);
    void publish_joint_torques(std::vector<double>& taus);

private:
    std::string _name_space;
    std::vector<ros::Publisher> _jnt_torque_pubs;
    ros::NodeHandle _nh;
};

[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
void matrix_mult(std::vector<double> &A, std::vector<double> &x, std::vector<double> &y);
[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
void matrix_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
[[deprecated("Use Eigen3 instead: cleaner code and less buggy")]]
std::vector<double> mat_transpose(const std::vector<double>& M, int ncol);

#endif //PROJECT_INFINITE_INTERACTION_LIB_H
