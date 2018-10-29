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

typedef std::vector<double> dVector;

/*! Discrete-time SISO filter.
 *
 * The filter is constructed by two arrays: the numeration and the denominator.
 * In the following, x[n] is the input while y[n] is the output.
 */
class DiscreteTimeFilter {
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
};

DiscreteTimeFilter::DiscreteTimeFilter(dVector b_in, dVector a_in, double y_initial_val) {
    n = 0;
    order = a_in.size() - 1;
    mem_sz = order + 1;
    assert(order >= 0); // at least a zeroth-order coefficient is expected
    assert(a_in[0] > 0);  // first coefficient needs to be positive, this might not be necessary though.

    x_mem.resize(mem_sz);
    y_mem.resize(mem_sz);
    for (int i = 0; i < mem_sz; ++i) {
        x_mem[i] = 0;
        y_mem[i] = y_initial_val;
    }

    // resize b_in to a vector of size (order + 1) and shift it rightward
    int size_a = b_in.size();
    for (int i=size_a; i <= order; ++i){
        b_in.push_back(0);
    }
    b = b_in;
    a = a_in;
    assert(b.size() == a.size());
}

double DiscreteTimeFilter::compute(double x_n) {
    n++;
    int deg = order + 1;
    x_mem[n % deg] = x_n;
    double sum = 0;
    for (int i = 0; i < deg; ++i) {
        sum += b[i] * x_mem[(n - i + deg) % deg];
    }
    for (int i = 1; i < deg; ++i){
        sum -= a[i] * y_mem[(n - i + deg) % deg];
    }
    y_mem[n % deg] = sum / a[0];
    return y_mem[n % deg];
}


/*! FT signal stream handler.
 *
 * A handler applies basic operations including collecting data and
 * filter them using the DiscreteTimeFilter class.
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
    FTSensorHandler(const dVector &wrench_offset_input, dVector b, dVector a); // with offset and low-pass filter
    void signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg);
    void get_latest_wrench(dVector &force, dVector &torque);
    void set_debug(ros::NodeHandle &nh);
    void log_latest_wrench(const std_msgs::Header &);
};

void FTSensorHandler::set_debug(ros::NodeHandle &nh) {
    debug=true;
    std::string topic_name = "debug/netft/filtered";
    wrench_pub_debug = nh.advertise<geometry_msgs::WrenchStamped>(topic_name, 5);
}

void FTSensorHandler::signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg) {
    fx = lp_filters[0].compute(msg->wrench.force.x - wrench_offset[0]);
    fy = lp_filters[1].compute(msg->wrench.force.y - wrench_offset[1]);
    fz = lp_filters[2].compute(msg->wrench.force.z - wrench_offset[2]);
    tx = lp_filters[3].compute(msg->wrench.torque.x - wrench_offset[3]);
    ty = lp_filters[4].compute(msg->wrench.torque.y - wrench_offset[4]);
    tz = lp_filters[5].compute(msg->wrench.torque.z - wrench_offset[5]);
    log_latest_wrench(msg->header);
}

void FTSensorHandler::get_latest_wrench(std::vector<double> &force, std::vector<double> &torque) {
    force.resize(3);
    torque.resize(3);
    force[0] = fx;
    force[1] = fy;
    force[2] = fz;
    torque[0] = tx;
    torque[1] = ty;
    torque[2] = tz;
}

FTSensorHandler::FTSensorHandler(const std::vector<double> &wrench_offset_input) :  fx(0), fy(0), fz(0), tx(0), ty(0), tz(0) {
    wrench_offset.resize(6);
    b = {1};
    a = {1};
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = wrench_offset_input[i];
        lp_filters.push_back(DiscreteTimeFilter(b, a, 0));  // set initial output to be the offset val
    }
}

FTSensorHandler::FTSensorHandler() : fx(0), fy(0), fz(0), tx(0), ty(0), tz(0) {
    wrench_offset.resize(6);
    b = {1};
    a = {1};
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = 0;
        lp_filters.push_back(DiscreteTimeFilter(b, a, 0));
    }
}

FTSensorHandler::FTSensorHandler(const dVector &wrench_offset_input, dVector b_in, dVector a_in) {
    b = b_in;
    a = a_in;
    wrench_offset.resize(6);
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = wrench_offset_input[i];
        lp_filters.push_back(DiscreteTimeFilter(b, a, 0));
    }
}

void FTSensorHandler::log_latest_wrench(const std_msgs::Header &header) {
    if (debug){
        geometry_msgs::WrenchStamped msg;
        msg.header = header;
        msg.wrench.force.x = fx;
        msg.wrench.force.y = fy;
        msg.wrench.force.z = fz;
        msg.wrench.torque.x = tx;
        msg.wrench.torque.y = ty;
        msg.wrench.torque.z = tz;
        wrench_pub_debug.publish(msg);
    }
}

class JointPositionHandler{
    std::vector<double> joint_position;
public:
    JointPositionHandler();
    void signal_callback(const sensor_msgs::JointStateConstPtr &msg);
    dVector get_latest_jnt_position();
    bool received_msg();
};

JointPositionHandler::JointPositionHandler() {
    joint_position.resize(6);
    for (int i = 0; i < 6; ++i) {
        joint_position[i] = -999; // coded, should changed if have received message
    }
}

void JointPositionHandler::signal_callback(const sensor_msgs::JointStateConstPtr &msg) {
    for (int i = 0; i < 6; ++i) {
        joint_position[i] = msg->position[i];
    }
}

dVector JointPositionHandler::get_latest_jnt_position() {
    return joint_position;
}

bool JointPositionHandler::received_msg() {
    for (int i=0; i < 6; ++i){
        if (joint_position[i] == -999){
            return false;
        }
    }
    return true;
}

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

JointPositionController::JointPositionController(std::string name_space, ros::NodeHandle &nh) : _name_space(name_space), _nh(nh) {
    for (int i=0; i < 6; i++){
        std::string jnt_topic = "/" + _name_space + "/j" + std::to_string(i + 1) + "/command";
        ros::Publisher jnt_pub = _nh.advertise<std_msgs::Float64>(jnt_topic, 5);
        _jnt_pubs.push_back(jnt_pub);
    }
}

void JointPositionController::set_joint_positions(std::vector<double>& jnt_positions) {
    for(int i=0; i < 6; i++){
        std_msgs::Float64 msg;
        msg.data = jnt_positions[i];
        _jnt_pubs[i].publish(msg);
    }
}


class ExternalTorquePublisher {
public:
    explicit ExternalTorquePublisher(std::string name_space, ros::NodeHandle& nh);
    void publish_joint_torques(std::vector<double>& taus);

private:
    std::string _name_space;
    std::vector<ros::Publisher> _jnt_torque_pubs;
    ros::NodeHandle _nh;
};

ExternalTorquePublisher::ExternalTorquePublisher(std::string name_space, ros::NodeHandle &nh): _name_space(name_space), _nh(nh)  {
    for (int i=0; i < 6; i++){
        std::string jnt_torque_topic = _name_space + "/tau" + std::to_string(i + 1);
        ros::Publisher jnt_torque_pub = _nh.advertise<std_msgs::Float64>(jnt_torque_topic, 5);
        _jnt_torque_pubs.push_back(jnt_torque_pub);
    }
}

void ExternalTorquePublisher::publish_joint_torques(std::vector<double> &taus) {
    for(int i=0; i < 6; i++){
        std_msgs::Float64 msg;
        msg.data = taus[i];
        _jnt_torque_pubs[i].publish(msg);
    }
}


void matrix_mult(std::vector<double> &A, std::vector<double> &x, std::vector<double> &y){
    unsigned int size_y = A.size() / x.size(), size_x = x.size();
    y.resize(size_y);
    assert(A.size() % x.size() == 0);
    for (int i = 0; i < size_y; ++i) {
        y[i] = 0;
        for (int j = 0; j < size_x; ++j) {
            y[i] += A[i * size_x + j] * x[j];
        }
    }
}

void matrix_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    assert(x.size() == y.size());
    z.resize(x.size());
    for (int i = 0; i < x.size(); ++i) {
        z[i] = x[i] + y[i];
    }
}

std::vector<double> mat_transpose(const std::vector<double>& M, int ncol){
    assert(M.size() % ncol == 0);
    std::vector<double> M_T (M.size());
    int nrow = static_cast<int>(M.size() / ncol);
    for(int i=0; i < nrow; i++){
        for(int j=0; j < ncol; j++){
            M_T[j * nrow + i] = M[i * ncol + j];
        }
    }
    return M_T;
}

#endif //PROJECT_INFINITE_INTERACTION_LIB_H
