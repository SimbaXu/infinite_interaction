//
// Created by hung on 8/9/18.
//

#ifndef PROJECT_INFINITE_INTERACTION_LIB_H
#define PROJECT_INFINITE_INTERACTION_LIB_H

#include <ros/ros.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float64.h>
#include <sensor_msgs/JointState.h>
#include <vector>

typedef std::vector<double> dVector;

class DiscreteTimeFilter {

    /*! Order of the filter. Is basically the order of the denominator. */
    int order;
    /*! Current index. */
    int n;
    /*! The last (order)-th input and output is stored. All are default to some initial condition. */
    dVector x_mem, y_mem;
    dVector cof_a, cof_b;
public:
    /*! Initialize the filter.
        *
        * \param cof_a_in numerator
        * \param cof_b_in denominator
        * \param y_initial_val initial condition of the output. initial condition of input is default to zero.
        */
    DiscreteTimeFilter(dVector &cof_a_in, dVector &cof_b_in, double y_initial_val=0);
    /*! Compute the output in the next time step.
     *
     * ccde: a0 x[n] + a1 x[n-1] + ... a_o x[n-o]  = b0 y[n] + b1 y[n-1] + ... + b_o y[n - o]
     *
     * @param x_n
     * @return y_n
     */
    double compute(double x_n);
};

DiscreteTimeFilter::DiscreteTimeFilter(dVector &cof_a_in, dVector &cof_b_in, double y_initial_val) {
    n = 0;
    order = cof_b_in.size() - 1;
    assert(order >= 0 and cof_b_in[0] > 0);

    x_mem.resize(order + 1);
    y_mem.resize(order + 1);
    for (int i = 0; i < order + 1; ++i) {
        x_mem[i] = 0;
        y_mem[i] = y_initial_val;
    }

    // resize cof_a_in to a vector of size (order + 1) and shift it rightward
    int order_a = cof_a_in.size() - 1;
    cof_a_in.resize(order + 1);
    for (int i=order; i >= 0; --i){
        if (i - (order - order_a) >= 0){
            cof_a_in[i] = cof_a_in[i - (order - order_a)];
        }
        else {
            cof_a_in[i] = 0;
        }
    }
    cof_a = cof_a_in;
    cof_b = cof_b_in;
}

double DiscreteTimeFilter::compute(double x_n) {
    n++;
    int deg = order + 1;
    x_mem[n % deg] = x_n;
    double sum = 0;
    for (int i = 0; i < deg; ++i) {
        sum += cof_a[i] * x_mem[(n - i + deg) % deg];
    }
    for (int i = 1; i < deg; ++i){
        sum -= cof_b[i] * y_mem[(n - i + deg) % deg];
    }
    y_mem[n % deg] = sum / cof_b[0];
    return y_mem[n % deg];
}


/*! A class for handling FT signal stream.
 *
 */
class FTSensorHandler {
    double fx, fy, fz, tx, ty, tz;
    std::vector<double> wrench_offset;
public:
    FTSensorHandler();;
    /*! Construct a FT Sensor Handler object.
     *
     * abc
     */
    explicit FTSensorHandler(const std::vector<double> &wrench_offset_input);;
    void signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg);
    void get_latest_wrench(std::vector<double> &force, std::vector<double> &torque);
//    void set_filter(); not implemented
};


void FTSensorHandler::signal_callback(const geometry_msgs::WrenchStampedConstPtr &msg) {
    fx = msg->wrench.force.x - wrench_offset[0];
    fy = msg->wrench.force.y - wrench_offset[1];
    fz = msg->wrench.force.z - wrench_offset[2];
    tx = msg->wrench.torque.x - wrench_offset[3];
    ty = msg->wrench.torque.y - wrench_offset[4];
    tz = msg->wrench.torque.z - wrench_offset[5];
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
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = wrench_offset_input[i];
    }
}

FTSensorHandler::FTSensorHandler() : fx(0), fy(0), fz(0), tx(0), ty(0), tz(0) {
    wrench_offset.resize(6);
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = 0;
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
