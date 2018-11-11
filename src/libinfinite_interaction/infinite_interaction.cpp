//
// Created by hung on 7/11/18.
//

#include "infinite_interaction/infinite_interaction_lib.h"
#include <ros/ros.h>

InfInteraction::TopicDebugger::TopicDebugger(std::string debug_ns_, ros::NodeHandle &node_handle_): debug_ns(debug_ns_), nh(node_handle_){ }

int InfInteraction::TopicDebugger::register_multiarray(std::string topic_name) {
    ros::Publisher pub_ = nh.advertise<MultiArrayMsg>(debug_ns + "/" + topic_name, 10);
    pub_vecs.push_back(std::make_shared<ros::Publisher> (std::move(pub_)));
    int topic_id = pub_vecs.size() - 1;
    return topic_id;
}

void InfInteraction::TopicDebugger::publish_multiarray(int topic_id, const dVector &data) {
    std_msgs::MultiArrayDimension dim;
    dim.label = "x";
    dim.size = static_cast<unsigned int>(data.size());
    dim.stride = 1;
    MultiArrayMsg msg;
    msg.layout.dim.push_back(dim);

    if (topic_id >= 0 and topic_id < pub_vecs.size()){
        for (double x : data) msg.data.push_back(x);
        pub_vecs[topic_id]->publish(msg);
    }
    else {
        ROS_WARN_STREAM("Not debugging topic found!");
    }
}

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

void FTSensorHandler::get_latest_wrench(std::vector<double> &wrench) {
    if (fx == 0.12543214218){
        ROS_FATAL_STREAM("FT handler has not received any signal. Stopping!");
        ros::shutdown();
    }
    wrench.resize(6);
    wrench[0] = fx;
    wrench[1] = fy;
    wrench[2] = fz;
    wrench[3] = tx;
    wrench[4] = ty;
    wrench[5] = tz;
}

FTSensorHandler::FTSensorHandler(const std::vector<double> &wrench_offset_input) :  fx(0.12543214218), fy(0), fz(0), tx(0), ty(0), tz(0) {
    // fx is initialized with a special constant, which is checked to make sure that the handler has received wrench reading.
    wrench_offset.resize(6);
    b = {1};
    a = {1};
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = wrench_offset_input[i];
        lp_filters.push_back(DiscreteTimeFilter(b, a, 0));  // set initial output to be the offset val
    }
}

FTSensorHandler::FTSensorHandler() : fx(0.12543214218), fy(0), fz(0), tx(0), ty(0), tz(0) {
    wrench_offset.resize(6);
    b = {1};
    a = {1};
    for (int i = 0; i < 6; ++i) {
        wrench_offset[i] = 0;
        lp_filters.push_back(DiscreteTimeFilter(b, a, 0));
    }
}

FTSensorHandler::FTSensorHandler(const dVector &wrench_offset_input, dVector b_in, dVector a_in): fx(0.12543214218), fy(0), fz(0), tx(0), ty(0), tz(0)  {
    b = b_in;
    a = a_in;
    wrench_offset.resize(6);
    for (unsigned int i = 0; i < 6; ++i) {
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


void FTSensorHandler::set_wrench_offset(dVector wrench_offset_) {
    wrench_offset_.resize(6);
     for (unsigned int i = 0; i < 6; ++i) {
        wrench_offset[i] = wrench_offset_[i];
    }
}


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



DiscreteTimeFilter::DiscreteTimeFilter(dVector b_in, dVector a_in, double y_initial_val) {
    n = 0;
    order = static_cast<int>(a_in.size() - 1);
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
    int size_a = static_cast<int>(b_in.size());
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

dVector DiscreteTimeFilter::compute(const dVector & x_n) {
    dVector y_n {compute(x_n[0])};
    return y_n;
}


void DiscreteTimeFilter::set_state(const dVector & x_n) {
    // do nothing
};


void FIRsrfb::init_filter(dVector uinit) {
    // check shape consistency
    if (T * ny * nu != L.size() or T * nu * nu != MB2.size()){
        throw std::invalid_argument("Wrong input shape");
    }
    // the 0-th impulse of MB2 is zero
    for (unsigned int i = 0; i < nu * nu; ++i) {
        if (MB2[i] != 0){
            throw std::invalid_argument("MB2[0] needs to be zero");
        }
    }
    // initiate ybank and ubank
    ybank.resize(T * ny);
    ubank.resize(T * nu);
    ybank_sz = ybank.size();
    ubank_sz = ubank.size();
    std::fill(ybank.begin(), ybank.end(), 0);
    std::fill(ubank.begin(), ubank.end(), 0);
    // fill ubank with initial value
    for(unsigned int i=0; i < T; ++i){
        for(unsigned int j=0; j < nu; ++j){
            ubank[j + i * nu] = uinit[j];
        }
    }
    // initiate temporary variable
    alpha.resize(nu);
    beta.resize(nu);
    un.resize(nu);
}


FIRsrfb::FIRsrfb(unsigned int T_, unsigned int ny_, unsigned int nu_, dVector L_, dVector MB2_, dVector uinit) :  T(T_), ny(ny_), nu(nu_), L(L_), MB2(MB2_) {
    init_filter(uinit);
}

FIRsrfb::FIRsrfb(unsigned int T_, unsigned int ny_, unsigned int nu_, dVector L_, dVector MB2_) :  T(T_), ny(ny_), nu(nu_), L(L_), MB2(MB2_) {
    dVector uinit (T_ * nu_);
    std::fill(uinit.begin(), uinit.end(), 0);
    init_filter(uinit);
}

dVector FIRsrfb::compute(const dVector & y_n){
    // store y_n to ybank
    for (unsigned int i = 0; i < ny; ++i) {
        yidx = (yidx - 1 + ybank_sz) % ybank_sz;
        ybank[yidx] = y_n[ny - 1 - i];
    }
    // compute alpha,
    std::fill(alpha.begin(), alpha.end(), 0);
    // (i,j,k) represent three indices of matrix L
    // since L, shaped (T, nu, ny) is stored row-order, id = k + j * ny + i * nu * ny
    for (unsigned int i=0; i < T; ++i){
        for (unsigned int j=0; j < nu; ++j){
            for (unsigned int k=0; k < ny; ++k){
                alpha[j] += L[k + j * ny + i * nu * ny] * ybank[(yidx + k + i * ny) % ybank_sz];
            }
        }
    }
    // beta
    std::fill(beta.begin(), beta.end(), 0);
    uidx = (uidx - nu + ubank_sz) % ubank_sz;  // shift current index on ubank
    for (unsigned int i=1; i < T; ++i){
        for (unsigned int j=0; j < nu; ++j){
            for (unsigned int k=0; k < nu; ++k){
                beta[j] += MB2[k + j * nu + i * nu * nu] * ubank[(uidx + k + i * nu) % ubank_sz];
            }
        }
    }
    // and finally u_n,
    for(unsigned int i=0; i < nu; ++i){
        un[i] = alpha[i] - beta[i];
        // store u_n to ubank
        ubank[(uidx + i) % ubank_sz] = un[i];
    }

    return un;
}


void FIRsrfb::set_state(const dVector & x_n) {
    // do nothing
};

// y = Ax
void matrix_mult(std::vector<double> &A, std::vector<double> &x, std::vector<double> &y){
    unsigned int size_y = A.size() / x.size(), size_x = x.size();
    y.resize(size_y);
    assert(A.size() % x.size() == 0);
    for (unsigned int i = 0; i < size_y; ++i) {
        y[i] = 0;
        for (unsigned  int j = 0; j < size_x; ++j) {
            y[i] += A[i * size_x + j] * x[j];
        }
    }
}

// z = x + y
void matrix_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    assert(x.size() == y.size());
    z.resize(x.size());
    for (unsigned int i = 0; i < x.size(); ++i) {
        z[i] = x[i] + y[i];
    }
}

std::vector<double> mat_transpose(const std::vector<double>& M, int ncol){
    assert(M.size() % ncol == 0);
    std::vector<double> M_T (M.size());
    auto nrow = static_cast<unsigned int>(M.size() / ncol);
    for(unsigned int i=0; i < nrow; i++){
        for(unsigned int j=0; j < ncol; j++){
            M_T[j * nrow + i] = M[i * ncol + j];
        }
    }
    return M_T;
}
