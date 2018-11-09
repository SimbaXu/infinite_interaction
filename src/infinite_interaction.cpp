//
// Created by hung on 7/11/18.
//

#include <include/infinite_interaction_lib.h>

#include "infinite_interaction_lib.h"


namespace InfInteraction {
    JointTorqueFromWrenchProjector::JointTorqueFromWrenchProjector(OpenRAVE::RobotBasePtr robot_ptr_, std::string ft_sensor_frame)
            : robot_ptr(robot_ptr_)
    {

        ft_sensor_ptr = robot_ptr->GetManipulator(ft_sensor_frame);
        if (!ft_sensor_ptr){
            ROS_ERROR_STREAM("Unable to find ft_sensor_frame [" + ft_sensor_frame + "]");
            ros::shutdown();
        }
        force.resize(3);
        torque.resize(3);
        tau1.resize(6);
        tau2.resize(6);
        tau.resize(6);
    }
    dVector JointTorqueFromWrenchProjector::compute(const dVector &u_n) {
        // Jacobians are (3x6) matrix despite the robot having 7 dof because
        // there are only 6 joints from the base link to the ft sensor.
        ft_sensor_ptr->CalculateJacobian(jacobian);
        ft_sensor_ptr->CalculateAngularVelocityJacobian(jacobian_rot);
        T_wee = ft_sensor_ptr->GetEndEffectorTransform();

        // transform measured force/torque to world frame
        temp_vec.x = u_n[0]; temp_vec.y = u_n[1]; temp_vec.z = u_n[2];
        rave_force = T_wee.rotate(temp_vec);
        temp_vec.x = u_n[3]; temp_vec.y = u_n[4]; temp_vec.z = u_n[5];
        rave_torque = T_wee.rotate(temp_vec);

        // project to joint space
        force[0] = rave_force.x; force[1] = rave_force.y; force[2] = rave_force.z;
        torque[0] = rave_torque.x; torque[1] = rave_torque.y; torque[2] = rave_torque.z;
        jacobian_T = mat_transpose(jacobian, 6);
        jacobian_rot_T = mat_transpose(jacobian_rot, 6);
        matrix_mult(jacobian_T, force, tau1);
        matrix_mult(jacobian_rot_T, torque, tau2);
        matrix_add(tau1, tau2, tau);
        return tau;
    }
    void JointTorqueFromWrenchProjector::set_state(const dVector &x_n) {
        robot_ptr->SetActiveDOFValues(x_n);
    }

    // ControllerCollection
    ControllerCollection::ControllerCollection(std::vector<std::shared_ptr<LTI> > controllers_):
    controllers(controllers_) {}
    dVector ControllerCollection::compute(const dVector &y_n) {
        dVector u_n(6);
        for (unsigned int i = 0; i < 6; ++i) {
            // run each filter over the given input individually
            dVector y_ = controllers[i]->compute(dVector {y_n[i]});
            u_n[i] = y_[0];
        }
        return u_n;
    }
    void ControllerCollection::set_state(const dVector &x_n) {
        // Do nothing
    }

    dVector SimpleOffset::compute(const dVector &y_n) {
        dVector y_offset = y_n;
        for(unsigned int i=0; i < y_offset.size() && i < offset.size(); i++){
            y_offset[i] += offset[i];
        }
        return y_offset;
    }
    void SimpleOffset::set_state(const dVector &x_n) {
        // Do nothing
    }
    SimpleOffset::SimpleOffset(const dVector &offset_) {
        offset.resize(offset_.size());
        for(unsigned int i=0; i < offset.size(); i++){
            offset[i] = offset_[i];
        }
    }

    Wrench2CartForceProjector::Wrench2CartForceProjector(OpenRAVE::RobotBasePtr robot_ptr_,
                                                         std::string ft_sensor_frame) {
        T_wee = robot_ptr_->GetManipulator(ft_sensor_frame)->GetEndEffectorTransform();
    }
    dVector Wrench2CartForceProjector::compute(const dVector &wrench) {
        OpenRAVE::RaveVector<double> force (wrench[0], wrench[1], wrench[2]), force_rotated;
        force_rotated = T_wee.rotate(force);
        return dVector {force_rotated.x, force_rotated.y, force_rotated.z};
    }
    void Wrench2CartForceProjector::set_state(const dVector &x_n) {
        // Do nothing
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
    wrench.resize(6);
    wrench[0] = fx;
    wrench[1] = fy;
    wrench[2] = fz;
    wrench[3] = tx;
    wrench[4] = ty;
    wrench[5] = tz;
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

std::vector<double> mat_transpose(const std::vector<double>& M, unsigned int ncol){
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
