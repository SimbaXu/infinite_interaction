//
// Created by hung on 11/11/18.
//

#include "infinite_interaction/infinite_interaction_lib.h"
#include <Eigen/Dense>
#include <qpOASES.hpp>
#include <eigen3/Eigen/src/Core/Matrix.h>


InfInteraction::JointTorqueFromWrenchProjector::JointTorqueFromWrenchProjector(OpenRAVE::RobotBasePtr robot_ptr, std::__cxx11::string ft_sensor_frame)
        : robot_ptr_(robot_ptr)
{

    InfInteraction::JointTorqueFromWrenchProjector::ft_sensor_ptr_ = InfInteraction::JointTorqueFromWrenchProjector::robot_ptr_->GetManipulator(ft_sensor_frame);
    if (!InfInteraction::JointTorqueFromWrenchProjector::ft_sensor_ptr_){
        ROS_ERROR_STREAM("Unable to find ft_sensor_frame [" + ft_sensor_frame + "]");
        ros::shutdown();
    }
    InfInteraction::JointTorqueFromWrenchProjector::force.resize(3);
    InfInteraction::JointTorqueFromWrenchProjector::torque.resize(3);
    InfInteraction::JointTorqueFromWrenchProjector::tau1.resize(6);
    InfInteraction::JointTorqueFromWrenchProjector::tau2.resize(6);
    InfInteraction::JointTorqueFromWrenchProjector::tau.resize(6);
}

dVector InfInteraction::JointTorqueFromWrenchProjector::compute(const dVector &wrench_measure) {
    // lock environment, then change robot dof
    boost::recursive_mutex::scoped_lock lock(robot_ptr_->GetEnv()->GetMutex());
    robot_ptr_->SetActiveDOFValues(joint_current_);
    // Jacobians are (3x6) matrix despite the robot having 7 dof because
    // there are only 6 joints from the base link to the ft sensor.
    ft_sensor_ptr_->CalculateJacobian(jacobian);
    ft_sensor_ptr_->CalculateAngularVelocityJacobian(jacobian_rot);
    T_wee_ = ft_sensor_ptr_->GetEndEffectorTransform();

    // transform measured force_/torque to world frame
    temp_vec_.x = wrench_measure[0]; temp_vec_.y = wrench_measure[1]; temp_vec_.z = wrench_measure[2];
    rave_force_ = T_wee_.rotate(temp_vec_);
    temp_vec_.x = wrench_measure[3]; temp_vec_.y = wrench_measure[4]; temp_vec_.z = wrench_measure[5];
    rave_torque_ = T_wee_.rotate(temp_vec_);

    // project to joint space
    force[0] = rave_force_.x; force[1] = rave_force_.y; force[2] = rave_force_.z;
    torque[0] = rave_torque_.x; torque[1] = rave_torque_.y; torque[2] = rave_torque_.z;
    jacobian_T = mat_transpose(jacobian, 6);
    jacobian_rot_T = mat_transpose(jacobian_rot, 6);
    matrix_mult(jacobian_T, force, tau1);
    matrix_mult(jacobian_rot_T, torque, tau2);
    matrix_add(tau1, tau2, tau);
    return tau;
}

void InfInteraction::JointTorqueFromWrenchProjector::set_state(const dVector &joint_measure) {
    joint_current_ = joint_measure;
}

// ControllerCollection
InfInteraction::ControllerCollection::ControllerCollection(std::vector<std::shared_ptr<SignalBlock> > controllers_):
controllers(controllers_) {
    nb_controllers = static_cast<unsigned int>(controllers.size());
    for(int i=0; i < nb_controllers; i++){
        if (controllers[i]->get_input_size() != 1 or controllers[i]->get_output_size() != 1){
            throw "Controller have more than 1 inputs or outputs";
        };
    }
}

dVector InfInteraction::ControllerCollection::compute(const dVector &x_n) {
    dVector y_n(nb_controllers), y_ ;
    for (unsigned int i = 0; i < nb_controllers; ++i) {
        // run each filter over the given input individually
        y_ = controllers[i]->compute(dVector {x_n[i]});
        y_n[i] = y_[0];
    }
    return y_n;
}

void InfInteraction::ControllerCollection::set_state(const dVector &x_n) {
    // Do nothing
}

InfInteraction::SimpleOffset::SimpleOffset(const dVector &offset_) {
    offset.resize(offset_.size());
    for(unsigned int i=0; i < offset.size(); i++){
        offset[i] = offset_[i];
    }
}

dVector InfInteraction::SimpleOffset::compute(const dVector &y_n) {
    dVector y_offset = y_n;
    for(unsigned int i=0; i < y_offset.size() && i < offset.size(); i++){
        y_offset[i] += offset[i];
    }
    return y_offset;
}

void InfInteraction::SimpleOffset::set_state(const dVector &x_n) {
    // Do nothing
}

InfInteraction::Wrench2CartForceProjector::Wrench2CartForceProjector(OpenRAVE::RobotBasePtr robot_ptr,
                                                                     std::__cxx11::string ft_frame_name) {
    T_wee_ = robot_ptr->GetManipulator(ft_frame_name)->GetEndEffectorTransform();
}

dVector InfInteraction::Wrench2CartForceProjector::compute(const dVector &wrench_measure) {
    force_.x = wrench_measure[0];
    force_.y = wrench_measure[1];
    force_.z = wrench_measure[2];

    force_rotated_ = T_wee_.rotate(force_);
    return dVector {force_rotated_.x, force_rotated_.y, force_rotated_.z};
}

void InfInteraction::Wrench2CartForceProjector::set_state(const dVector &dummy_var) {
    // Do nothing
}

InfInteraction::CartPositionTracker::CartPositionTracker(OpenRAVE::RobotBasePtr robot_ptr_, std::__cxx11::string manip_frame, dVector jnt_pos_init, double gam, double gam2):
        robot_ptr(robot_ptr_), _jnt_pos_current(jnt_pos_init), _jnt_pos_init(jnt_pos_init), _gam2(gam2), _gam(gam)
{
    // lock
    boost::recursive_mutex::scoped_lock lock(robot_ptr->GetEnv()->GetMutex());
    robot_ptr_->SetActiveDOFValues(_jnt_pos_current);

    manip_ptr = robot_ptr->GetManipulator(manip_frame);
    if(!manip_ptr){
        ROS_ERROR_STREAM("Manipulator [" << manip_frame << "] not found!");
        ros::shutdown();
    }
    else {
        OpenRAVE::Transform T_wee = manip_ptr->GetTransform();
        _quat_init = T_wee.rot;
        _pos_init = T_wee.trans;
    }
}

dVector InfInteraction::CartPositionTracker::compute(const dVector &pos_n) {
    // lock
    boost::recursive_mutex::scoped_lock lock(robot_ptr->GetEnv()->GetMutex());
    robot_ptr->SetActiveDOFValues(_jnt_pos_current);

    // get current position and quaternion
    OpenRAVE::Transform T_wee_cur = manip_ptr->GetTransform();
    OpenRAVE::geometry::RaveVector<double> pos_n_rave (pos_n[0], pos_n[1], pos_n[2]);
    // pos_cur = _pos_init + pos_n
    auto dpos_rave = _pos_init + pos_n_rave - T_wee_cur.trans;
    auto dquat_rave = _quat_init - T_wee_cur.rot;
    dVector J_trans_arr, J_rot_arr;
    manip_ptr->CalculateJacobian(J_trans_arr);
    manip_ptr->CalculateRotationJacobian(J_rot_arr);

    // Eigen::Map provides a view on the contigous memory array. Very convenience.
    Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::RowMajor> > J_trans(J_trans_arr.data());
    Eigen::Map<Eigen::Matrix<double, 4, 6, Eigen::RowMajor> > J_rot(J_rot_arr.data());
    Eigen::Matrix<double, 3, 1> dpos; dpos << dpos_rave.x, dpos_rave.y, dpos_rave.z;
    Eigen::Matrix<double, 4, 1> dquat; dquat << dquat_rave.x, dquat_rave.y, dquat_rave.z, dquat_rave.w;


    // Solve an optimization to find the next best action:
    // min (dpos_rave - J_trans_arr dq)^2 + (dquat_rave - J_rot dq) ^ 2 + gam * dq^2 (t3) + (q_init - q_cur - dq)^2 * gam2 (t4)
    // s.t.    dqmin <= dq <= dqmax
    //         dqmin - q_cur <= dq <= qmax - q_cur
    //
    // NOTE: the following equalities are helpful in deriving coefficients for qpOASES solver
    //
    //   (dpos_rave - J_trans_arr dq)^2  = 0.5 dq.T (2 J^T J) dq - (2 J^T dpos) dq + dpos^T dpos
    //
    //   (q_i - q_c - dq)^2 = dq^T dq - 2 (q_i - q_c)^T dq + (q_i - q_c)^2
    //
    // NOTE: There are 6 constraints, and 6 bounds on dq.
    // Coefficient t3 is for regularization. gam is set to 0.01
    //
    // Coefficient t4 is to encourgate the robot to keep its
    // initial configuration. A side-effect of this term is that the robot would try to move back
    // to the initial configuration when there is no force_ acting on it. gam2 is set to 0.01

    // Form Quadratic objective: 0.5 x^T H x + g^T x
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> H;
    H = 2 * J_trans.transpose() * J_trans + 2 * J_rot.transpose() * J_rot + 2 * Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Identity() * (_gam + _gam2);
    Eigen::Matrix<double, 6, 1> g;
    Eigen::Map<Eigen::Matrix<double, 6, 1> > eigen_jnt_pos_init (_jnt_pos_init.data());
    Eigen::Map<Eigen::Matrix<double, 6, 1> > eigen_jnt_pos_current (_jnt_pos_current.data());
    g = - 2 * J_trans.transpose() * dpos - 2 * J_rot.transpose() * dquat - 2 * _gam2 * (eigen_jnt_pos_init - eigen_jnt_pos_current);

    // Form fixed bounds -0.1 <= dq <= 0.1 and joint limits constraints
    dVector dqmin(6), dqmax(6), qlow, qhigh;
    robot_ptr->GetActiveDOFLimits(qlow, qhigh);
    for(int i=0; i < 6; i++){
        dqmin[i] = MAX(qlow[i] - _jnt_pos_current[i], -0.1);
        dqmax[i] = MIN(qhigh[i] - _jnt_pos_current[i], 0.1);
    }

    // Form optimization problem and solve
    qpOASES::Options options;
    options.printLevel = qpOASES::PL_LOW;
    qpOASES::SQProblem qp_instance(6, 0);
    qp_instance.setOptions(options);
    int nWSR = 100;
    qpOASES::returnValue ret = qp_instance.init(H.data(), g.data(), NULL, dqmin.data(), dqmax.data(), NULL, NULL, nWSR, 0);
    dVector dq(6);

    // check result and set joint command
    dVector jnt_pos_cmd;
    if (ret == qpOASES::returnValue::SUCCESSFUL_RETURN){
        qpOASES::real_t dq_opt[6];
        qp_instance.getPrimalSolution(dq_opt);
        for(int i=0; i < 6; ++i) dq[i] = dq_opt[i];
        for(unsigned int i=0; i < 6; i ++){
            jnt_pos_cmd.push_back(_jnt_pos_current[i] + dq[i]);
        }
        // joint limits are satisfied by the optimization parameters

        // check for collision

    }
    else {
        ROS_WARN_STREAM("[Optimization fails] Keeping the robot at its current position.");
        for(unsigned int i=0; i < 6; i ++){
            jnt_pos_cmd.push_back(_jnt_pos_current[i]);
        }
    }

    return jnt_pos_cmd;
}

void InfInteraction::CartPositionTracker::set_state(const dVector &jnt_pos_n) {
    _jnt_pos_current = jnt_pos_n;
}


dVector InfInteraction::DelayFeedback::compute(const dVector &x_n) {
    dVector alpha, beta(_input_size), y_n;
    alpha = _fb_block->compute(y_last);
    for(int i; i < _input_size; i++){
        beta[i] = x_n[i] + _sign * alpha[i];
    }
    y_n = _fwd_block->compute(beta);
    y_last = y_n;
    return y_n;
}

InfInteraction::DelayFeedback::DelayFeedback(std::shared_ptr<SignalBlock> fwd_block,
                                             std::shared_ptr<SignalBlock> fb_block, int sign): _sign(sign) {
    _fwd_block = std::move(fwd_block);
    _fb_block = std::move(fb_block);
    y_last.resize(_fwd_block->get_output_size());
    _input_size = _fwd_block->get_output_size();

    std::fill(y_last.begin(), y_last.end(), 0);

    // check consistent input/output dimensions
    if(_fwd_block->get_output_size() != _fwd_block->get_input_size() or _fb_block->get_input_size() != _fb_block->get_output_size() or _fwd_block->get_output_size() != _fb_block->get_input_size()){
        throw "Inconsistent input/output dimensions";
    }
}
