//
// Created by hung on 11/11/18.
//

#include "infinite_interaction/infinite_interaction_lib.h"
#include <Eigen/Dense>
#include <qpOASES.hpp>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <include/infinite_interaction/infinite_interaction_lib.h>


InfInteraction::JointTorqueFromWrenchProjector::JointTorqueFromWrenchProjector(OpenRAVE::RobotBasePtr robot_ptr_, std::__cxx11::string ft_sensor_frame)
        : robot_ptr(robot_ptr_)
{

    InfInteraction::JointTorqueFromWrenchProjector::ft_sensor_ptr = InfInteraction::JointTorqueFromWrenchProjector::robot_ptr->GetManipulator(ft_sensor_frame);
    if (!InfInteraction::JointTorqueFromWrenchProjector::ft_sensor_ptr){
        ROS_ERROR_STREAM("Unable to find ft_sensor_frame [" + ft_sensor_frame + "]");
        ros::shutdown();
    }
    InfInteraction::JointTorqueFromWrenchProjector::force.resize(3);
    InfInteraction::JointTorqueFromWrenchProjector::torque.resize(3);
    InfInteraction::JointTorqueFromWrenchProjector::tau1.resize(6);
    InfInteraction::JointTorqueFromWrenchProjector::tau2.resize(6);
    InfInteraction::JointTorqueFromWrenchProjector::tau.resize(6);
}

dVector InfInteraction::JointTorqueFromWrenchProjector::compute(const dVector &u_n) {
    // lock environment, then change robot dof
    boost::recursive_mutex::scoped_lock lock(robot_ptr->GetEnv()->GetMutex());
    robot_ptr->SetActiveDOFValues(jnt_pos_current);
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

void InfInteraction::JointTorqueFromWrenchProjector::set_state(const dVector &x_n) {
    jnt_pos_current = x_n;
}

// ControllerCollection
InfInteraction::ControllerCollection::ControllerCollection(std::vector<std::shared_ptr<LTI> > controllers_):
controllers(controllers_) {}

dVector InfInteraction::ControllerCollection::compute(const dVector &y_n) {
    dVector u_n(6);
    for (unsigned int i = 0; i < 6; ++i) {
        // run each filter over the given input individually
        dVector y_ = controllers[i]->compute(dVector {y_n[i]});
        u_n[i] = y_[0];
    }
    return u_n;
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

InfInteraction::Wrench2CartForceProjector::Wrench2CartForceProjector(OpenRAVE::RobotBasePtr robot_ptr_,
                                                                     std::__cxx11::string ft_sensor_frame) {
    T_wee = robot_ptr_->GetManipulator(ft_sensor_frame)->GetEndEffectorTransform();
}

dVector InfInteraction::Wrench2CartForceProjector::compute(const dVector &wrench) {
    OpenRAVE::geometry::RaveVector<double> force (wrench[0], wrench[1], wrench[2]), force_rotated;
    force_rotated = T_wee.rotate(force);
    return dVector {force_rotated.x, force_rotated.y, force_rotated.z};
}

void InfInteraction::Wrench2CartForceProjector::set_state(const dVector &x_n) {
    // Do nothing
}

InfInteraction::CartPositionTracker::CartPositionTracker(OpenRAVE::RobotBasePtr robot_ptr_, std::__cxx11::string manip_frame,
                                                         dVector jnt_pos_init): robot_ptr(robot_ptr_), jnt_pos_current(jnt_pos_init) {
    // lock
    boost::recursive_mutex::scoped_lock lock(robot_ptr->GetEnv()->GetMutex());
    robot_ptr_->SetActiveDOFValues(jnt_pos_current);

    manip_ptr = robot_ptr->GetManipulator(manip_frame);
    if(!manip_ptr){
        ROS_ERROR_STREAM("Manipulator [" << manip_frame << "] not found!");
        ros::shutdown();
    }
    else {
        OpenRAVE::Transform T_wee = manip_ptr->GetTransform();
        quat_init = T_wee.rot;
        pos_init = T_wee.trans;
    }
}

dVector InfInteraction::CartPositionTracker::compute(const dVector &pos_n) {
    // lock
    boost::recursive_mutex::scoped_lock lock(robot_ptr->GetEnv()->GetMutex());
    robot_ptr->SetActiveDOFValues(jnt_pos_current);

    // get current position and quaternion
    OpenRAVE::Transform T_wee_cur = manip_ptr->GetTransform();
    OpenRAVE::geometry::RaveVector<double> pos_n_rave (pos_n[0], pos_n[1], pos_n[2]);
    // pos_cur = pos_init + pos_n
    auto dpos_rave = pos_init + pos_n_rave - T_wee_cur.trans;
    auto dquat_rave = quat_init - T_wee_cur.rot;
    dVector J_trans_arr, J_rot_arr;
    manip_ptr->CalculateJacobian(J_trans_arr);
    manip_ptr->CalculateRotationJacobian(J_rot_arr);

    // Eigen::Map provides a view on the contigous memory array. Very convenience.
    Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::RowMajor> > J_trans(J_trans_arr.data());
    Eigen::Map<Eigen::Matrix<double, 4, 6, Eigen::RowMajor> > J_rot(J_rot_arr.data());
    Eigen::Matrix<double, 3, 1> dpos; dpos << dpos_rave.x, dpos_rave.y, dpos_rave.z;
    Eigen::Matrix<double, 4, 1> dquat; dquat << dquat_rave.x, dquat_rave.y, dquat_rave.z, dquat_rave.w;


    // Solve an optimization to find the next best action:
    // min (dpos_rave - J_trans_arr dq)^2 + (dquat_rave - J_rot dq) ^ 2 + gam * dq^2 (for regularization)
    // s.t.    dqmin <= dq <= dqmax
    //         dqmin - q_cur <= dq <= qmax - q_cur
    //
    // NOTE: the following equality is every helpful:
    //   (dpos_rave - J_trans_arr dq)^2  = 0.5 dq.T (2 J^T J) dq - (2 J^T dpos) dq + dpos^T dpos
    // NOTE: There are 6 constraints, and 6 bounds on dq.

    // Form Quadratic objective: 0.5 x^T H x + g^T x
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> H;
    H = 2 * J_trans.transpose() * J_trans + 2 * J_rot.transpose() * J_rot + 2 * Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Identity() * 0.01;
    Eigen::Matrix<double, 6, 1> g;
    g = - 2 * J_trans.transpose() * dpos - 2 * J_rot.transpose() * dquat;

    // Form constraints: TODO: include constraints
    dVector dqmin(6), dqmax(6);
    std::fill(dqmin.begin(), dqmin.end(), -0.1);
    std::fill(dqmax.begin(), dqmax.end(),  0.1);

    // Form problem
    qpOASES::Options options;
    options.printLevel = qpOASES::PL_LOW;
    qpOASES::SQProblem qp_instance(6, 0);
    qp_instance.setOptions(options);
    int nWSR = 100;
    qpOASES::returnValue ret = qp_instance.init(H.data(), g.data(), NULL, dqmin.data(), dqmax.data(), NULL, NULL, nWSR, 0);
    dVector dq(6);
    if (ret == qpOASES::returnValue::SUCCESSFUL_RETURN){
        qpOASES::real_t dq_opt[6];
        qp_instance.getPrimalSolution(dq_opt);
        for(int i=0; i < 6; ++i) dq[i] = dq_opt[i];
    }
    else {
        ROS_WARN_STREAM("Optimization with qpOASES fails. Setting dq to zero.");
        std::fill(dq.begin(), dq.end(), 0);
    }
    dVector jnt_pos_cmd;
    for(unsigned int i=0; i < 6; i ++){
        jnt_pos_cmd.push_back(jnt_pos_current[i] + dq[i]);
    }
    return jnt_pos_cmd;
}

void InfInteraction::CartPositionTracker::set_state(const dVector &jnt_pos_n) {
    jnt_pos_current = jnt_pos_n;
}

