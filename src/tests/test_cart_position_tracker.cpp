//
// Created by hung on 10/11/18.
//

#include "gtest/gtest.h"
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

void print_6dvec(dVector v){
    std::cout
            << v[0] << ", "
            << v[1] << ", "
            << v[2] << ", "
            << v[3] << ", "
            << v[4] << ", "
            << v[5] << ", " << std::endl;
}

TEST(CartTracker, visual){
    std::string scenefilename = "robots/denso_ft_sensor_gripper.robot.xml";
    std::string robot_name = "denso_ft_sensor_gripper";
    std::string manip_name = "gripper";
    // Create an OpenRAVE instance for kinematic computations (Jacobian and stuffs)
    OpenRAVE::RaveInitialize(true); // start openrave core
    OpenRAVE::EnvironmentBasePtr env_ptr = OpenRAVE::RaveCreateEnvironment(); // create the main environment
    OpenRAVE::RaveSetDebugLevel(OpenRAVE::Level_Info);
//    boost::thread thviewer(boost::bind(SetViewer,env_ptr, "qtosg"));  // create viewer
    env_ptr->Load(scenefilename); // load the scene
    OpenRAVE::RobotBasePtr robot_ptr;
    robot_ptr = env_ptr->GetRobot(robot_name);
    robot_ptr->SetActiveDOFs(std::vector<int> {0, 1, 2, 3, 4, 5});
    auto manip_ptr = robot_ptr->GetManipulator(manip_name);

    // Fed the same position 100 times
    std::vector<dVector > cart_poss; // desired cart positions
    double dx = -0.1, dy = 0.05, dz=-0.035;
    for (int i=0; i < 20; ++i){
        cart_poss.push_back(dVector {dx, dy, dz});
    }

    // Desired value
    dVector jnt_pos_init = {0, 0.57, .9, 0.2, 0, 0}; // initial position
    robot_ptr->SetActiveDOFValues(jnt_pos_init);
    auto T_wee_init = manip_ptr->GetTransform();
    auto pose_init = T_wee_init.rot;
    robot_ptr->SetActiveDOFValues(dVector {0, 0, 0, 0, 0, 0}); // reset the robot to the zero psotion

    // Cartesian Tracker
    std::shared_ptr<InfInteraction::CartPositionTracker> position_map_ptr =
            std::make_shared<InfInteraction::CartPositionTracker>(robot_ptr, manip_name, jnt_pos_init);
    dVector jnt_pos_cur = jnt_pos_init;
    for (int i=0; i < cart_poss.size(); ++i){
        position_map_ptr->set_state(jnt_pos_cur);
        jnt_pos_cur = position_map_ptr->compute(cart_poss[i]);
    }
    // Test
    robot_ptr->SetActiveDOFValues(jnt_pos_cur);
    auto T_wee_final = manip_ptr->GetTransform();

    double eps = 1e-3, eps_quat=1e-2;

    EXPECT_NEAR(T_wee_init.trans.x + dx, T_wee_final.trans.x, eps) << "pos[x] wrong";
    EXPECT_NEAR(T_wee_init.trans.y + dy, T_wee_final.trans.y, eps) << "pos[y] wrong";
    EXPECT_NEAR(T_wee_init.trans.z + dz, T_wee_final.trans.z, eps) << "pos[z] wrong";
    EXPECT_NEAR(pose_init.x, T_wee_final.rot.x, eps_quat) << "quat[x] wrong";
    EXPECT_NEAR(pose_init.y, T_wee_final.rot.y, eps_quat) << "quat[y] wrong";
    EXPECT_NEAR(pose_init.z, T_wee_final.rot.z, eps_quat) << "quat[z] wrong";
    EXPECT_NEAR(pose_init.w, T_wee_final.rot.w, eps_quat) << "quat[w] wrong";

//    thviewer.join(); // wait for the viewer thread to exit
    env_ptr->Destroy(); // destroy
}

