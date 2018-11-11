//
// Created by hung on 11/11/18.
//

#include "gtest/gtest.h"
#include <infinite_interaction/infinite_interaction_lib.h>


TEST(TopicDebugger, Basic){
    ros::NodeHandle nh;
    InfInteraction::TopicDebugger debugger("/debugging", nh);
    int t1Id = debugger.register_multiarray("topic1");
    int t2Id = debugger.register_multiarray("topic2");

    ros::Rate rate(10);
    for(int i=0; i < 10 and ros::ok(); i++){
        debugger.publish_multiarray(0, dVector {0, 1, 2, 3});
        debugger.publish_multiarray(1, dVector {4, 5, 6, 7});
        rate.sleep();
    }

    EXPECT_EQ(0, t1Id);
    EXPECT_EQ(1, t2Id);
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  ros::init(argc, argv, "tester");
  ros::NodeHandle nh;
  return RUN_ALL_TESTS();
}

