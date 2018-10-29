//
// Created by hung on 3/10/18.
//

#include "gtest/gtest.h"
#include <infinite_interaction_lib.h>


TEST(IIR, case1){
    dVector num {0}, den{1, -0.99};
    DiscreteTimeFilter filter (num, den, 10);
    double y;
    // testing
    y = filter.compute(123.321);  // random number
    EXPECT_NEAR(y, 9.9, 1e-3);
    y = filter.compute(3.321);
    EXPECT_NEAR(y, 9.9 * 0.99, 1e-3);
}

TEST(IIR, case2){
    // implement the filter y[n] - y[n-1] = x[n]
    dVector b{1}, a{1, -1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    y = filter.compute(1);
    EXPECT_EQ(y, 1);
    y = filter.compute(1.5);
    EXPECT_EQ(y, 2.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 0.0);
}


TEST(IIR, order_zero){
    // implement the filter y[n] = -x[n]
    dVector b{-1}, a{1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    y = filter.compute(1);
    EXPECT_EQ(y, -1);
    y = filter.compute(1.5);
    EXPECT_EQ(y, -1.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 2.5);
}
