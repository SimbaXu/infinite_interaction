//
// Created by hung on 3/10/18.
//

#include "gtest/gtest.h"
#include <infinite_interaction/infinite_interaction_lib.h>


/*! Script to generate the desired result
 *
import control as co
ff = co.tf([0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0], 1)
fb = co.tf([2, 0], [1, -1], 1) / co.tf([1, 0], [1], 1)  # div z
cl = co.feedback(ff, fb, sign=-1)
 
xin = [0, 1, 2, 1, 2, 3, 4, -2, -3]
_, y, _ = co.forced_response(cl, range(len(xin)), xin)
y.flatten()
 */
TEST(DelayFeedback, general_case){
    DoubleVector taps{0, 1, 1, 0, 0, 0};
    std::shared_ptr<SignalBlock> ff_ptr = std::make_shared<DiscreteTimeFilter>(taps);

    DoubleVector num {2}, den{1, -1};
    std::shared_ptr<SignalBlock> fb_ptr = std::make_shared<DiscreteTimeFilter>(num, den);
    InfInteraction::DelayFeedback fb(ff_ptr, fb_ptr, -1);
    DoubleVector xin{0, 1, 2, 1, 2, 3, 4, -2, -3};
    DoubleVector y_correct { 0.,   0.,   1.,   3.,   1.,  -7., -13.,   1.,  36.};

    for (int i=0; i < xin.size(); i++){
        DoubleVector y = fb.compute(DoubleVector {xin[i]});
        EXPECT_EQ(y_correct[i], y[0])<< "Fail at i=" << i;
    }
}


/*! Script to generate the desired result
 *
import control as co
ff = co.tf([0, 1, 2, -1, 0, 0], [1, 0, 0, 0, 0, 0], 1)
fb = co.tf([0], [1, -0.5], 1) / co.tf([1, 0], [1])  # div z
cl = co.feedback(ff, fb, sign=-1)

xin = [0, 1, 2, 1, 2, 3, 4, -2, -3]
_, y, _ = co.forced_response(cl, range(len(xin)), xin)
print(y)
 */
TEST(DelayFeedback, zero_feedback){
    DoubleVector taps{0, 1, 2, -1, 0, 0};
    std::shared_ptr<SignalBlock> ff_ptr = std::make_shared<DiscreteTimeFilter>(taps);

    DoubleVector num {0}, den{1, -0.5};
    std::shared_ptr<SignalBlock> fb_ptr = std::make_shared<DiscreteTimeFilter>(num, den);
    InfInteraction::DelayFeedback fb(ff_ptr, fb_ptr, -1);
    DoubleVector xin{0, 1, 2, 1, 2, 3, 4, -2, -3};
    DoubleVector y_correct {0., 0., 1., 4., 4., 2., 6., 8., 3.};

    for (int i=0; i < xin.size(); i++){
        DoubleVector y = fb.compute(DoubleVector {xin[i]});
        EXPECT_EQ(y_correct[i], y[0])<< "Fail at i=" << i;
    }
}

TEST(DiscreteTimeFilter_fir_scalar, case1){
    DoubleVector taps{0, 1, 2, -1, 0, 0};
    DiscreteTimeFilter filter (taps);
    DoubleVector xin{0, 1, 2, 1, 0};
    DoubleVector y_correct{0, 0, 1, 4, 4, 0, -1, 0, 0};

    // test
    for(int i=0; i < xin.size(); i++){
        double y = filter.compute(xin[i]);
        EXPECT_EQ(y_correct[i], y)<< "Fail at i=" << i;
    }
}

TEST(DiscreteTimeFilter_ccde_scalar, case1){
    DoubleVector num {0}, den{1, -0.99};
    DiscreteTimeFilter filter (num, den, 10);
    double y;
    // testing
    y = filter.compute(123.321);  // random number
    EXPECT_NEAR(y, 9.9, 1e-3);
    y = filter.compute(3.321);
    EXPECT_NEAR(y, 9.9 * 0.99, 1e-3);
}

TEST(DiscreteTimeFilter_ccde_scalar, case2){
    // implement the filter y[n] - y[n-1] = x[n]
    DoubleVector b{1}, a{1, -1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    y = filter.compute(1);
    EXPECT_EQ(y, 1);
    y = filter.compute(1.5);
    EXPECT_EQ(y, 2.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 0.0);
}


TEST(DiscreteTimeFilter_ccde_scalar, num_longer){
    DoubleVector b {0, 2, -2}, a{1, 3};
    DiscreteTimeFilter filter (b, a);
    DoubleVector xin{0, 1, 2, 3, 4, 3};
    DoubleVector y_correct{0.,   0.,   2.,  -4.,  14., -40.};

    // test
    for(int i=0; i < xin.size(); i++){
        double y = filter.compute(xin[i]);
        EXPECT_EQ(y_correct[i], y) << "Failed at i = " << i;
    }
}

TEST(DiscreteTimeFilter_ccde_scalar, order_zero){
    // implement the filter y[n] = -x[n]
    DoubleVector b{-1}, a{1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    y = filter.compute(1);
    EXPECT_EQ(y, -1);
    y = filter.compute(1.5);
    EXPECT_EQ(y, -1.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 2.5);
}

TEST(IIR_vec, case1){
    DoubleVector num {0}, den{1, -0.99};
    DiscreteTimeFilter filter (num, den, 10);
    DoubleVector y {0};
    DoubleVector u {0};
    // testing
    u = filter.compute(DoubleVector {132.321});  // random number
    EXPECT_NEAR(u[0], 9.9, 1e-3);
    u = filter.compute(DoubleVector {3.321});
    EXPECT_NEAR(u[0], 9.9 * 0.99, 1e-3);
}

TEST(IIR_vec, case2){
    // implement the filter y[n] - y[n-1] = x[n]
    DoubleVector b{1}, a{1, -1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    DoubleVector yvec;
    y = filter.compute(1);
    EXPECT_EQ(y, 1);
    // alternate computation between scalar and vector interface.
    yvec = filter.compute(DoubleVector {1.5});
    EXPECT_EQ(yvec[0], 2.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 0.0);
}


TEST(IIR_vec, order_zero){
    // implement the filter y[n] = -x[n]
    DoubleVector b{-1}, a{1};
    DiscreteTimeFilter filter (b, a, 0);
    DoubleVector y;
    y = filter.compute(DoubleVector {1});
    EXPECT_EQ(y[0], -1);
    y = filter.compute(DoubleVector{1.5});
    EXPECT_EQ(y[0], -1.5);
    y = filter.compute(DoubleVector{-2.5});
    EXPECT_EQ(y[0], 2.5);
}
