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
    dVector taps{0, 1, 1, 0, 0, 0};
    std::shared_ptr<SignalBlock> ff_ptr = std::make_shared<DiscreteTimeFilter>(taps);

    dVector num {2}, den{1, -1};
    std::shared_ptr<SignalBlock> fb_ptr = std::make_shared<DiscreteTimeFilter>(num, den);
    InfInteraction::DelayFeedback fb(ff_ptr, fb_ptr, -1);
    dVector xin{0, 1, 2, 1, 2, 3, 4, -2, -3};
    dVector y_correct { 0.,   0.,   1.,   3.,   1.,  -7., -13.,   1.,  36.};

    for (int i=0; i < xin.size(); i++){
        dVector y = fb.compute(dVector {xin[i]});
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
    dVector taps{0, 1, 2, -1, 0, 0};
    std::shared_ptr<SignalBlock> ff_ptr = std::make_shared<DiscreteTimeFilter>(taps);

    dVector num {0}, den{1, -0.5};
    std::shared_ptr<SignalBlock> fb_ptr = std::make_shared<DiscreteTimeFilter>(num, den);
    InfInteraction::DelayFeedback fb(ff_ptr, fb_ptr, -1);
    dVector xin{0, 1, 2, 1, 2, 3, 4, -2, -3};
    dVector y_correct {0., 0., 1., 4., 4., 2., 6., 8., 3.};

    for (int i=0; i < xin.size(); i++){
        dVector y = fb.compute(dVector {xin[i]});
        EXPECT_EQ(y_correct[i], y[0])<< "Fail at i=" << i;
    }
}

TEST(DiscreteTimeFilter_fir_scalar, case1){
    dVector taps{0, 1, 2, -1, 0, 0};
    DiscreteTimeFilter filter (taps);
    dVector xin{0, 1, 2, 1, 0};
    dVector y_correct{0, 0, 1, 4, 4, 0, -1, 0, 0};

    // test
    for(int i=0; i < xin.size(); i++){
        double y = filter.compute(xin[i]);
        EXPECT_EQ(y_correct[i], y)<< "Fail at i=" << i;
    }
}

TEST(DiscreteTimeFilter_ccde_scalar, case1){
    dVector num {0}, den{1, -0.99};
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


TEST(DiscreteTimeFilter_ccde_scalar, num_longer){
    dVector b {0, 2, -2}, a{1, 3};
    DiscreteTimeFilter filter (b, a);
    dVector xin{0, 1, 2, 3, 4, 3};
    dVector y_correct{0.,   0.,   2.,  -4.,  14., -40.};

    // test
    for(int i=0; i < xin.size(); i++){
        double y = filter.compute(xin[i]);
        EXPECT_EQ(y_correct[i], y) << "Failed at i = " << i;
    }
}

TEST(DiscreteTimeFilter_ccde_scalar, order_zero){
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

TEST(IIR_vec, case1){
    dVector num {0}, den{1, -0.99};
    DiscreteTimeFilter filter (num, den, 10);
    dVector y {0};
    dVector u {0};
    // testing
    u = filter.compute(dVector {132.321});  // random number
    EXPECT_NEAR(u[0], 9.9, 1e-3);
    u = filter.compute(dVector {3.321});
    EXPECT_NEAR(u[0], 9.9 * 0.99, 1e-3);
}

TEST(IIR_vec, case2){
    // implement the filter y[n] - y[n-1] = x[n]
    dVector b{1}, a{1, -1};
    DiscreteTimeFilter filter (b, a, 0);
    double y;
    dVector yvec;
    y = filter.compute(1);
    EXPECT_EQ(y, 1);
    // alternate computation between scalar and vector interface.
    yvec = filter.compute(dVector {1.5});
    EXPECT_EQ(yvec[0], 2.5);
    y = filter.compute(-2.5);
    EXPECT_EQ(y, 0.0);
}


TEST(IIR_vec, order_zero){
    // implement the filter y[n] = -x[n]
    dVector b{-1}, a{1};
    DiscreteTimeFilter filter (b, a, 0);
    dVector y;
    y = filter.compute(dVector {1});
    EXPECT_EQ(y[0], -1);
    y = filter.compute(dVector{1.5});
    EXPECT_EQ(y[0], -1.5);
    y = filter.compute(dVector{-2.5});
    EXPECT_EQ(y[0], 2.5);
}
