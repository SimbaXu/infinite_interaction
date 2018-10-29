//
// Created by hung on 29/10/18.
//

#include "gtest/gtest.h"
#include <infinite_interaction_lib.h>


// impulse L is zero, expect zero output
TEST(FIR, zeroL_scalar){
    dVector L {0, 0, 0, 0}, MB2{0, -2, 3, 5};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector u{0, 2, 3, 9, 2, 3};
    double yi;
    for (int i = 0; i < u.size(); ++i) {
        yi = controller.compute(u[i]);
        EXPECT_EQ(yi, 0);
    }
}

// impulse L is zero, expect zero output
TEST(FIR, zeroL_multi){
    dVector L {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, MB2{0, 0, 0, 0, -2, 3, 5, 2, 3, 4, 5, 2};

    FIRsrfb controller (3, 2, 2, L, MB2);  // 4-steps horizon
    dVector u1{0, 2, 3, 9, 2, 3}, u2{2, 3, 9, 0, 2, 3}, yi;
    for (int i = 0; i < u1.size(); ++i) {
        dVector u{u1[i], u2[i]};
        yi = controller.compute(u);
        EXPECT_EQ(yi[0], 0);
        EXPECT_EQ(yi[1], 0);
    }
}

// impulse MB2 is zero, expect feedthrough by L
TEST(FIR, zeroMB2_scalar){
    // if the L impulse is zero, expect zero output
    dVector L {2, 0, 0, 0}, MB2{0, 0, 0, 0};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector u{0, 2, 3, 9, 2, 3};
    double yi;
    for (int i = 0; i < u.size(); ++i) {
        yi = controller.compute(u[i]);
        EXPECT_EQ(yi, 0);
    }
}


// impulse MB2 is zero, expect
TEST(FIR, zeroMB2_multi){
}


// general case, computed using the below python script

// # python script: generate response
// import SLSsyn as Ss
// import numpy as np
// import control as co
// L = np.array([0, 1, 2, 3], dtype=float).reshape(-1, 1, 1)
// MB2 = np.array([0, -2, 3, 5], dtype=float).reshape(-1, 1, 1)
// Lz = Ss.impulse2tf(L, 1)
// MB2z = Ss.impulse2tf(MB2, 1)
// CLz = Lz * co.feedback(1, MB2z, sign=-1)
// u = [0, 2, 3, 9, 2, 3]
// t, y, x = co.forced_response(CLz, np.arange(len(u)), u)
// print(y)  # got this: [[ 0.  0.  2. 11. 37. 60.]]
TEST(FIR, general_scalar){
    dVector L {0, 1, 2, 3}, MB2{0, -2, 3, 5};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector u{0, 2, 3, 9, 2, 3}, ydesired{0, 0, 2, 11, 37, 60};
    double yi;
    for (int i = 0; i < u.size(); ++i) {
        yi = controller.compute(u[i]);
        EXPECT_EQ(yi, ydesired[i]);
    }
}


TEST(FIR, general_multi){}


TEST(FIR, badinput_MB2_0_nonzero){
    dVector L {0, 1, 2, 3}, MB2{0, -2, 3, 5};
     EXPECT_THROW(
            {
                try {
                    FIRsrfb controller (2, 2, 1, L, MB2);  // 4-steps horizon
                }
                catch (const std::invalid_argument & e){
                    EXPECT_STREQ("MB2[0] needs to be zero.", e.what());
                    throw;  // simply rethrow the catched expression
                }
            }, std::invalid_argument);
}


// throw exception if the input vectors have wrong shape
TEST(FIR, badinput_wrongshape){
    dVector L {0, 1, 2, 3}, MB2{0, -2, 3, 5};
    EXPECT_THROW(
            {
                try {
                    FIRsrfb controller (5, 1, 1, L, MB2);  // 4-steps horizon
                }
                catch (const std::invalid_argument & e){
                    EXPECT_STREQ("Wrong input shape", e.what());
                    throw;  // simply rethrow the catched expression
                }
            }, std::invalid_argument);
}





