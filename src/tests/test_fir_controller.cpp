//
// Created by hung on 29/10/18.
//

#include "gtest/gtest.h"
#include <infinite_interaction_lib.h>


// impulse L is zero, expect zero output
TEST(FIR, zeroL_scalar){
    dVector L {0, 0, 0, 0}, MB2{0, -2, 3, 5};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector yin{0, 2, 3, 9, 2, 3};
    dVector u_n;
    for (int i = 0; i < yin.size(); ++i) {
        dVector y_n{yin[i]};
        u_n = controller.compute(y_n);
        EXPECT_EQ(u_n[0], 0);
    }
}

// impulse L is zero, expect zero output
TEST(FIR, zeroL_multi){
    dVector L {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, MB2{0, 0, 0, 0, -2, 3, 5, 2, 3, 4, 5, 2};

    FIRsrfb controller (3, 2, 2, L, MB2);  // 4-steps horizon
    dVector yin{0, 2, 3, 9, 2, 3}, yin2{2, 3, 9, 0, 2, 3}, u_n;
    for (int i = 0; i < yin.size(); ++i) {
        dVector y_n{yin[i], yin2[i]};
        u_n = controller.compute(y_n);
        EXPECT_EQ(u_n[0], 0);
        EXPECT_EQ(u_n[1], 0);
    }
}

// impulse MB2 is zero, expect feedthrough by L
TEST(FIR, zeroMB2_scalar){
    // if the L impulse is zero, expect zero output
    dVector L {2, 0, 0, 0}, MB2{0, 0, 0, 0};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector yin{0, 2, 3, 9, 2, 3};
    dVector u_n;
    for (int i = 0; i < yin.size(); ++i) {
        dVector y_n{yin[i]};
        u_n = controller.compute(y_n);
        EXPECT_EQ(u_n[0], 2 * yin[i]);
    }
}


// impulse MB2 is zero, expect feedthrough by L
TEST(FIR, zeroMB2_scalar_case2){
    // if the L impulse is zero, expect zero output
    dVector L {1, 1, 0, 0}, MB2{0, 0, 0, 0};

    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector yin{0, 2, 3, 9, 2, 3}, udesired{0, 2, 5, 12, 11, 5};
    dVector u_n;
    for (int i = 0; i < yin.size(); ++i) {
        dVector y_n{yin[i]};
        u_n = controller.compute(y_n);
        EXPECT_EQ(udesired[i], u_n[0]);
    }
}


// impulse MB2 is zero, expect
TEST(FIR, zeroMB2_multi){

}


// general case, computed using the below python script

//# python script: generate response
//import SLSsyn as Ss
//import numpy as np
//import control as co
//L = np.array([0, 1, 2, 3], dtype=float).reshape(-1, 1, 1)
//MB2 = np.array([0, -2, 3, 5], dtype=float).reshape(-1, 1, 1)
//Lz = Ss.impulse2tf(L, 1)
//MB2z = Ss.impulse2tf(MB2, 1)
//CLz = Lz * co.feedback(1, MB2z, sign=-1)
//y = [0, 2, 3, 9, 2, 3, 10, 20, 100]
//t, u, x = co.forced_response(CLz, np.arange(len(y)), y)
//print(u)  # got this:  [[   0.    0.    2.   11.   37.   60.  -12. -367. -949.]]
TEST(FIR, general_scalar){
    dVector L {0, 1, 2, 3}, MB2{0, -2, 3, 5};
    FIRsrfb controller (4, 1, 1, L, MB2);  // 4-steps horizon
    dVector y{0, 2, 3, 9, 2, 3, 10, 20, 100}, udesired{0, 0, 2, 11, 37, 60, -12, -367, -949};
    dVector ui;
    for (int i = 0; i < y.size(); ++i) {
        ui = controller.compute(dVector {y[i]});
        EXPECT_EQ(udesired[i], ui[0]) << "Error at i=" << i;
    }
}


// general case with multi-dimensional input and output
// the desired output (udesired) is generated using the below script

// # python script: generate response
// import SLSsyn as Ss
// import numpy as np
// import control as co
// np.random.seed(0)
// # T = 4, nu = 3, ny = 2
// L = np.random.randint(0, 10, size=(4, 3, 2))
// MB2 = np.random.randint(0, 10, size=(4, 3, 3))
// MB2[0, :] = 0
// Lz = Ss.mtf2ss(Ss.impulse2tf(L, 1), minreal=True)
// MB2z = Ss.mtf2ss(Ss.impulse2tf(MB2, 1), minreal=True)
// I3 = co.ss([], [], [], np.eye(3), 1)
// CLz = co.feedback(I3, MB2z, sign=-1) * Lz
// yin = np.random.randint(0, 10, size=(5, 2))
// t, uout, x = co.forced_response(CLz, np.arange(yin.shape[0]), yin, transpose=True)
//
// # reproduce results
// alpha0 = np.dot(L[0], yin[0])
// alpha1 = np.dot(L[0], yin[1]) + np.dot(L[1], yin[0])
// alpha2 = np.dot(L[0], yin[2]) + np.dot(L[1], yin[1]) + np.dot(L[0], yin[2])
// beta0 = np.zeros(3)
// u0 = alpha0 - beta0
// beta1 = MB2[1].dot(u0)
// u1 = alpha1 - beta1
// beta2 = MB2[1].dot(u1) + MB2[2].dot(u0)
// u2 = alpha2 - beta2

TEST(FIR, general_multi){
    dVector L {5, 0, 3, 3, 7, 9, 3, 5, 2, 4, 7, 6, 8, 8, 1, 6, 7, 7, 8, 1, 5, 9, 8, 9},
            MB2 {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 3, 7, 0, 1, 9, 9, 0, 4, 7, 3,
                 2, 7, 2, 0, 0, 4, 5, 5, 6, 8, 4, 1, 4, 9};
    std::vector<std::vector<double > > yin = {{8, 1}, {1, 7}, {9, 9}, {3, 6}, {7, 2}},
    udesired = {{ 4.000000e+01,  2.700000e+01,  6.500000e+01},
                {-2.820000e+02, -2.650000e+02, -7.360000e+02},
                { 2.877000e+03,  2.170000e+03,  9.467000e+03},
                {-3.198000e+04, -1.783600e+04, -1.074580e+05},
                { 3.392790e+05,  1.484960e+05,  1.162159e+06}};
    FIRsrfb controller (4, 2, 3, L, MB2);
    // (T, nu, ny) = (4, 3, 2)
    for (int i=0; i < yin.size(); ++i){
        dVector ui = controller.compute(yin[i]);
        for(int j=0; j < 3; ++j){
            EXPECT_DOUBLE_EQ(udesired[i][j], ui[j]) << "Fail at i=" << i << ", j=" << j;
        }
    }
}


TEST(FIR, badinput_MB2_0_nonzero){
    dVector L {0, 1, 2, 3}, MB2{1, -2};
     EXPECT_THROW(
            {
                try {
                    FIRsrfb controller (2, 2, 1, L, MB2);  // 4-steps horizon
                }
                catch (const std::invalid_argument & e){
                    EXPECT_STREQ("MB2[0] needs to be zero", e.what());
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



// throw exception if the input vectors have wrong shape
TEST(FIR, badinput_wrongshape_L){
    dVector L (2 * 3 * 2 - 1), MB2 (2 * 2 * 2);
    int T = 2, ny = 3, nu = 2;
    EXPECT_THROW(
            {
                try {
                    FIRsrfb controller (T, ny, nu, L, MB2);  // 4-steps horizon
                }
                catch (const std::invalid_argument & e){
                    EXPECT_STREQ("Wrong input shape", e.what());
                    throw;  // simply rethrow the catched expression
                }
            }, std::invalid_argument);
}



