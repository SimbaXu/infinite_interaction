import sys
sys.path.insert(0, "/home/hung/catkin_ws/src/infinite_interaction/SLS-scripts")
import SLSsyn as Ss
import pytest
import control as co
import numpy as np
import matplotlib.pyplot as plt
from conftest import random_dtime_sys


@pytest.mark.parametrize("i", [0])  # output index
@pytest.mark.parametrize("j", [0])  # input index
def test_Ny1_Nu1(i, j, Ts=0.008):
    """ Correct impulse response for y in R and u in R.
    """
    P = random_dtime_sys(5, 5, Ts)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(P, 1, 1)
    weight = np.array([0.1, 0.4, 0, 0.2])
    Tsim = np.arange(100) * Ts
    z = co.tf([1, 0], [1], Ts)
    Q = Ts * (0.1 + 0.4 * z ** (-1) + 0 + 0.2 * z ** (-3))

    # expected result
    H_expected = Pzw[i, j] + Pzu[i, 0] * Q * Pyw[0, j]
    _, imp_expected = co.impulse_response(H_expected, Tsim)
    imp_expected = imp_expected[0, :]

    # actual result
    imp_actual = Ss.Qsyn.obtain_time_response_var(weight, (i, j), Tsim, Pzw, Pzu, Pyw)

    # # debug plot
    # plt.plot(imp_actual)
    # plt.plot(imp_expected, '--')
    # plt.show()

    np.testing.assert_allclose(imp_expected, imp_actual, atol=1e-5)


@pytest.mark.parametrize("i", [0])  # output index
@pytest.mark.parametrize("j", [0])  # input index
@pytest.mark.parametrize("zero_Pzw", [True, False])  # input index
def test_Ny2_Nu1(i, j, zero_Pzw, Ts=0.008):
    """ Correct impulse response for y in R2 and u in R.
    """
    P = random_dtime_sys(5, 5, Ts)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(P, 1, 2)
    if zero_Pzw:
        Pzw = Ss.tf_blocks(np.zeros((3, 4)), dt=Ts)
    weight = np.random.randn(6)
    # weight = np.array([0, 0, 0., 1., 0.0, 0.0])  # for debug
    Tsim = np.arange(30) * Ts
    z = co.tf([1, 0], [1], Ts)
    Q00 = Ts * (weight[0] + weight[1] * z ** (-1) + weight[2] * z ** (-2))
    Q01 = Ts * (weight[3] + weight[4] * z ** (-1) + weight[5] * z ** (-2))

    # expected result
    H_expected = Pzw[i, j] + (Pzu[i, 0] * Q00 * Pyw[0, j]
                              + Pzu[i, 0] * Q01 * Pyw[1, j])
    _, imp_expected = co.impulse_response(H_expected, Tsim)
    imp_expected = imp_expected[0, :]

    # actual result
    imp_actual = Ss.Qsyn.obtain_time_response_var(weight, (i, j), Tsim, Pzw, Pzu, Pyw)

    # # debug plot
    # plt.plot(imp_actual)
    # plt.plot(imp_expected, '--')
    # plt.show()

    np.testing.assert_allclose(imp_expected, imp_actual)


@pytest.mark.parametrize("i", [0, 2])  # output index
@pytest.mark.parametrize("j", [0, 1])  # input index
def test_Ny2_Nu2(i, j, zero_Pzw=True, Ts=0.008):
    """ Correct impulse response for y in R and u in R.
    """
    P = random_dtime_sys(5, 5, Ts)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(P, 2, 2)
    if zero_Pzw:
        Pzw = Ss.tf_blocks(np.zeros((3, 3)), dt=Ts)
    weight = np.random.randn(12)
    # weight = np.array([0, 0, 0., 1., 0.0, 0.0])  # for debug
    Tsim = np.arange(30) * Ts
    z = co.tf([1, 0], [1], Ts)
    Q00 = Ts * (weight[0] + weight[1] * z ** (-1) + weight[2] * z ** (-2))
    Q01 = Ts * (weight[3] + weight[4] * z ** (-1) + weight[5] * z ** (-2))
    Q10 = Ts * (weight[6] + weight[7] * z ** (-1) + weight[8] * z ** (-2))
    Q11 = Ts * (weight[9] + weight[10] * z ** (-1) + weight[11] * z ** (-2))

    # expected result
    H_expected = Pzw[i, j] + (Pzu[i, 0] * Q00 * Pyw[0, j]
                              + Pzu[i, 0] * Q01 * Pyw[1, j]
                              + Pzu[i, 1] * Q10 * Pyw[0, j]
                              + Pzu[i, 1] * Q11 * Pyw[1, j])
    _, imp_expected = co.impulse_response(H_expected, Tsim)
    imp_expected = imp_expected[0, :]

    # actual result
    imp_actual = Ss.Qsyn.obtain_time_response_var(weight, (i, j), Tsim, Pzw, Pzu, Pyw)

    # # debug plot
    # plt.plot(imp_actual)
    # plt.plot(imp_expected, '--')
    # plt.show()

    np.testing.assert_allclose(imp_expected, imp_actual, atol=1e-8)


def test_wrong_weight_dim():
    """ Given weight have wrong dimension.
    """
    P = random_dtime_sys(5, 5, 0.008)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(P, 2, 3)

    weight = np.random.randn(11)  # 11 = 6 + 5
    Tsim = np.arange(100) * 0.008

    with pytest.raises(ValueError) as exinfo:
        Ss.Qsyn.obtain_time_response_var(weight, (0, 0), Tsim, Pzw, Pzu, Pyw)
    assert str(exinfo.value) == "Inconsistent weight dimension"


def test_wrong_Tsim():
    """ Given sim Time has wrong increment.
    """
    P = random_dtime_sys(2, 2, 0.008)
    Pzw, Pzu, Pyw, Pyu = Ss.get_partitioned_transfer_matrices(P, 1, 1)

    with pytest.raises(ValueError) as exinfo:
        Ss.Qsyn.obtain_time_response_var([0, 1, 2, 3], (0, 0), np.arange(5), Pzw, Pzu, Pyw)
    assert str(exinfo.value) == "Inconsistent simulation time `Tsim`"
