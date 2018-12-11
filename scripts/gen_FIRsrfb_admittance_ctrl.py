"""This script generates decoupled classical Admittance controller as
two FIR banks L and MB2:

y (measured force)             u (output)
---> [L] ---o------------------->
           -|              |
            \----[MB2]-----/
"""
import control as co
import sys
import numpy as np
import yaml

def main():
    m, b, k = map(float, sys.argv[1:4])
    nu, ny = map(int, sys.argv[4:])
    # m = 1.0
    # b = 3.0
    # k = 5.0
    # ny = 3
    # nu = 3

    tf_c = co.tf([1], [m, b, k])
    tf_d = co.c2d(tf_c, 0.008)
    num = tf_d.num[0][0].tolist()
    den = tf_d.den[0][0].tolist()
    # balance num and den
    num = [0] * (len(den) - len(num)) + num

    T = len(den)
    L = np.zeros((T, nu, ny))
    MB2 = np.zeros((T, nu, nu))

    for i in range(T):
        L[i, :] = np.eye(nu) * num[i]
        if i != 0:
            MB2[i, :] = np.eye(nu) * den[i]

    admittance3D = {
        'L': L.flatten().tolist(),
        'MB2': MB2.flatten().tolist(),
        'T': T,
        'nu': nu,
        'ny': ny
    }

    string = yaml.dump({'admittance3D': admittance3D})
    print(string)


if __name__ == '__main__':
    main()
