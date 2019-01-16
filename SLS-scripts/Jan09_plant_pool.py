""" Contain several plants. A plant pool.

Must be ran with python 3.6

See associate org-node at 'system design notes' for more
details. See Mujin, page 10.
"""
import sympy as sym
import control as co
import numpy as np
import SLSsyn as Ss


class PlantV1:
    gen_exp_file = "PlantV1_gen_exp"
    z_idx_map = {
        0: 'fm', 1: 'x', 2: 'xro'
    }
    w_idx_map = {
        0: 'fd', 1: 'xe', 2: 'n', 3: 'fro'
    }

    @staticmethod
    def plant(tau_R1=0.0437, ke=10e3, br=10, kr=500, mr=0.05, Ts=0.008):
        """
        This model is generated by running Jan09_model_derivation.py
        outputs = [fm, x, y1, y2]
        inputs = [fd, xe, n, u1]
        """
        import SLSsyn as Ss
        # constant transfer function
        z = co.tf([1, 0], [1], Ts)
        s = co.tf([1, 0], [1])

        def dist(exp):
            if type(exp) == float:
                return exp
            return co.c2d(exp, Ts)
        R1 = 1 / (1 + tau_R1 * s)

        with open(PlantV1.gen_exp_file, 'r') as f:
            expr_string = f.read()
        exec(expr_string)
        P = Ss.tf_blocks(locals()["Plist"])
        return P

    @staticmethod
    def derive():
        """See Muji, p.10, Model force-control
        """
        # dynamic symbols
        x, v, c, xe, ke, xr, vr = sym.symbols('x, v, c, xe, ke, xr, vr')
        # tool position, velocity and command
        # environment position, stiffness
        # robot position and velocity
        kr, br, mr, s = sym.symbols('kr, br, mr, s')  # tool stiffness, damping, mass and Laplace' s
        fm, fd, n = sym.symbols('fm fd n')  # measured force, desired force, noise
        f, R1 = sym.symbols('f R1')  # acting force, p2p transfer function
        Ts = sym.symbols('Ts')  # sampling time
        fro, xro = sym.symbols('fro, xro')  # force and displacement for robust stability

        # annotated symbols
        y1, y2, z1, z2, z3, u1, w1, w2, w3 = sym.symbols('y1, y2, z1, z2, z3, u1, w1, w2, w3')

        eqs = []
        eqs.extend([
            # robot dynamics
            xr - R1 * c * sym.exp(- 2 * Ts),
            vr - s * xr,
            v - s * x,
            fm - (kr * (x - xr) + br * (v - vr)) * sym.exp(- 2 * Ts),
            - (kr * (x - xr) + br * (v - vr)) + f - mr * v * s,

            xro - x + xe,  # xro = x - xe
            # loop
            ke * (x - xe) + f + fro,

            # input/output categorization
            fd - y1,
            y2 - (fm + n),
            c - u1,
            fd - w1,
            xe - w2,
            n - w3,
            fm - z1,
        ]
        )

        symdyn = [vr, xr, v, x, fm, f, c, fd, xe, n, xro]
        symannot_out = [z1, xr, xro, y1, y2]
        symannot_in = [w1, w2, w3, fro, u1]

        sym2solve = symdyn + symannot_out
        res_ = sym.solve(eqs, sym2solve)
        # form overall tranfer matrix
        P = []
        for sym_out in symannot_out:
            t_ = sym.expand(res_[sym_out])
            terms = sym.collect(t_, symannot_in, evaluate=False)
            P.append([])
            for sym_in in symannot_in:
                if sym_in in terms:
                    P[-1].append(terms[sym_in].simplify())
                else:
                    P[-1].append(0)
                assert 1 not in terms, "there should not be any constant"
        P = sym.Matrix(P)
        P

        string_expr = print_symmatrix_as_list(P)
        with open(PlantV1.gen_exp_file, 'w') as f:
            f.write(string_expr)

        print(string_expr)


def contain_exp(expr):
    """ Return True if expression contains the exponential operator.
    """
    for arg in sym.preorder_traversal(expr):
        # check if there is an exponential function
        if arg.func == sym.exp:
            return True
    return False


def print_symmatrix_as_list(P):
    """Print code to evaluate a continous-time transfer matrix to a
    discrete time one.
    """
    Ts = sym.symbols('Ts')  # sampling time
    fail = False
    output = ""
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            output += "P{:d}{:d} = ".format(i, j)
            # best case, no exponential
            if not contain_exp(P[i, j]):
                if P[i, j].is_constant():
                    output += "{:}\n".format(P[i, j])
                else:
                    output += "dist({:})\n".format(P[i, j])
                continue

            # check if the exponential form can be factored out
            factored_args = P[i, j].collect(sym.exp(Ts), evaluate=False)
            if any([contain_exp(arg) for arg in factored_args.values()]):
                output += "fail!!\n"
                fail = True
            else:
                for exp_part in factored_args:
                    # compute the power raised by the exponential part
                    power = round(sym.log(exp_part, sym.exp(Ts)).evalf(subs={Ts: 1}))
                    output += "z ** ({:d}) * dist({:})\n".format(power, factored_args[exp_part])

    output += "Plist = [\n"
    for i in range(P.shape[0]):
        output += "["
        for j in range(P.shape[1]):
            output += "P{:d}{:d}, ".format(i, j)
        output += "],\n"
    output += "]"

    if fail:
        output = "CONVERSION FAILS. DO NOT USE\n" + output
    return output

class Controllers:
    """ Controller pool.
    """

    @staticmethod
    def PI_v1(Kp=1e-3, Ki=20e-3):
        """ integral controller

        Parameters selected for stiff environment with stifness = 10e3.
        """
        s = co.tf([1, 0], [1])
        K_ctime_single = Kp + Ki / s
        K_dtime_single = co.c2d(K_ctime_single, 0.008)
        PI_v1 = - Ss.tf_blocks([[K_dtime_single, - K_dtime_single]])
        return PI_v1

    @staticmethod
    def Qsyn(filename="Jan09_K_ss_params.npz"):
        data = np.load(filename)
        K = co.StateSpace(data['A'], data['B'], data['C'], data['D'], data['dt'])
        return K


if __name__ == "__main__":
    PlantV1.derive()
    # PlantV1.plant()


