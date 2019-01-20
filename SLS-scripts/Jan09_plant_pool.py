"""This source file contains LTI plant models.

Compatible with python 3.6. See associate org-node at 'system design
notes' for more details. See Mujin, page 10.

I use sympy to derive analytical expressions of the model. The whole
process has two steps. First, one draws block/model and identifies
symbolic expression. Second, one uses sympy to derive the expressions
of interest.

Each class is a model.
"""
import sympy as sym
import control as co
import numpy as np
import SLSsyn as Ss


class PlantV1:
    """ Point robot with attached spring.

    The robot has flexible link. At the tip of the link there is a FT sensor. The tool attached to
    the FT sensor is also flexible.

    Signal descriptions:
        inputs = [f_desired, x_env, f_noise, f_robust, w_add, u]
        outputs = [f_measure, x_env_diff, x_robot, z_add, y1, y2]
    """
    gen_exp_file = "PlantV1_gen_exp"
    z_idx_map = {
        0: 'fm', 1: 'x', 2: 'xro'
    }
    w_idx_map = {
        0: 'fd', 1: 'xe', 2: 'n', 3: 'fro'
    }

    input_descriptions = ['f_desired', 'x_env', 'f_noise', 'f_robust', 'w_add', 'u']
    # output_descriptions = ['f_measure', 'x_env_diff', 'x_robot', 'z_add', 'y1', 'y2']
    output_descriptions = ['f_measure', 'x_tip', 'x_link', 'x_robot']
    @staticmethod
    def plant(tau_R1=0.0437, K_env=10e3, omega_add=1,
              b_link=200, k_link=100e3, m_link=13,
              k_tip=100e4, b_tip=5000, m_tip=0.1, Ts=0.008):
        """Apply dynamic parameters.

        k_link, b_link and m_link are tuned so that the step response
        of the measure force to robot postion command is almost identical to
        the trajectories recorded in this bag: `2018-01-18-e1.bag`.
        """

        # constant transfer function
        z = co.tf([1, 0], [1], Ts)
        s = co.tf([1, 0], [1])
        R1 = 1 / (1 + tau_R1 * s)  # robot first-order response

        with open(PlantV1.gen_exp_file, 'r') as f:
            expr_string = f.read()
        exec(expr_string)
        P = Ss.tf_blocks(locals()["Plist"])
        return P

    @staticmethod
    def derive():
        """ Derive the symbolic expression of the plant, then print code.

        See Muji, p.10, Model force-control
        """
        # dynamic symbols
        symbols_dynamic_string = 'v_tip, x_tip, v_link, x_link, v_robot, x_robot, x_cmd, f_measure, f_noise, f_act, f_desired, w_add, z_add, x_env, f_env, x_env_diff, u, y1, y2, f_robust'
        v_tip, x_tip, v_link, x_link, v_robot, x_robot, x_cmd, f_measure, f_noise, f_act, f_desired, w_add, z_add, x_env, f_env, x_env_diff, u, y1, y2, f_robust = sym.symbols(symbols_dynamic_string)
        symbols_dynamic_all = sym.symbols(symbols_dynamic_string)

        # physical symbols
        K_env, omega_add, m_link, k_link, b_link, R1 = sym.symbols("K_env, omega_add, m_link, k_link, b_link, R1")
        m_tip, k_tip, b_tip = sym.symbols("m_tip, k_tip, b_tip ")

        # misc symbols
        s, Ts = sym.symbols('s Ts')

        # symbols_input = [f_desired, x_env, f_noise, f_robust, w_add, u]
        symbols_input = [sym.symbols(_x) for _x in PlantV1.input_descriptions]
        symbols_output = [sym.symbols(_x) for _x in PlantV1.output_descriptions]
        # symbols_output = [f_measure, x_env_diff, x_robot, z_add, y1, y2]

        # dynamic equation
        Eq = sym.Eq
        eqs = []
        eqs.extend([
            # robot dynamics
            Eq(x_robot, R1 * x_cmd * sym.exp(- 2 * Ts)),
            Eq(x_cmd, u),
            Eq(v_robot, s * x_robot),

            # tip link dynamics
            Eq(v_link, s * x_link),
            Eq(v_tip, s * x_tip),
            Eq(m_tip * x_tip * s ** 2, f_act + k_tip * (x_link - x_tip) + b_tip * (v_link - v_tip)),
            Eq(m_link * x_link * s ** 2, k_tip * (x_tip - x_link) + b_tip * (v_tip - v_link)
               + k_link * (x_robot - x_link) + b_link * (v_robot - v_link)),

            # force measurement
            Eq(f_measure, k_tip * (x_tip - x_link) + b_tip * (v_tip - v_link)),
            Eq(y2, f_measure + f_noise),

            # others
            Eq(f_desired, y1),

            # env dynamics
            Eq(x_env_diff, w_add + x_link - x_env),
            Eq(z_add, omega_add * (x_link - x_env)),
            Eq(f_env, K_env * x_env_diff),
            Eq(f_act, - f_env - f_robust),
        ]
        )

        symbols_to_solve = []
        for symbol_out in symbols_dynamic_all:
            if symbol_out not in symbols_input:
                symbols_to_solve.append(symbol_out)

        assert len(symbols_to_solve) == len(eqs), "Error: The number of equalities and the number of symbols to solve are different."

        print("Start solving for dynamic equations")
        solve_result = sym.solve(eqs, symbols_to_solve)
        # form overall tranfer matrix
        plant = []
        for symbol_out in symbols_output:
            transfer_function = sym.expand(solve_result[symbol_out])
            terms_dictionary = sym.collect(transfer_function, symbols_input, evaluate=False)
            # This error occurs if the relevant transfer function is not a sum of
            # terms that are linear in all input symbols.
            assert 1 not in terms_dictionary, "A constant appears. Error occurs."
            # extract individual component transfer function
            plant.append([])
            for symbol_in in symbols_input:
                if symbol_in in terms_dictionary:
                    plant[-1].append(terms_dictionary[symbol_in].simplify())
                else:
                    plant[-1].append(0)
        plant = sym.Matrix(plant)

        # print symbolic transfer matrix as string.
        plant_string_expr = print_symmatrix_as_list(plant)
        with open(PlantV1.gen_exp_file, 'w') as f:
            f.write(plant_string_expr)
        print("Obtain the following string expression for this plant\n----------------------------------------\n\n")
        print(plant_string_expr)


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
                    output += "co.c2d({:}, Ts)\n".format(P[i, j])
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
                    output += "z ** ({:d}) * co.c2d({:}, Ts)\n".format(power, factored_args[exp_part])

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
    plant = PlantV1.plant(K_env=1, b_link=200, k_link=100e3, m_link=13, k_tip=100e4, b_tip=5000)
    Ss.plot_step_responses(plant, 1, PlantV1.input_descriptions, PlantV1.output_descriptions)
    import IPython
    if IPython.get_ipython() is None:
        IPython.embed()

