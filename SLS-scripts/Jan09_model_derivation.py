# python 3. See associate org-node at 'system design notes' for more
# details. See Mujin, page 10.
import sympy as sym

# dynamic symbols
x, v, c, xe, ke, xr, vr = sym.symbols('x, v, c, xe, ke, xr, vr')
# tool position, velocity and command
# environment position, stiffness
# robot position and velocity
kr, br, mr, s = sym.symbols('kr, br, mr, s')  # tool stiffness, damping, mass and Laplace' s
fm, fd, n = sym.symbols('fm fd n')  # measured force, desired force, noise
f, R1 = sym.symbols('f R1')  # acting force, p2p transfer function
Ts = sym.symbols('Ts') # sampling time

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

    # loop
    ke * (x - xe) + f,

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

symdyn = [vr, xr, v, x, fm, f, c, fd, xe, n]
symannot_out = [z1, xr, y1, y2]
symannot_in = [w1, w2, w3, u1]

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
