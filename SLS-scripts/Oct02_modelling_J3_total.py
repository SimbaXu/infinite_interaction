import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# In this script, the multiplicative uncertainty wI of the loop gain
# transfer function L:
# Lp = L (1 + wI Delta)
# is identified.

s = np.complex(0, 1) * np.logspace(-3, 4, 50)

R1nom = (-s + 55.56) / (s + 55.56) / (0.0437 * s + 1)
Larmnom = 0.4
Hnom = (50 + 5 * s)
Senom = 2 * np.pi * 73 / (s + 2 * np.pi * 73)
R2nom = 0.475 * s ** 2 / (s / 100 + 1)**2 * 35 / (s + 35) * (-s + 66) / (s + 66)
Lnom = (R1nom * Larmnom * Hnom + R2nom) * Senom * Larmnom

# perb
perbs = []
Lperbs = []
for Larm, kH, bH, KR2, tauR2 in product(np.linspace(0.35, 0.45, 2),
                                        np.linspace(0, 80, 3),
                                        np.linspace(0, 10, 3),
                                        np.linspace(0.455, 0.495, 3),
                                        np.linspace(0.025, 0.035)):
    R1perb = np.exp(-0.036 * s) / (0.0437 * s + 1)
    Larmperb = float(Larm)
    Hperb = (kH + bH * s)
    Seperb = np.copy(Senom)
    R2perb = KR2 * s ** 2 / (s / 100 + 1)**2 * 35 / (s + 35) * np.exp(-tauR2 * s)
    Lperb = (R1perb * Larmperb * Hperb + R2perb) * Seperb * Larmperb
    lI = (Lperb - Lnom) / Lnom
    perbs.append(lI)
    Lperbs.append(Lperb)

# fit
fig = plt.figure()
ax = plt.subplot(1, 1, 1)
for Lperb in Lperbs:
    ax.plot(np.abs(s), np.abs(Lperb), '--', c='gray')
ax.plot(np.abs(s), np.abs(Lnom), 'x-', label="nom")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
plt.show()

# fit
# wI = (s * (4.7 / 300)**(1.0 / 3) + np.power(1.025, 1.0 / 3))**3 / (s / np.power(300, 1.0/3) + 1)**3
wI = (s * 0.25 + 1.0083)**3 / (s * 0.149 + 1)**3;
fig = plt.figure()
ax = plt.subplot(1, 1, 1)
for perb in perbs:
    ax.plot(np.abs(s), np.abs(perb), '--', c='gray')
ax.plot(np.abs(s), np.abs(wI), '-', c='blue')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
plt.show()


# fit
fig = plt.figure()
ax = plt.subplot(2, 1, 1)
ax.plot(np.abs(s), np.abs(Lnom), 'x-', label="nom")
for Lperb in Lperbs:
    ax.plot(np.abs(s), np.abs(Lperb), '--')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()

ax = plt.subplot(2, 1, 2)
ax.plot(np.abs(s), np.angle(Lnom), 'x-', label="model")
ax.legend()
ax.set_xscale('log')
ax.grid()
plt.show()

