import numpy as np
import matplotlib.pyplot as plt

# data: q3 cmd -> tau
d_omega = [1, 4, 8, 12, 16, 20, 24]
d_mag = [0.2, 3.42 ,   13.9 ,     30 ,   49.8 ,  72.9 ,   98.1]
d_angle = [0 , -0.297 , -0.586 , -0.826 , -1.118 , -1.27 ,  -1.50 ]

# model
s = np.complex(0, 1) * np.logspace(-1, 5, 1000)
K = 100
# freq_resp = - 0.19 * s**2 * np.exp(- 0.036 * s) / (0.0237 * s + 1) / ((s / K) ** 2 + 2 * s/K + 1**2)

freq_resp = - (0.19
               * s / (s / K + 1)
               * s / (s / K + 1)
               * 35 / (s + 35)
               * np.exp(- 0.03 * s)
)

# fit
fig = plt.figure()
ax = plt.subplot(2, 1, 1)
ax.plot(np.abs(s), np.abs(freq_resp), 'x-', label="model")
ax.plot(d_omega, d_mag, 'o--', label='data')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()

ax = plt.subplot(2, 1, 2)
ax.plot(np.abs(s), np.angle(freq_resp), 'x-', label="model")
ax.plot(d_omega, d_angle, 'o--', label='data')
ax.legend()
ax.set_xscale('log')
ax.grid()
plt.show()

