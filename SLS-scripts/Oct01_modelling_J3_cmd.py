import numpy as np
import matplotlib.pyplot as plt

# data: q3 cmd -> tau
d_omega = [1, 4, 8, 12, 16, 20, 24]
d_mag = [1 ,   0.98 ,   0.94 ,   0.88 ,   0.81 ,  0.75 ,   0.69]
d_angle = [-0.08 , -0.292 , -0.630 ,  -0.85 , -1.183 , -1.44 , -1.665]

# model
s = np.complex(0, 1) * np.linspace(1, 24, 100)
freq_resp = np.exp(- 0.036 * s) / (0.0437 * s + 1)

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

