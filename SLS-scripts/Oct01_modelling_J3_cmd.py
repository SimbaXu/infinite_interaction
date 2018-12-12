import numpy as np
import matplotlib.pyplot as plt
import control as co

# data: q3 cmd -> tau
d_omega = [1, 4, 8, 12, 16, 20, 24]
d_mag = [1, 0.98, 0.94, 0.88, 0.81, 0.75, 0.69]
d_angle = [-0.08, -0.292, -0.630, -0.85, -1.183, -1.44, -1.665]

# model
omegas = np.linspace(1, 24, 100)

# # continuous model:
# # e^(-0.036 s) / (1 + 0.0437 s)
# freq_resp = np.exp(- 0.036 * 1j * omegas) / (0.0437 * 1j * omegas + 1)

# discrete-time model:
# z^-4 * c2d(1 / (1 + 0.0437 s))
z = co.tf([1, 0], [1], 0.008)
model = co.c2d(co.tf([1], [0.0437, 1]), 0.008) * z**(-4)
mag, phase, _ = model.freqresp(omegas)
freq_resp = mag[0, 0] * np.exp(1j * phase[0, 0])
print(model)

# fit
fig = plt.figure()
ax = plt.subplot(2, 1, 1)
ax.plot(omegas, np.abs(freq_resp), 'x-', label="model")
ax.plot(d_omega, d_mag, 'o--', label='data')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()

ax = plt.subplot(2, 1, 2)
ax.plot(omegas, np.angle(freq_resp), 'x-', label="model")
ax.plot(d_omega, d_angle, 'o--', label='data')
ax.legend()
ax.set_xscale('log')
ax.grid()
plt.show()
