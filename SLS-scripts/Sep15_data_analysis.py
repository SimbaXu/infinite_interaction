import matplotlib.pyplot as plt
import numpy as np
import re

file_name = "data/admittance_control.log"
# file_name = "data/admittance_control_loopshape_try2_stable.log"

# for removing ANSI escape seq
ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')

with open(file_name) as f:
    doc = f.read()
lines = doc.split('\n')
lines = filter(lambda l: 'ros.infinite_interaction' not in l, lines)
lines = filter(lambda l: 'tau' in l or 'cmd' in l, lines)
lines = map(lambda l: ansi_escape.sub("", l), lines)
tau_data = map(lambda l: l.split(':')[-1], filter(lambda l0: 'tau' in l0, lines))
cmd_data = map(lambda l: l.split(':')[-1], filter(lambda l0: 'cmd' in l0, lines))

tau_idx = []
tau = []
for l in tau_data:
    tau_idx.append(int(l.split(',')[0]))
    tau.append(np.array(map(float, l.split(',')[1:])))

cmd_idx = []
cmd = []
for l in cmd_data:
    cmd_idx.append(int(l.split(',')[0]))
    cmd.append(np.array(map(float, l.split(',')[1:])))

tau = np.array(tau)
tau_idx = np.array(tau_idx)
cmd = np.array(cmd)
cmd_idx = np.array(cmd_idx)

fig = plt.figure()
plt.plot(tau_idx, tau[:, 2], label="tau[2]")
plt.plot(cmd_idx, cmd[:, 2] - 1, label="cmd[2]")
plt.legend()
plt.show()
