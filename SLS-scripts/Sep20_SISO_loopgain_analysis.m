figure(1)
clf()
L1 = RH_cmb * A1 * Se;
L2 = RH_cmb * A2 * Se;
L3 = RH_cmb * A3 * Se;
L4 = RH_cmb * A4 * Se;

% plot
bode(L1, L2, L3, L4)
axs = get(gcf(), 'children')
grid()
legend(axs(2), 'classical A', 'loop-shaping', 'H-inf')
legend(axs(3), 'classical A', 'loop-shaping', 'H-inf')
xlim(axs(2), [1e-1, 1e3])
% ylim(axs(3), [-20, 20])
% ylim(axs(2), [-360, 0])
title('loop gain')