G1 = exp(- s * 0.036) / (0.0437 * s + 1) * (1.25 * s + 20) / (0.00218 * s + 1)
G2 = exp(- s * 0.03) * (-5.5e6 * s^2) / (s+1000)^2 / (s + 22.88) / (0.00218 * s + 1)

% bode
figure(1)
bode(G1, 'r-', G2, 'b--', A2*(G1 + G2), 'g')
xlim([1e-1, 1e2])
grid()

figure(2)
nyquist(G1, 'r-', G2, 'b--', G1 + G2, 'g', {1e-1, 2e1})
legend()