A = [0.9 0.5; 0 0.99];
B = [1 0; 5 1];
C = [1 1; 1 0];
D = [0 1; 0 0];
sys = ss(A, B, C, D, 1);
[K, CL, GAM] = h2syn(sys, 1, 1);
step(CL);

B1 = B(:, 1);
B2 = B(:, 2);
C1 = C(1, :);
C2 = C(2, :);
D11 = 0;
D12 = 1;
D21 = 0;
D22 = 0;

z = tf('z', 1);
R = (z * eye(2) - A - B2 * K * C2)^(-1);
M = K * C2 * R;
N = R * B2 * K;
L = K + K * C2 * R * B2 * K;

Rs = impulse(R, 10);
Ms = impulse(M, 10);
Ns = impulse(N, 10);
Ls = impulse(L, 10);