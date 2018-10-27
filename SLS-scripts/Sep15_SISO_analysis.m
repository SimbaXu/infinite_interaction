% this script analyzes admittance control for a single robot joint
% (the third joint)
% note: human dynamic model:
%   tau = rH kH L - q L^2 (kH + s bH)
% where L is the moment arm length from the point of contact to the
% joint. (kH, bH) is the catersian stiffness of the human
% arm. Should be
% date: 15/9/2018

%% Model parameters for controller design (H-infinity)
TR = 0.032; % See org-note: Closed-loop characteristics of the Denso VS-060
TS = 1 / (2 * pi * 73);  % filter at Nyquist freq
LH = 0.7;
kH = 80; % human dynamic coefficients
bH = 5;
tauR = 0.05;
tauS = 0;
tauH = 0;
%% the blocks
s = tf('s');
% Robot
R = 1 / (TR * s + 1);
R.InputDelay = 0;
R = 1 / (TR * s + 1) * pade(exp(- tauR * s), 1);
% Sensor
Se = 1 / (TS * s + 1);
Se.InputDelay = tauS;
% Human
H2 = LH ^ 2 * (kH + bH * s);
H2.InputDelay = tauH;
H1 = LH * kH * tf(1);
H1.InputDelay = tauH;
% block diagram
%                     rH          +----------+
%                     ----------->|   H1     |----------+
%                                 +----------+          |
%                                                       |
%   r   +      +----------+        +----------+         v +
%   ---->o-----|  R       +--------|   H2     |-------->o
%        ^-    +----------+        +----------+        -|
%        |                                              |
%        |                                              |  f
%        |     +----------+        +----------+         |
%        +-----|   A      |--------+   Se     +---------+
%              +----------+   	   +----------+

%% Hinf synthesis
Wp1 = (s / 2 + 8) / (s + 8 * 0.5);  % S < 1 / Wp1
Wp2 = (s / 500 + 0.14) / (s + 0.14 * 1e-4);  % T_frH < 1 / Wp2
Wp3 = (s / 0.02 + 100) / (s + 100 * 1);  % T < 1 / Wp3

G = [Wp1, - R * H2 / H1 * Wp1;
     H1 / LH * Wp2, - R * H2 / LH * Wp2;
     0, R * H2 / H1 * Wp3;
     H1 * Se, - R * H2 * Se]
[K, CL, GAM, INFO] = hinfsyn(G, 1, 1, 'GMIN', 3, 'GMAX', 100);
K =reduce(K, 3);

%% Controller selection
% 1. classic
A1 = tf([1], [1, 6, 5]);
% 2. loop-shaping
A2 = 1 / 5 * 1 / (0.3 * s + 1) / (0.5 * s + 1) * (0.032 * s + 1) * (0.0021 * s + 1)
% 3. H-infty
A3 = tf(K);

A = A3;

%% Blocks for testing
LH = 0.7;
% Robot
R = 1 / (TR * s + 1);
R.InputDelay = tauR;
% Sensor
Se = 1 / (TS * s + 1);
Se.InputDelay = tauS;
% Human
H2 = LH ^ 2 * (kH + bH * s);
H2.InputDelay = tauH;
H1 = LH * kH * tf(1);
H1.InputDelay = tauH;
%% Closed-Loop
L = R * Se * H2 * A;
S = (1 + L) ^ (-1);
T = L * S; % tf from rH to y
TfrH = H1 * S / LH; % tf from rH to f
figure(1); clf;
margin(L);
hold on;
bode(S, 'b-');
bode(1 / Wp1, 'b--');
bode(T, 'g');
bode(1 / Wp3, 'g--');
bode(TfrH, 'r-');
bode(1 / Wp2, 'r--');
xlim([10e-1, 10e2])
legend()

% Time responses
TfrH1 = H1 * (1 + R * Se * H2 * A1)^(-1) / LH;
TfrH2 = H1 * (1 + R * Se * H2 * A2)^(-1) / LH;
TfrH3 = H1 * (1 + R * Se * H2 * A3)^(-1) / LH;
figure(2); clf;
step(0.1 * TfrH1, 5)
hold on;
step(0.1 * TfrH2, 5)
step(0.1 * TfrH3, 5)
legend('classical A', 'loop-shaping', 'H-inf')
title('Step response');
% figure(3);
% nyquist(R * Se * H2);

Ad = c2d(A, 0.008);
figure(3); clf;
step(1 / s * 0.1 * TfrH1, 5)
hold on;
step(1 / s * 0.1 * TfrH2, 5)
step(1 / s * 0.1 * TfrH3, 5)
legend('classical A', 'loop-shaping', 'H-inf')
title('Ramp response');
