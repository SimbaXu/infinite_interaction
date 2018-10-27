%% Select controller
A = A1;
%% uncertainties: unit value is nominal
delta1 = 1.0; 
kH = 100;
bH = 10;
tauH = 0.00;
%% Construct blocks for testing
R1 = 1 / (TR * s + 1);
R1.InputDelay = tauR1;
% cmd -> tau
R2 = - delta1 * 0.24 * 1e6 * s^2 / (s^2 + 2 * 1e3 * s + 1e6) / (0.0437 * s + 1)
R2.InputDelay = tauR2;

% Sensor
Se = 1 / (TS * s + 1);
Se.InputDelay = tauS;
% Human
H2 = LH ^ 2 * (kH + bH * s);
H2.InputDelay = tauH;
H1 = LH * kH * tf(1);
H1.InputDelay = tauH;
RH_cmb = R1 * H2 + R2;

%% Closed-Loop
L = RH_cmb * A * Se;
S = (1 + L) ^ (-1);
KS = pade(A * S);
T = L * S; % tf from rH to y
TfrH = H1 * S / LH; % tf from rH to f
figure(1); clf;
hold on;
bode(S, 'b-');
bode(1 / Wp1, 'b--');
bode(T, 'g');
bode(1 / Wp3, 'g--');
bode(TfrH, 'r-');
bode(1 / Wp2, 'r--');
bodeplot(KS, 'y');
margin(L);
allAxesInFigure = findall(gcf(),'type','axes')
xlim([1e-2, 1e4]);
ylim(allAxesInFigure(3), [-80, 20])
grid on;
legend();

% Time responses
TfrH1 = H1 * (1 + Se * RH_cmb * A1)^(-1) / LH;
TfrH2 = H1 * (1 + Se * RH_cmb * A2)^(-1) / LH;
TfrH3 = H1 * (1 + Se * RH_cmb * A3)^(-1) / LH;
TfrH4 = H1 * (1 + Se * RH_cmb * A4)^(-1) / LH;
figure(2); clf;
step(0.1 * TfrH1, 10)
hold on;
step(0.1 * TfrH2, 10)
step(0.1 * TfrH3, 10)
legend('classical A', 'loop-shaping', 'H-inf')
title('Step response');
% figure(3);
% nyquist(R * Se * H2);

figure(3); clf;
% ramp -- step response
t = (0:1e-3:5);
u = 0.2 * t.^2; u(1000:5001) = 0.2;
lsim(TfrH1, u, t);
hold on;
lsim(TfrH2, u, t);
lsim(TfrH3, u, t);

% step(1 / s * 0.1 * TfrH1, 5)
% hold on;
% step(1 / s * 0.1 * TfrH2, 5)
% step(1 / s * 0.1 * TfrH3, 5)

legend('classical A', 'loop-shaping', 'H-inf')
title('Ramp response');
