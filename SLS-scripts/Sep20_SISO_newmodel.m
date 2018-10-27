% this script designs and analyzes H-infinity controller for
% admittance control.  

% novelty: the mapping from robot command/motion to ft sensor reading
% is now considered explicitly. See org_imgs/Admittance-model-new.png
% for a drawing of this new model.

%% Model parameters for controller design (H-infinity)
TR = 0.0437; % See org-note: Closed-loop characteristics of the Denso VS-060
TS = 1 / (2 * pi *73);  % filter at Nyquist freq
LH = 0.4;
kH = 80; % human dynamic coefficients
bH = 10; % cannot set to one, badly conditioned
tauR1 = 0.036;
tauR2 = 0.03;
tauS = 0.0;
tauH = 0;

%% blocks
s = tf('s');
% Robot: cmd -> position
R1 = 1 / (TR * s + 1);
R1.InputDelay = tauR1;
R1 = pade(R1, 1);
% cmd -> tau
R2 = - 1.0 * 0.24 * 1e6 * s^2 / (s^2 + 2 * 1e3 * s + 1e6) / (0.0437 * s + 1)
R2.InputDelay = tauR2;
R2 = pade(R2, 1);

% Sensor
Se = 1 / (TS * s + 1);
Se.InputDelay = tauS;
% Human
H2 = LH ^ 2 * (kH + bH * s);
H2.InputDelay = tauH;
H1 = LH * kH * tf(1);
H1.InputDelay = tauH;

RH_cmb = R1 * H2 + R2;
% RH_cmb = R1 * H2;  % wout considering position -> tau

%% Hinf synthesis
Wp1 = (s / 2^(1/3) + 3.0)^3 / (s + 3.0 * 0.16^(1/3))^3;  % S < 1 / Wp1
Wp2 = (s / 1000 + 0.0001) / (s + 0.0001 * 20);  % T_frH < 1 / Wp2
Wp3 = (s / sqrt(0.3) + 60)^2 / (s + 60 * sqrt(1.2))^2;  % T < 1 / Wp3

G = [Wp1, - RH_cmb / H1 * Wp1;
     % H1 / LH * Wp2, - RH_cmb / LH * Wp2;
     0, RH_cmb / H1 * Wp3;
     H1 * Se, - RH_cmb * Se
    ];
[K, CL, GAM] = hinfsyn(G, 1, 1, [2.5,10], 'Display', 'On', 'Method', 'lmi');
% [K, CL, GAM] = h2syn(G, 1, 1);
figure(4); sigmaplot(CL, tf(GAM)); xlim([1e-3, 1e4]);
K =reduce(K, 5);

%% Controller selection
% 1. classic
A1 = tf([1], [1, 6, 5]);
% 2. loop-shaping
A2 = 1 / 5 * 1 / (0.3 * s + 1) / (0.5 * s + 1) * (0.032 * s + 1) * (0.0021 * s + 1);
% 3. H-infty
A3 = tf(K);
A = A3;

% 4. imc
A4 = (s + 458) * (s + 23) / 13121 / (s + 16) / ((0.3 * s + 1) / 0.8 - exp(- ...
                                                  0.036 * s));

controllers = containers.Map;
controllers('classic-A') = A1;
controllers('loop-shaping') = A2;
controllers('H-inf') = A3;
controllers('imc') = A4;

Ad = c2d(A3, 0.008);

