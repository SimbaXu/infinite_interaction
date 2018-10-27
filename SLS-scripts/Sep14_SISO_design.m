% This script is made while I was trying to design an admittance
% controller with the specification "ease of use" (via mapping from rh
% to q) and "robust stability" (via gain and phase margin).
% date: 14/9/2018

%% Model parameters
TR = 0.032; % See org-note: Closed-loop characteristics of the Denso VS-060
TS = 1 / (2 * pi * 62.5);  % filter at Nyquist freq
mH = 5;
kH = 10; % human dynamic coefficients
bH = 1;
tauR = 0.05;
tauS = 0;
tauH = 0;
%% Blocks
s = tf('s');
% Robot
R = 1 / (TR * s + 1);
R.InputDelay = tauR;
% Sensor
Se = 1 / (TS * s + 1);
Se.InputDelay = tauS;
% Human
H2 = mH * kH + mH * bH * s;
H2.InputDelay = tauH;
H1 = mH * kH * tf(1);
H1.InputDelay = tauH;

%% Controller
A = 1 / (0.05 * s ^ 2 + 1.0 * s + 5);  % desired static gain from f to displacement is 1/5
A = 1 / 5 * 1 / (0.3 * s + 1) / (0.1 * s + 1) * (0.032 * s + 1) * ...
    (0.0025 * s + 1)

%% Closed-loop tf
L = R * Se * H2 * A;
S = (1 + L) ^ (-1);
L1 = H1 * R * A * Se;
T = L1 * S;
figure(1); clf;
margin(L);
hold on;
bode(S);
bode(T, {0, 5e2});
legend()
figure(2);
step(S * L1);
stepinfo(S * L1)
title('Step response for human guidance.')
figure(3);
nyquist(R * Se * H2);


