%% Note

% In the previous design script, I was trying to shape the sensitivity
% transfer function S and the complementary transfer function T to
% design H-inf controllers. This approach does not seem to work very
% well.

% In this design script, I will try to optimize directly the transfer
% function from the input fH to the error which is defined as the
% robot's current position w.r.t to its desired value. All quantities
% are properly scaled.

s = tf('s');

%% Parameters

Wpe = (s / sqrt(1.2) + 0.5)^2 / (s + 0.5 * sqrt(0.05))^2; % error weight
Wpd = (s / sqrt(0.01) + 60)^2 / (s * 1e-2 + 60 * 1.2)^2;  % displacement weight


%% Nominal model
Larm = 0.4;
TR = 0.0437; % See org-note: Closed-loop characteristics of the Denso VS-060
TS = 1 / (2 * pi * 73);  % filter at Nyquist freq
tauR1 = 0.036;
tauR2 = 0.03;
tauS = 0.0;  % sensor
tauH = 0;
kH = 80;
bH = 10;
deltaH = 0.2;
P = Sep26_SISO_model(Larm, TR, TS, tauR1, tauR2, tauS, tauH, kH, bH, deltaH);

%% Scaling
Df = 5; Dd = 1 / 0.106;
Dq = 0.3; Dm = 1 / 2.5;

Din = diag([Df, Dq]);
Dout = diag([Dd, Dm]);

Phat_ = Dout * tf(pade(P)) * Din;
%% Model for H-inf design: control output is error:= dhat - fHhat
Phat = [
    Wpe * (Phat_(1, 1) - 1), Wpe * Phat_(1, 2);
    Wpd * Phat_(1, 1), Wpd * Phat_(1, 2);
    Phat_(2, 1), Phat_(2, 2)
       ]

%% H-inf
[A, CL, GAM] = hinfsyn(Phat, 1, 1, [1, 40], 'Display', 'On', 'Method', 'MAXE');


%% controllers
A1 = tf([1], [1, 6, 5]);
% 2. loop-shaping
A2 = 1 / 5 * 1 / (0.3 * s + 1) / (0.5 * s + 1) * (0.032 * s + 1) * (0.0021 * s + 1);
% 3. H-infty
A3 = Dq * tf(A) * Dm;

controllers = containers.Map;
controllers('classic-A') = A1;
% controllers('loop-shaping') = A2;
controllers('H-inf') = reduce(A3, 5);

figure(10);clf;
Phat_ = pade(Dout * P * Din);
stepplot(lft(Phat_, A1 / (Dq * Dm)), lft(Phat_, A3 / (Dq * Dm)), 10)
legend('classic-A', 'H-inf')

figure(11); clf;
Phat_scaled = [
    Phat_(1, 1) - 1, Phat_(1, 2);
    Phat_(1, 1), Phat_(1, 2);
    Phat_(2, 1), Phat_(2, 2)]
CL3_scaled = lft(Phat_scaled, A3 / (Dq * Dm));
CL1_scaled = lft(Phat_scaled, A1 / (Dq * Dm));
bodeplot(CL3_scaled(1, 1), 'r-o', CL3_scaled(2, 1), 'b-o', ...
         CL1_scaled(1, 1), 'r-', CL1_scaled(2, 1), 'b-', ...
         1 / Wpe, 'r--', 1 / Wpd, 'b--');
legend('e:=fH - d', 'd')

% figure(12); clf;
% CL0 = lft(Phat, 0);
% CL3 = lft(Phat, A3 / (Dq * Dm));
% CL1 = lft(Phat, A1 / (Dq * Dm));
% sigmaplot(CL1, 'r', CL3, 'g', tf(GAM));
% legend('classic-A', 'H-inf-man')

