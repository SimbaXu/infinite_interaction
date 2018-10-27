%% Note

% In this script, measurement/sensor noise is accounted for, in
% addition to the design in Oct02_SISO_design.m

s = tf('s');
P_ = Oct02_SISO_model();

%% scaling
Df = 5; Dd = 0.0615384615;
Dm  = 2.0; Dq = 0.25;

Din = diag([Df, Dq]);
Dout = diag([Dd, Df, Dm]);

P = inv(Dout) * P_ * Din;
% agument qu as another output (this is for shaping KS)
[bn, an] = butter(3, 80, 'high', 's')
Gd = tf(bn, an);
P = [P(1,1), P(1, 2);
     P(2, 1), P(2, 2);
     0, Gd;    % <---- augmented output
     P(3, 1), P(3, 2)]

%% model for design
P(1, 1) = P(1, 1) - 1;
omegaP = 4; zetaP = 0.8;
wP = 1 / tf([1, 2 * zetaP * omegaP, 1e-1], [1, 2 * zetaP * omegaP, omegaP^2]); % performance: error weight
wI = (s * 0.25 + 1.0083)^3 / (s * 0.149 + 1)^3;  % robust stability: multiplicative err
wn = 1 / tf(2e-3);  % actual weight
% wn = 1 / tf(1);  % "loose" weight
%% H-inf design:
Pinf = [wP * P(1, 1), wP * P(1, 2); 
        wI * P(2, 1), wI * P(2, 2); 
        wn * P(3, 1), wn * P(3, 2);
        P(4, 1), P(4, 2)]
[A_, CL, GAM] = hinfsyn(Pinf, 1, 1, [1, 40], 'Display', 'On', 'Method', 'MAXE');
A_ = Dq * A_ * inv(Dm);  % scalj
% A_ = reduce(A_, 5);
Ad_ = c2d(tf(A_), 0.008);


%% Controllers performance
% admittance controller
figure(1)
A1_ = tf([1], [1, 6, 5]);
% A1_ = tf([1], [0.7, 10, 5]);
A1 = inv(Dq) * A1_ * Dm;
CL1 = lft(P, A1);
bodeplot(CL1(1, 1), 'r', CL1(2, 1), 'b', 1 / wP, 'r--', 1 / wI, 'b--', ...
         CL1(3, 1), 'c', 1 / wn, 'c--', tf([1, 2 * 3, 0], [1, 2 * 3, 3^2]), 'g--')
legend('S', 'T', '1/wP', '1/wI', 'KS', '1/wn', 'ideal performance')
grid on;
title('Admittance controller Performance')

% % loop shaping
% figure(2)
% A1_ = 1 / 5 * 1 / (0.3 * s + 1) / (0.5 * s + 1) * (0.032 * s + 1) * (0.0021 * s + 1);
% A1 = inv(Dq) * A1_ * Dm;
% CL2 = lft(P, A1);
% bodeplot(CL2(1, 1), 'r', CL2(2, 1), 'b', 1 / wP, 'r--', 1 / wI, 'b--', ...
%          tf([1, 2 * 3, 0], [1, 2 * 3, 3^2]), 'g--')
% legend('S', 'T', '1/wP', '1/wI', 'ideal performance')
% grid on;
% title('Performance: Loop shaping controller ')

% H inf
figure(3)
A = inv(Dq) * A_ * Dm;
CL3 = lft(P, A);
bodeplot(CL3(1, 1), 'r', CL3(2, 1), 'b', 1 / wP, 'r--', 1 / wI, 'b--', ...
         CL3(3, 1), 'c', 1 / wn, 'c--', tf([1, 2 * 3, 0], [1, 2 * 3, 3^2]), 'g--')
legend('S', 'T', '1/wP', '1/wI', 'KS', '1/wn', 'ideal performance')
grid on;
title('Performance: H-inf controller ')


% % compare steps responses
% figure(4)
% step(CL1(1, 1), CL2(1, 1), CL3(1, 1));
% legend('admittance', 'l-shaping', 'H-inf')