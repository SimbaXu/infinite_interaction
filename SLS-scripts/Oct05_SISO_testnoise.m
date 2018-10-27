s = tf('s');
P_ = Oct02_SISO_model();
% noise weight
% [bn, an] = cheby1(3, 3, 80, 'high', 's')
[bn, an] = butter(3, 80, 'high', 's')
Gd = tf(bn, an);
Wn = Gd * 1.0 / 1e-4;
% augument qu as another output
P_ = [P_(1, 1), P_(1, 2);
      P_(2, 1), P_(2, 2);
      0, Wn * 1;
      P_(3, 1), P_(3, 2);
      P_(3, 1), P_(3, 2)]

%% Controllers performance
% admittance controller
figure(1)
A1_ = tf([1], [1, 6, 5]);
% A1_ = tf([1], [0.7, 10, 5]);
CL1 = lft(P_, A1_);
bodeplot(CL1(1, 1), 'r', CL1(2, 1), 'b', CL1(3, 1), 'g', CL1(4, 1), 'c')
legend('1-S(!=T)', 'T', 'KS*Wn', 'tau')
grid on;
title('Admittance controller Performance')


% H inf
figure(3)
CL3 = lft(P_, A_);
bodeplot(CL3(1, 1), 'r', CL3(2, 1), 'b', CL3(3, 1), 'g', CL3(4, 1), 'c')
legend('1-S(!=T)', 'T', 'KS*Wn', 'tau')
grid on;
title('Performance: H-inf controller ')
