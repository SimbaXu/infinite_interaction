function P = Oct02_SISO_model(Larm, TR, TS, tauR1, tauR2, tauS, tauH, kH, bH, KR2)
% Return a linear model of the admittance control robot. 
% The model has the following form:

% [d  ] = [P11 P12] [fH]
% [fHT]   [P21 P22] [qu]
% [mm ]   [P31 P32]

% y:   position of the end-effector
% fHT: force output, does not have physical meaning, only use to
%      form the complementary sensitivity transfer function T.
% mm: effective torque measured,

% fH: force exerted by human,
% qu: control generated by the admittance controller,

s = tf('s');
% if no parameter is passed, the nominal model is returned
if nargin == 0
    Larm = 0.4;
    R1 = (-s + 55.56) / (s + 55.56) / (0.0437 * s + 1);
    H = (50 + 5 * s);
    Se = 2 * pi * 73 / (s + 2 * pi * 73);
    R2 = 0.475 * s^2 / (s / 100 + 1)^2 * 35 / (s + 35) * (-s + 66) ...
         / (s + 66);
else
    R1 = 1 / (TR * s + 1);
    R1.InputDelay = tauR1;

    R2 = KR2 * s^2 * (s / 100 + 1)^2 * 35 / (s + 35); 
    R2.InputDelay = tauR2;
    
    Se = 1 / (TS * s + 1);
    Se.InputDelay = tauS;

    H = kH + bH * s;
end

P = [0, R1 * Larm;
     0, (R1 * Larm * H + R2);
     Larm * Se, - (R1 * Larm * H + R2) * Se * Larm];
end



    

