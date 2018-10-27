%% parameters / model
Larm = 0.4;
TR = 0.0437; % See org-note: Closed-loop characteristics of the Denso VS-060
TS = 1 / (2 * pi * 73);  % filter at Nyquist freq
tauR1 = 0.036;
tauR2 = 0.03;
tauS = 0.0;  % sensor
tauH = 0.0;
kH = 80;
bH = 10;
deltaH = 1.1;
P = Sep26_SISO_model(Larm, TR, TS, tauR1, tauR2, tauS, tauH, kH, bH, deltaH);

scale = 10 / 0.3; % scale from 10 N -> 1 and 0.3 meter in deviation
                  % to 1;

%% closed-loop response
Gs = containers.Map;  % closed-loop tf: from fH -> d
Ls = containers.Map;  % loop gain
for key = controllers.keys
    G = lft(P, controllers(key{1}));
    Gs(key{1}) = G;
    Ls(key{1}) = P(2, 2) * controllers(key{1});
end

%% plots step response (10Nm step)
fig = figure(1); clf;
ax = subplot(2, 1, 1);
for key = controllers.keys
    stepplot(ax, Gs(key{1}), 10);  hold on;
end
grid;
legend(controllers.keys);
ylim(ax, [0, 0.05]);

%% ramp-then-step response
ax = subplot(2, 1, 2);
t = (0:1e-3:5);
u = 10 * t; u(1000:5001) = 10;
for key = controllers.keys
    lsimplot(ax, Gs(key{1}), u, t); hold on;
end
legend(controllers.keys);
ylim(ax, [0, 0.5]);
grid;

% %% frequency-response (force->position, scaled)
% figure(2); clf;
% scale = 1 / 0.03;
% for key = controllers.keys
%     bodeplot(1 - scale * Gs(key{1}));  hold on;
% end
% legend(controllers.keys);
% grid;
% xlim([1e-1, 1e2]);
% title('Scaled Error response')


% %% frequency-response (force->position, scaled)
% figure(3); clf;
% scale = 1 / 0.03;
% for key = controllers.keys
%     bodeplot(scale * Gs(key{1}));  hold on;
% end
% legend(controllers.keys);
% grid;
% xlim([1e-1, 1e2]);
% title('Input-Output response')


%% frequency-response (open loop margin)
figure(4); clf
margin(- Ls('H-inf')); hold on;
for key = controllers.keys
    if ~strcmp(key{1}, 'H-inf')
        bode(- Ls(key{1}));  hold on;
    end
end
legend(controllers.keys);
grid;
xlim([1e-1, 1e3]);
% title('Margin')