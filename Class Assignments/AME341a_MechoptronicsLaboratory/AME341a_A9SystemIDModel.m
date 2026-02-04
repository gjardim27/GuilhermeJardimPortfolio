clear; clc; close all;

% Prepare Data Import
numTrials = 6;

baseName = 'trial'; 
fileSuffix = '_5000hz_50000n.lvm';

% Manually Found Vibration Beginning or Time at x0

t_release_values = [...
    1.8046;...
    1.2574; ...
    1.1174; ...
    1.3496; ...
    1.6078; ...
    1.2438];

i0     = 50;             % first steady peak 
n_min  = 50;             % Plateau start
n_max  = 100;             % Plateau end

% Setup Data Storage 
zeta   = zeros(numTrials,1);
omega0 = zeros(numTrials,1);

% Zeta vs n Loop over all trials
for trial = 1:numTrials

    filename = [baseName num2str(trial) fileSuffix];

    % Load Data
    T = readtable(filename,'FileType','text');
    time   = T{:,1};
    signal = T{:,2};

    % Peak Detection
    [pks, locs] = findpeaks(signal, time, ...
        'MinPeakDistance', 0.025, ... 
        'MinPeakHeight', 0); % Positive Peaks Only

    valid = locs > t_release_values(trial); % Sets Valid Region
    pks   = pks(valid);
    locs  = locs(valid);

    % Zeta Calculation
    x0 = pks(i0); % Initial Peak Amplitude
    Nmax = length(pks) - i0;

    n_values    = (1:Nmax)';
    zeta_values = zeros(Nmax,1);

    for k = 1:Nmax
        xn = pks(i0 + k);
        zeta_values(k) = log(x0/xn)/(2*pi*k);
    end

    % Plateau Index
    idx = (n_values >= n_min) & (n_values <= n_max);

    zeta_trial_mean = mean(zeta_values(idx)); % Mean of Zeta Values in Plateau
    zeta_trial_std  = std(zeta_values(idx)); % Std. Dev. of Zeta Values in Plateau

    zeta(trial) = zeta_trial_mean; 

    % Omega Calculation
    Td = mean(diff(locs(i0:i0+5)));    % Averaged period across 5 peaks
    omega_d = 2*pi / Td; % Damping Frequency
    omega_0 = omega_d / sqrt(1 - zeta_trial_mean^2); % Fundamental Frequency

    omega0(trial) = omega_0;

    fprintf('Trial %d:  zeta = %.5f,  omega0 = %.2f rad/s\n', ...
            trial, zeta_trial_mean, omega_0);

    if trial == 1
        figure('Color','w');
        plot(n_values, zeta_values, 'bo'); hold on;
        plot(n_values(idx), zeta_values(idx), 'ro','MarkerFaceColor','r');
        yline(zeta_trial_mean,'r--','LineWidth',2);
        grid on;
        xlabel('n', 'FontSize', 14); 
        ylabel('\zeta_n','FontSize', 14);
        xlim([0 100])
        set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
        'XColor','k','YColor','k');
        L = legend('\zeta_n', ...
           'Plateau Region Values', ...
           ['Mean \zeta = ' num2str(zeta_trial_mean,'%0.2g')], ...
           'Location','best', ...
           'Box','on', 'Color','w','TextColor','k');
        set(L,'TextColor','k');
    end
end

% Final System ID Values
zeta_final = mean(zeta);
zeta_unc   = 2.571 * std(zeta);

omega0_final = mean(omega0);
omega0_unc   = 2.571 * std(omega0);

fprintf('Final Zeta     = %.6f ± %.6f\n', zeta_final, zeta_unc);
fprintf('Final Omega_0  = %.2f ± %.2f rad/s\n', omega0_final, omega0_unc);

% Second Order Model
T = readtable([baseName '6' fileSuffix],'FileType','text'); % 6 -> Representative Trace
time   = T{:,1};
signal = T{:,2};

[pks, locs] = findpeaks(signal, time, ...
    'MinPeakDistance', 0.025, ...
    'MinPeakHeight', 0);

valid = locs > t_release_values(trial);
pks   = pks(valid);
locs  = locs(valid);

x0_model = pks(i0);
t0 = locs(i0);

Td = mean(diff(locs(i0:i0+5)));
omega_d = 2*pi / Td;

time_adjustment = time - t0;

model = x0_model .* exp(-zeta_final*omega0_final*time_adjustment).*cos(omega_d*time_adjustment);

% Envelope Calculation
env = x0_model * exp(-zeta_final * omega0_final * time_adjustment);

x0_unc   = 1e-6;        % given uncertainty in x0

% Partial derivatives
denv_dx0   = env ./ x0_model;
denv_dzeta = -omega0_final .* time_adjustment .* env;
denv_dw0   = -zeta_final .* time_adjustment .* env;

% Kline–McClintock propagation
env_unc = sqrt( ...
    (denv_dx0   .* x0_unc  ).^2 + ...
    (denv_dzeta .* zeta_unc).^2 + ...
    (denv_dw0   .* omega0_unc  ).^2 );

% Upper and lower bounds
env_upper = env + env_unc;
env_lower = env - env_unc;

% Plot Model and Experimental Time Trace
figure('Color','w'); hold on;

plot(time_adjustment, signal, 'b');
plot(time_adjustment, model, 'r');

% Envelope
plot(time_adjustment,  env, 'k-', 'LineWidth', 2);

% Envelope Uncertainty
plot(time_adjustment,  env_upper, 'k:','LineWidth',1.5);
plot(time_adjustment,  env_lower, 'k:','LineWidth',1.5);

xlim([0 5]);
xlabel('t [s]', 'FontSize', 14);
ylabel('e_o [V]', 'FontSize', 14);
set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');

L = legend('Experimental', ...
       '2nd Order Model', ...
       'Nominal Envelope', ...
       'Envelope Uncertainty',...
       'Location','northeast','Box','off', 'Color','w','TextColor','k');
set(L, 'TextColor', 'k'); 

% Overplot
t_zoom = [4 4.3];     % Time Window
y_zoom = [0.14 0.2];  % Voltage Window

ax1 = gca;                    
ax2 = axes('Position',[0.71 0.16 0.18 0.18]); 
box on; hold on;

% Plot zoomed region
plot(time_adjustment, signal, 'b');
plot(time_adjustment, model, 'r');
plot(time_adjustment, env, 'k-', 'LineWidth', 2);
plot(time_adjustment, env_upper, 'k:','LineWidth',1.5);
plot(time_adjustment, env_lower, 'k:','LineWidth',1.5);

xlim(t_zoom)
ylim(y_zoom)

set(gca,'FontSize',10,'Color','w', ...
    'XGrid','off','YGrid','off','Box','on', ...
    'XColor','k','YColor','k');

% Zoomed in Plots for Phase Analysis

% Phase Difference 1

win1 = (time_adjustment >= 0) & (time_adjustment <= 0.3); % Choose Window

signal1 = signal(win1);
model1 = model(win1);
t1   = time_adjustment(win1);

[pks_exp1, t_exp1] = findpeaks(signal1, t1, 'MinPeakDistance', 0.025);
[pks_model1, t_model1] = findpeaks(model1, t1, 'MinPeakDistance', 0.025);

N1 = min(length(t_exp1), length(t_model1));
delta_t1 = t_exp1(1:N1) - t_model1(1:N1); 
delta_t1_mean = mean(delta_t1); 

T = 2*pi / omega_d;
phi1 = 360 * delta_t1_mean / T; % Calculate Phase Differnce

fprintf('Phase Difference (0–0.3 s) = %.3f degrees\n', phi1);


% Phase Difference 2

win2 = (time_adjustment >= 4.5) & (time_adjustment <= 4.8);

signal2 = signal(win2);
model2 = model(win2);
t2   = time_adjustment(win2);

[pks_exp2, t_exp2] = findpeaks(signal2, t2, 'MinPeakDistance', 0.025);
[pks_model2, t_model2] = findpeaks(model2, t2, 'MinPeakDistance', 0.025);

N2 = min(length(t_exp2), length(t_model2));
delta_t2 = t_exp2(1:N2) - t_model2(1:N2);
delta_t2_mean = mean(delta_t2);

phi2 = 360 * delta_t2_mean / T;

fprintf('Phase Difference (4.5–4.8 s) = %.3f degrees\n', phi2);

figure('Color', 'w'); hold on;
subplot(2,1,1)
plot(time_adjustment, signal, 'b'); hold on;
plot(time_adjustment, model, 'r');
title('(a)', 'Interpreter','latex','FontSize',14,'Color','k');
xlim([0 0.3]);
xlabel('t [s]', 'FontSize', 14);
ylabel('e_o [V]', 'FontSize', 14);
set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');

L = legend('Experimental', ...
       '2nd Order Model', ...
       'Location','Best','Box','on', 'Color','w','TextColor','k');
set(L, 'TextColor', 'k'); 

subplot(2,1,2)
plot(time_adjustment, signal, 'b'); hold on;
plot(time_adjustment, model, 'r');
title('(b)', 'Interpreter','latex','FontSize',14,'Color','k');
xlim([4.5 4.8]);
xlabel('t [s]', 'FontSize', 14);
ylabel('e_o [V]', 'FontSize', 14);
set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');