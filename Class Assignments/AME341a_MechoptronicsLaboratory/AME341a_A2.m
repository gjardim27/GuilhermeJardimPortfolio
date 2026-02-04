clear; clc;

%% Load Data 
table = readtable('temperature_data25.txt', ...
    'FileType','text','Delimiter','\t','ReadVariableNames',true, ...
    'CommentStyle','#');

t  = table{:,1}; % Time [h]
T  = table{:,2}; % Temperature [K]

Tavg = mean(T); % Average Temperature

%% Create Fit Curve

omega = 2*pi/24; % [rad/h], daily frequency
stdT = std(T); % Standard Deviation of Temperature
A = stdT*sqrt(2); % Amplitude 
tmodel = linspace(0, 24, 400); 
Tmodel = A*sin(omega*tmodel) + Tavg; % Sinusoidal Model

%% Plot Data
figure('Color','w','Position',[100 100 680 480]); 
hold on;

ax = gca;
ax.Color = 'w'; % White plot area
ax.XColor = 'k'; % Black x-axis
ax.YColor = 'k'; % Black y-axis
ax.Box   = 'off'; 
ax.Layer = 'top';
grid on;

% Raw Data Plot
plot(t, T, 'o', ...
    'MarkerSize',6, 'LineWidth',1.1, ...
    'MarkerEdgeColor',[0 0.45 0.74], ...   
    'MarkerFaceColor','none');
box off

% Sinusoidal fit
plot(tmodel, Tmodel, '--', 'Color','k', 'LineWidth',1.2);

% Average line
yline(Tavg, '-', 'Color', 'b', 'LineWidth', 1);

% Min/Max lines
yline(292, ':', 'Color', 'r')
yline(294, ':', 'Color', 'r')

% Set Plot Parameters
xlim([0 24]);
xticks(0:6:24);
ylim([291.5 294]);
yticks(291.5:0.5:294.0);
yticklabels( compose('%.1f', yticks) );
xlabel('$t[h]$', 'Interpreter','latex');
ylabel('$T[K]$', 'Interpreter','latex');
grid on; box on;

% Create Legend
lg = legend( ...
    '$T$', ...
    '$A \sin(\omega t) + C$', ...
    '$T_{\text{avg}}$', ...
    '$T_{\text{min/max}}$', ...
    'Location','southwest', ...
    'Interpreter','latex');
lg.Box = 'on';
lg.Color = [1 1 1];
lg.FontSize = 9;
lg.TextColor = 'k';    
lg.EdgeColor = 'k';    