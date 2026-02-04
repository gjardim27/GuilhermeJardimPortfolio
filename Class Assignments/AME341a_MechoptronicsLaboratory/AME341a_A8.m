clear; clc;

% Collect Files
thisFolder = fileparts(mfilename('fullpath'));
dataFolder = fullfile(thisFolder, 'A8_Data_Files');
files = dir(fullfile(dataFolder, 'Station*.xlsx'));

N = numel(files);
dP       = zeros(N,1);
Q        = zeros(N,1);
Mass     = zeros(N,1);
Pressure = zeros(N,1);

% Constants
rho  = 1000; % [kg/m^3] Density of Water     
dt   = 30; % [s] Experiment Time       
mu   = 0.000890; % [Pas] Dynamic Viscovity of Water

%  Tube Geometry 
D_in  = 0.0229; % [in] Diameter
D     = D_in * 0.0254; % [m] Diameter
R     = D / 2; % [m] Radius
L     = 12 * 0.0254; % [m] Tube Length


for ii = 1:N
    fname    = files(ii).name;
    fullpath = fullfile(files(ii).folder, fname);

    % Pressure
    tankP_str = extractBetween(fname, "TankP_", "_");
    tankP_str = strrep(tankP_str, "p", ".");
    Pressure(ii) = str2double(tankP_str);

    % Mass
    mass_str = extractBetween(fname, "Mass_", ".xlsx");
    mass_str = strrep(mass_str, "p", ".");
    Mass_g   = str2double(mass_str);
    Mass(ii) = Mass_g;

    T = readtable(fullpath);

    % Time
    if isnumeric(T.Time)
        time_sec = T.Time * 86400;
    else
        time_sec = seconds(T.Time - T.Time(1));
    end

    PT01 = T.PT01_psig_;
    PT02 = T.PT02_psig_;

    dP(ii) = mean(PT01) - mean(PT02);
    
    mass_kg = Mass_g / 1000;
    Q(ii)   = mass_kg / (rho * dt);
end

% dP vs Q Figure
figure('Color','w'); hold on;

scatter(Q, dP, 50, [1 0.5 0], 'filled');
Q_th   = linspace(min(Q), max(Q), 200);
dP_th  = (8 * mu * L ./ (pi * R^4)) .* Q_th;
dP_th_psi = dP_th / 6894.76;

plot(Q_th, dP_th_psi, 'b-', 'LineWidth', 1.5);

xlabel('Q [m^3/s]');
ylabel('Δp [psi]');

set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');

Lh = legend('Orange Station Data','Hagen-Poiseuille Theory','Location','NorthWest','Box','off');
set(Lh, 'TextColor', 'k'); 

% Q for Re = 2300
Vcrit   = 2300 * mu / (rho * D);
Qcrit   = Vcrit * (pi * R^2);
minor_losses = 0.6e-6;

xline(Qcrit, '--k','LineWidth',1, 'HandleVisibility', 'off');
xline(minor_losses, '--k', 'LineWidth', 1, 'HandleVisibility', 'off')

text(Qcrit*1.02, min(dP) + 0.05*(max(dP)-min(dP)), ...
    'Re ≈ 2300', 'Color','k','FontSize',9);

text(minor_losses * 1.05, ...              
     max(dP)*0.95, ...  
     'Minor losses', ...
     'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left');

text(max(Q)*0.7, max(dP)*0.95, 'Minor losses + turbulent region', ...
    'Color','b','FontSize',10,'FontWeight','bold');

hold off;


% Time Trace

repName = 'Station03_TankP_27p8_Mass_020p35Q.xlsx'; 
repPath = fullfile(dataFolder, repName);
T_rep   = readtable(repPath);

if isnumeric(T_rep.Time)
    t_rep = T_rep.Time * 86400;
else
    t_rep = seconds(T_rep.Time - T_rep.Time(1));
end

figure('Color','w'); hold on;

plot(t_rep, T_rep.PT01_psig_, 'b-', 'LineWidth',1.2);
plot(t_rep, T_rep.PT02_psig_, 'r-', 'LineWidth',1.2);

xlabel('t [s]');
ylabel('p [psi]');

set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');

L = legend('PT01 (upstream)','PT02 (downstream)','Location','Best','Box','off');
set(L, 'TextColor', 'k');  

xline(10, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(30, '--k', 'LineWidth', 1, 'HandleVisibility','off');
 
t_mid = (10 + 30) / 2;
y_top = max([T_rep.PT01_psig_; T_rep.PT02_psig_]) + ...
           0.03 * range([T_rep.PT01_psig_; T_rep.PT02_psig_]);

text(t_mid, y_top, 'Steady-State Region', ...
    'HorizontalAlignment','center','Color','k','FontSize',10);

hold off;
