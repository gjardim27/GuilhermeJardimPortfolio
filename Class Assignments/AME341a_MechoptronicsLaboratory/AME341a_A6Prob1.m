clear; close all; clc;

%% ---------------- Power Spectra 1 (Sine Wave) ----------------
figure('Name','Sine Wave Spectra', 'Color', 'w'); hold on;

% Files
p1_files = {...
    'p1.1-f200-V2-fs1000-n512.txt',...
    'p1.2-f200-V2-fs1000-n512.txt',...
    'p1.3-f200-V2-fs1000-n512.txt',...
    'p1.4-f200-V2-fs1000-n512.txt',...
    'p1.5-f200-V2-fs1000-n512.txt'
    };

% Plots
for i = 1:1
    data = readmatrix(p1_files{i}, 'NumHeaderLines', 2);
    f = data(:,1); P = data(:,2);
    f(1) = []; P(1) = []; 
    plot(f, P, 'LineWidth', 1.2, 'Color', 'r');
end

xlabel('$f$ [Hz]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'k');
ylabel('$P$ [dB]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'k');
lgd = legend('$f_s = 1000 \mathrm{[Hz]}, n = 512$', ...
    'Interpreter','latex','Location','best');
set(lgd, 'Color', 'w', 'Box', 'on', 'TextColor', 'k');

set(gca, 'Color', 'w', 'Box', 'off', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Arial', 'FontSize', 11);

peakF = zeros(1, numel(p1_files));
peakP = zeros(1, numel(p1_files));

for i = 1:numel(p1_files)
    data = readmatrix(p1_files{i}, 'NumHeaderLines', 2);
    [pk, idx] = max(data(:,2));
    peakF(i) = data(idx,1);
    peakP(i) = pk;
end

N = numel(peakF);
tcrit = tinv(0.975, N-1);  
freq_mean = mean(peakF);
freq_unc = tcrit * std(peakF) / sqrt(N);
pow_mean = mean(peakP);
pow_unc = tcrit * std(peakP) / sqrt(N);

fprintf('Sine Wave → Mean Freq = %.2f ± %.2f Hz (95%% CI),  Mean Power = %.2f ± %.2f dB (95%% CI)\n', ...
    freq_mean, freq_unc, pow_mean, pow_unc);

xline(200, '--b', '$200~\pm~2~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');
hold off;


%% ---------------- Square Wave (p2.1–p2.5) ----------------
figure('Name','Square Wave Spectra','Color','w'); hold on;
p2_files = {...
    'p2.1-f200-V2-fs3000-n1024-SQR.txt',...
    'p2.2-f200-V2-fs3000-n1024-SQR.txt',...
    'p2.3-f200-V2-fs3000-n1024-SQR.txt',...
    'p2.4-f200-V2-fs3000-n1024-SQR.txt',...
    'p2.5-f200-V2-fs3000-n1024-SQR.txt',...
    };

% Plots
for i = 1:1
    data = readmatrix(p2_files{i}, 'NumHeaderLines', 2);
    f = data(:,1); P = data(:,2);
    f(1) = []; P(1) = [];
    plot(f, P, 'LineWidth', 1.2, 'Color', 'b');
end

xlabel('$f$ [Hz]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'k');
ylabel('$P$ [dB]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', 'k');
lgd2 = legend('$f_s = 3000 \mathrm{[Hz]}, n = 1024$','$\mathrm{Test~2}$','$\mathrm{Test~3}$','$\mathrm{Test~4}$','$\mathrm{Test~5}$', ...
    'Interpreter','latex','Location','southwest');
set(lgd2, 'Color', 'w', 'Box', 'on', 'TextColor', 'k');

set(gca, 'Color', 'w', 'Box', 'off', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Arial', 'FontSize', 11);

hold on;
for h = 1:2:7
    xline(200*h, '--r', sprintf('$%d~\\pm~3~\\mathrm{Hz}$', 200*h), ...
        'LabelOrientation','horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');
end
hold off;

peakF2 = zeros(1, numel(p2_files));
peakP2 = zeros(1, numel(p2_files));

for i = 1:numel(p2_files)
    data = readmatrix(p2_files{i}, 'NumHeaderLines', 2);
    [pk, idx] = max(data(:,2));
    peakF2(i) = data(idx,1);
    peakP2(i) = pk;
end

N2 = numel(peakF2);
tcrit2 = tinv(0.975, N2-1);
freq_mean2 = mean(peakF2);
freq_unc2 = tcrit2 * std(peakF2) / sqrt(N2);
pow_mean2 = mean(peakP2);
pow_unc2 = tcrit2 * std(peakP2) / sqrt(N2);

fprintf('Square Wave → Mean Freq = %.2f ± %.2f Hz (95%% CI),  Mean Power = %.2f ± %.2f dB (95%% CI)\n', ...
    freq_mean2, freq_unc2, pow_mean2, pow_unc2);
