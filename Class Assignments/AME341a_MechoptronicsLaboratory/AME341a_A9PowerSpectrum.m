clear; clc; close all;


T = readtable('ps_8_5000hz_5000n.lvm','FileType','text'); % Import Chosen Power Spectra
frequency   = T{:,1};
amplitude = T{:,2};

% Plot Power Spectra
figure('Color', 'w')
plot(frequency,amplitude)

xlim([0 1000]);
ylim([-120 0])
xlabel('f [Hz]', 'FontSize',14);
ylabel('|H| [-]', 'FontSize',14);
set(gca,'Color','w','XGrid','off','YGrid','off','Box','off',...
    'XColor','k','YColor','k');

xline(170.2/(2*pi), 'c--', 'LineWidth', 1) % System ID Fundamental Frequency

xline(27, 'r--', 'LineWidth', 1) % First Peak Frequency

xline(168, 'g--', 'LineWidth', 1) % Second Peak Frequency

xline(471, 'k--', 'LineWidth', 1) % Third Peak Frequency

xline(917, 'm--', 'LineWidth', 1) % Fourth Peak Frequency

L = legend('Power Spectrum', ...
       'f_{0,ID} = 27.1 Hz', ...
       'f_0 = 27 Hz', ...
       'f_1 = 168 Hz', ...
       'f_2 = 471 Hz',...
       'f_3 = 917 Hz',...   
       'Location','Best','Box','on', 'Color','w','TextColor','k');

set(L, 'TextColor', 'k'); 