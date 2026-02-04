clear; close all; clc;

%% ---------------- Time Domain ----------------
tdata = readmatrix('t4-fmys-Vmys-fs100000-n4096.txt','NumHeaderLines',2);
time = tdata(:,1);
voltage = tdata(:,2);

figure('Name','Mystery Signal – Time Domain','Color','w');
plot(time, voltage, 'g','LineWidth',1.2);

xlabel('$t$ [s]', 'Interpreter','latex','FontSize',12,'Color','k');
ylabel('$e(t)$ [V]', 'Interpreter','latex','FontSize',12,'Color','k');
set(gca,'Color','w','Box','off','XColor','k','YColor','k',...
    'FontName','Arial','FontSize',11);

legend({'$f_s = 100~\mathrm{kHz}, n = 4096'}, ...
    'Interpreter','latex','Location','southwest','TextColor','k','Color','w','Box','on');

xline(0.00185, '--r', '$T_{sine} = 0.00185\mathrm{[s]}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xline(0.033, '--r', '$T_{square} = 0.033\mathrm{[s]}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

axis tight;
hold off;


%% ---------------- Frequency Domain ----------------
% Files
file1 = 'p5-fmys-Vmys-fs2000-n4096.txt';
file2 = 'p6-fmys-Vmys-fs1000-n4096.txt';

figure('Name','Mystery Signal – Power Spectra','Color','w');

data1 = readmatrix(file1,'NumHeaderLines',2);
f1 = data1(:,1); P1 = data1(:,2);
f1(1) = []; P1(1) = [];

subplot(2,1,1);
plot(f1, P1, 'Color', 'b', 'LineWidth', 1.4); hold on;

xline(30, '--r', '$30.0~\pm~0.5~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xline(545, '--r', '$540.0~\pm~0.5~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xlabel('$f$ [Hz]', 'Interpreter','latex','FontSize',12,'Color','k');
ylabel('$P$ [dB]', 'Interpreter','latex','FontSize',12,'Color','k');
title('(a)', 'Interpreter','latex','FontSize',12,'Color','k');

legend({'$f_s = 2~\mathrm{kHz}$, n = 4096'}, ...
    'Interpreter','latex','Location','northeast','TextColor','k','Color','w','Box','on');

set(gca,'Color','w','Box','off','XColor','k','YColor','k', ...
    'FontName','Arial','FontSize',11);
hold off;

data2 = readmatrix(file2,'NumHeaderLines',2);
f2 = data2(:,1); P2 = data2(:,2);
f2(1) = []; P2(1) = [];

subplot(2,1,2);
plot(f2, P2, 'Color', 'r', 'LineWidth', 1.4); hold on;

xline(30, '--b', '$30.0~\pm~0.3~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xline(90, '--b', '$90.0~\pm~0.3~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xline(150, '--b', '$150.0~\pm~0.3~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xline(455, '--b', '$455.0~\pm~0.3~\mathrm{Hz}$', ...
    'LabelOrientation', 'horizontal', 'Interpreter', 'latex', 'HandleVisibility','off');

xlabel('$f$ [Hz]', 'Interpreter','latex','FontSize',12,'Color','k');
ylabel('$P$ [dB]', 'Interpreter','latex','FontSize',12,'Color','k');
title('(b)', 'Interpreter','latex','FontSize',12,'Color','k');

legend({'$f_s = 1~\mathrm{kHz}$, n = 4096'}, ...
    'Interpreter','latex','Location','northeast','TextColor','k','Color','w','Box','on');

set(gca,'Color','w','Box','off','XColor','k','YColor','k', ...
    'FontName','Arial','FontSize',11);
hold off;
