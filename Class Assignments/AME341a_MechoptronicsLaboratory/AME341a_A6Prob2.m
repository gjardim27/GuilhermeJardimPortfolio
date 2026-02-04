clear; close all; clc;

%% ---------------- Power Spectra ----------------
% Files
p3 = readmatrix('p3-f2000-V2-fs10000-n2048.txt','NumHeaderLines',2);  
p4 = readmatrix('p4-f2000-V2-fs2500-n2048.txt','NumHeaderLines',2);   
f3 = p3(:,1); P3 = p3(:,2);
f4 = p4(:,1); P4 = p4(:,2);

% Plot 
figure('Name','Aliasing Power Spectra','Color','w');
plot(f3, P3, 'b', 'LineWidth', 1.4); hold on;
plot(f4, P4, 'r', 'LineWidth', 1.4);

% Annotate
xline(2000, '--b', '$2000~\pm~5~\mathrm{Hz}$', ...
    'Interpreter','latex','FontSize',10, ...
    'LabelVerticalAlignment','top', 'LabelOrientation','horizontal', ...
    'HandleVisibility','off');
xline(500, '--r', '$500~\pm~2~\mathrm{Hz}$', ...
    'Interpreter','latex','FontSize',10, ...
    'LabelVerticalAlignment','top', 'LabelOrientation','horizontal', ...
    'HandleVisibility','off');


xlabel('$f$ [Hz]', 'Interpreter','latex','FontSize',12,'Color','k');
ylabel('Power [dB]', 'Interpreter','latex','FontSize',12,'Color','k');
legend({'$f_s = 10~\mathrm{kHz},~n=2048$', ...
        '$f_s = 2.5~\mathrm{kHz},~n=2048$'}, ...
        'Interpreter','latex','Location','best','Color','w','TextColor','k');
set(gca, 'Color','w','Box','off','XColor','k','YColor','k', ...
         'FontName','Arial','FontSize',11);
axis tight;
hold off;


%% ---------------- Time-Domain Traces ----------------
t1 = readmatrix('t1-f2000-V2-fs60000-n1024.txt','NumHeaderLines',2);
t2 = readmatrix('t2-f2000-V2-fs2000-n1024.txt','NumHeaderLines',2);
t3 = readmatrix('t3-f2001-V2-fs2000-4096.txt','NumHeaderLines',2);

figure('Name','Aliasing in Time Domain','Color','w');

% (a)
subplot(3,1,1)
plot(t1(:,1), t1(:,2), 'b','LineWidth',1.2);
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',11,'Color','k');
ylabel('$e(t)$ [V]', 'Interpreter','latex','FontSize',11,'Color','k');
legend({'$f_s = 60~\mathrm{kHz},~n=1024$'}, 'Interpreter','latex', ...
       'Location','best','Color','w','TextColor','k');
title('(a)', 'Interpreter','latex','FontSize',12,'Color','k');
set(gca, 'Color','w','Box','off','XColor','k','YColor','k', ...
         'FontName','Arial','FontSize',10);

% (b)
subplot(3,1,2)
plot(t2(:,1), t2(:,2), 'r','LineWidth',1.2);
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',11,'Color','k');
ylabel('$e(t)$ [V]', 'Interpreter','latex','FontSize',11,'Color','k');
legend({'$f_s = 2~\mathrm{kHz},~n=1024$'}, 'Interpreter','latex', ...
       'Location','best','Color','w','TextColor','k');
title('(b)', 'Interpreter','latex','FontSize',12,'Color','k');
set(gca, 'Color','w','Box','off','XColor','k','YColor','k', ...
         'FontName','Arial','FontSize',10);

% (c)
subplot(3,1,3)
plot(t3(:,1), t3(:,2), 'g','LineWidth',1.2);
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',11,'Color','k');
ylabel('$e(t)$ [V]', 'Interpreter','latex','FontSize',11,'Color','k');
legend({'$f_s = 2~\mathrm{kHz},~n=4096$'}, 'Interpreter','latex', ...
       'Location','best','Color','w','TextColor','k');
title('(c)', 'Interpreter','latex','FontSize',12,'Color','k');
set(gca, 'Color','w','Box','off','XColor','k','YColor','k', ...
         'FontName','Arial','FontSize',10);
