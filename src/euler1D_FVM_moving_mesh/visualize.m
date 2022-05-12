% Plots results from run.m
clc
clearvars

%% read data

resultsPath = '../euler1D_FVM_moving_mesh/Results_Vcol168/';
V = readmatrix('V_0_upwind_tantheta_0');
t = readmatrix(strcat(resultsPath,'time.txt'));
wdinit = readmatrix(strcat(resultsPath,'baseflow.txt'));
xdata = dataRead(resultsPath,'x');
wdata = dataRead(resultsPath,'rstrt');

%% calculate percent variation in IC for each variable - primitive variables

variation = zeros(3,1);
for var = 1:3
    variation(var) = max(abs(wdinit(var)-wdata{1}(:,var)))/wdinit(var);
end

fprintf('variation of p: %.2f%%\n',variation(1)*100)
fprintf('variation of u: %.2f%%\n',variation(2)*100)
fprintf('variation of RT: %.2f%%\n',variation(3)*100)


%% plot variables at specific time - primitive variables
% close all

% specify what time to plot
tplot = 0;

% find index of t closest to tplot
idx = find( abs(t-tplot) == min(abs(t-tplot)) );
tplot = t(idx);

figure('Position',[0 10000 700 500]);
timeplot = tiledlayout(2,2);
title(timeplot,'t = ',tplot)
nexttile
stairs(xdata{idx},cat(1,wdata{idx}(:,1),wdata{idx}(end,1)))
title('p')
nexttile
stairs(xdata{idx},cat(1,wdata{idx}(:,2),wdata{idx}(end,2)))
title('u')
nexttile
stairs(xdata{idx},cat(1,wdata{idx}(:,3),wdata{idx}(end,3)))
title('RT')


%% plot decay rate of variable

var = 3; % specify variable to plot
tRange = [0,7]; % specify range to fit in units of time

varMin = zeros(length(t),1);
for i = 1:length(t)
    varMin(i) = min(wdata{i}(:,var));
end
varMean = wdinit(var);
varDiff = abs(varMin-varMean);
tidx1 = find( abs(t-tRange(1)) == min(abs(t-tRange(1))) ); % find index closest to tRange(1)
tidx2 = find( abs(t-tRange(2)) == min(abs(t-tRange(2))) ); % find index closest to tRange(2)
tRange = t(tidx1:tidx2);
varDiffRange = varDiff(tidx1:tidx2);
[estimates,varDiffFit] = exponentialfit(tRange',varDiffRange');

figure
semilogy(t,varDiff)
grid on
hold on
semilogy(tRange,varDiffFit,'ro')

title('Decay of var ',var)
xlabel('time')
ylabel('')
legend('data','exp fit')

%% plot eigenvector

Vcol = 168;

figure
for var = 1:4
    realVtemp = real(V(var:4:end-1,Vcol)); % stores eigenvector for current variable
    nexttile
    plot(realVtemp)
    title('var ',var)
    xlabel('y')
end

%% plot shock position over time

xshock = zeros(size(xdata,2),1);
for i = 1:size(xdata,2)
    xshock(i) = xdata{1}(1);
end
figure('Position',[0 0 700 366]);
plot(t,xshock,'ko')
title('shock position vs time')
xlabel('t')
ylabel('x')