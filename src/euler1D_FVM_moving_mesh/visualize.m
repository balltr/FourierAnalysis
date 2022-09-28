% Plots results from run.m
clc
clear

% read data
% NOTE - REMOVE ALL OTHER RESULTS FOLDERS FROM MATLAB PATH BEFORE CHANGING resultsPath

resultsPath = '../euler1D_FVM_moving_mesh/kdx_0/Results_Vcol2/';
% resultsPath = '../euler1D_FVM_moving_mesh/Results/';
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
fprintf('variation of v: %.2f%%\n',variation(2)*100)
fprintf('variation of RT: %.2f%%\n',variation(3)*100)


%% plot variables at specific time - primitive variables
close all

% 0: one tiled plot
% 1: one plot for each var
multPlot = 1;

% specify what time to plot
tplot = 0;

% find index of t closest to tplot
idx = find( abs(t-tplot) == min(abs(t-tplot)) );
tplot = t(idx);


if multPlot == 0
    figure('Position',[0 10000 700 500]);
    timeplot = tiledlayout(2,2);
    % title(timeplot,['t = ' num2str(tplot)])
    title(timeplot,'Initial Condition');
    nexttile
    stairs(xdata{idx},cat(1,wdata{idx}(:,1),wdata{idx}(end,1)),'k')
    grid on
    title('p')
    
    nexttile
    stairs(xdata{idx},-1*cat(1,wdata{idx}(:,2),wdata{idx}(end,2)),'k')
    grid on
    title('v')
    
    nexttile
    stairs(xdata{idx},cat(1,wdata{idx}(:,3),wdata{idx}(end,3)),'k')
    grid on
    title('RT')
else
    figure
    stairs(-1*xdata{idx},cat(1,wdata{idx}(:,1),wdata{idx}(end,1)),'k')
    grid on
    title('Initial p')
    xlabel('y')
    ylabel('p')
    hold on
    plot(-1*xdata{idx},wdinit(1)*ones(1,length(xdata{idx})),'r--');
    hold off
    legend('Perturbed Flow','Mean Flow','Location','southeast')
    
    figure
    stairs(-1*xdata{idx},cat(1,-1*wdata{idx}(:,2),-1*wdata{idx}(end,2)),'k')
    grid on
    title('Initial v')
    xlabel('y')
    ylabel('v')
    hold on
    plot(-1*xdata{idx},-1*wdinit(2)*ones(1,length(xdata{idx})),'r--');
    hold off
    legend('Perturbed Flow','Mean Flow','Location','northeast')
    
    figure
    stairs(-1*xdata{idx},cat(1,wdata{idx}(:,3),wdata{idx}(end,3)),'k')
    grid on
    title('Initial RT')
    xlabel('y')
    ylabel('RT')
    hold on
    plot(-1*xdata{idx},wdinit(3)*ones(1,length(xdata{idx})),'r--');
    hold off
    legend('Perturbed Flow','Mean Flow','Location','southeast')
end


%% plot decay rate of flow variable
% close all
clc
clearvars -except resultsPath t V wdata wdinit xdata
close all


% 1 to plot exponential fit
plotfit = 1;
% 1 to include legend
includeLegend = 0;

var = 2; % variable to plot
tRange = t(1:17); % range to plot

varMin = zeros(length(t),1);
for i = 1:length(t)
    varMin(i) = min(wdata{i}(:,var));
end
varMean = wdinit(var);
varDiff = abs(varMin-varMean);
fitRange = [0,tRange(end)]; % range to fit in units of time

tidx1 = find( abs(t-fitRange(1)) == min(abs(t-fitRange(1))) ); % find index closest to tRange(1)
tidx2 = find( abs(t-fitRange(2)) == min(abs(t-fitRange(2))) ); % find index closest to tRange(2)
fitRange = t(tidx1:tidx2);
varDiffRange = varDiff(tidx1:tidx2);
[estimates,varDiffFit] = exponentialfit(fitRange',varDiffRange');

figure
semilogy(tRange,varDiff(1:length(tRange)),'k')
grid on
if var == 1
    title('Decay of p Perturbation')
    ylabel('$||p-\overline{p}||_\infty$','Interpreter','latex')
elseif var == 2
    title('Decay of v Perturbation')
    ylabel('$||v-\overline{v}||_\infty$','Interpreter','latex')
elseif var == 3
    title('Decay of RT Perturbation')
    ylabel('$||RT-\overline{RT}||_\infty$','Interpreter','latex')
end
xlabel('t')

if plotfit == 1
    hold on
    semilogy(fitRange,varDiffFit,'r--','LineWidth',1)
%     text(11,10^-7,[num2str(estimates(1)) ' \cdot ' 'e^{' num2str(estimates(2)) '}'], 'FontSize',12)
    
    if includeLegend == 1
        if var == 1
            legend('Difference of p From Mean','Exponential Fit', 'FontSize',12)
        elseif var == 2
            legend('Difference of v From Mean','Exponential Fit', 'FontSize',12)
        elseif var == 3
            legend('Difference of RT From Mean','Exponential Fit', 'FontSize',12)
        end
    end
end



%% plot decay of shock position

fitRange = [2,10]; % specify range to fit in units of time


xshock = zeros(length(t),1);
for i = 1:length(t)
    xshock(i) = xdata{i}(1);
end
xDiff = abs(xshock);
tidx1 = find( abs(t-fitRange(1)) == min(abs(t-fitRange(1))) ); % find index closest to tRange(1)
tidx2 = find( abs(t-fitRange(2)) == min(abs(t-fitRange(2))) ); % find index closest to tRange(2)
fitRange = t(tidx1:tidx2);
xDiffRange = xDiff(tidx1:tidx2);
[estimates,xDiffFit] = exponentialfit(fitRange',xDiffRange');


figure
semilogy(t,xshock,'k')
grid on
hold on
semilogy(fitRange,xDiffFit,'r--','LineWidth',2)
title('Shock Position vs. Time')
xlabel('Time')
ylabel('Shock Position')
text(15,10^-8.75,[num2str(estimates(1)) ' \cdot ' 'e^{' num2str(estimates(2)) '}'], 'FontSize',12)


%% plot eigenvector

Vcol = 2;

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
