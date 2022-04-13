close all
clear
clc

tbl1 = readtable('decayRates.xlsx');

numericalData = tbl1.numerical;
analyticalData = tbl1.analytical;
kdx = zeros(length(tbl1.kdx),1);
for i = 1:length(tbl1.kdx)
    kdx(i) = str2num(tbl1.kdx{i});  %#ok<ST2NM>
end

figure
p1 = plot(kdx,numericalData,'o',kdx,analyticalData,'x');
p1(1).LineWidth = 1.5;
p1(2).LineWidth = 1.5;
title('Wavenumber Decay Rates')
xlabel('kdx'),ylabel('Decay Rate')
legend('FEM','Fourier')%,'location','southeast')
grid on


tbl2 = readtable('oscillationRates.xlsx');

numericalData = tbl2.numerical;
analyticalData = tbl2.analytical;
kdx = zeros(length(tbl2.kdx),1);
for i = 1:length(tbl2.kdx)
    kdx(i) = str2num(tbl2.kdx{i});  %#ok<ST2NM>
end

figure
p2 = plot(kdx,numericalData,'o',kdx,analyticalData,'x');
p2(1).LineWidth = 1.5;
p2(2).LineWidth = 1.5;
title('Wavenumber Oscillation Rates')
xlabel('kdx'),ylabel('Oscillation Frequency [Hz]')
legend('FEM','Fourier')%,'location','southeast')
grid on