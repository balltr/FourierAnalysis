function [estimates,fitdata] = powerfitcurve(dx, ydata)

% Call this function like this [estimates,fitdata] = powerfitcurve(dx,ydata)
%
% dx should be your list of grid resolutions (some estimate of delta x of grid)
% starting with the biggest and moving to the finest
%
% ydata should be your list of calculated values for each mesh
%
% General form for curve fit
% ydata = y_exact +C dx^p
%
% output is:
%	estimates(1): predicted exact value (y_exact in curve fit)
% 	estimates(2): error constant (C constant for curve fit)
%   esitmates(3): order of accuracy (p constant for curve fit)
% 
% 	fitdata: values of curve fit evaluated at dx locations
%

% This part solves the minimization problem for the curve fit
%Need initial guess for y_exact, C, p
n = size(dx,2);
if (n <= 2) 
	disp('Need more data to fit');
	return;
end

% CALCULATE FIRST GUESS WITH APPROXIMATE SOLUTION FROM LAST 3 POINTS (ASSUMES LARGE dx RATIO)
p = log((ydata(n)-ydata(n-1))/(ydata(n-1)-ydata(n-2)))/log(dx(n-1)/dx(n-2));
C = (ydata(n-1)-ydata(n))/dx(n-1)^p;
y_exact = ydata(n)-C*dx(n)^p;
estimates = [y_exact,C,p];
disp(['First guess for fit parameters: y_exact = ' num2str(estimates(1)) ', C = ' num2str(estimates(2)) ', p = ' num2str(estimates(3))]);


% GOING TO USE THE fminshearch FUNCTION TO FIT DATA
% FINDS THE THREE CONSTANTS TO MINIMIZE THE ERROR IN FIT
estimates = fminsearch(@powerfun, estimates);

disp(['Final answer for fit parameters: y_exact = ' num2str(estimates(1)) ', C = ' num2str(estimates(2)) ', p = ' num2str(estimates(3))]);

% CALCULATE datapoints of fit for return
[sse, fitdata] = powerfun(estimates);


% powefun accepts curve parameters as inputs, and outputs
% the sum of the square error for y_exact +C dx.^p - ydata (sse) 
% it also outputs the FittedCurve data points so you can plot it
    function [sse, FittedCurve] = powerfun(params)
        y_exact = params(1);
        C = params(2);
        p = params(3);
        FittedCurve = y_exact +C*dx.^p;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end
