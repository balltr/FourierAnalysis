function [estimates,fitdata] = exponentialfit(x, ydata)

% Call this function like this [estimates,fitdata] = powerfitcurve(dx,ydata)
%
% dx should be your list of grid resolutions (some estimate of delta x of grid)
% starting with the biggest and moving to the finest
%
% ydata should be your list of calculated values for each mesh
%
% General form for curve fit
% ydata = C*exp(a*x)
%
% output is:
%	estimates(1): predicted exact value (y_exact in curve fit)
% 	estimates(2): error constant (C constant for curve fit)
%   esitmates(3): order of accuracy (p constant for curve fit)
% 
% 	fitdata: values of curve fit evaluated at dx locations
%

% This part solves the minimization problem for the curve fit
% Need initial guess for C and a
n = size(x,2);
if (n <= 2) 
	disp('Need more data to fit');
	return;
end

% CALCULATE FIRST GUESS WITH APPROXIMATE SOLUTION FROM RANGE PROVIDED
a = (log(ydata(n))-log(ydata(1)))/(x(n)-x(1));
C = ydata(1)/exp(a*x(1));
estimates = [C,a];
disp(['First guess for fit parameters: C = ' num2str(estimates(1)) ', a = ' num2str(estimates(2))]);


% GOING TO USE THE fminshearch FUNCTION TO FIT DATA
% FINDS THE CONSTANTS TO MINIMIZE THE ERROR IN FIT
estimates = fminsearch(@expfun, estimates);

disp(['Final answer for fit parameters: C = ' num2str(estimates(1)) ', a = ' num2str(estimates(2))]);

% CALCULATE datapoints of fit for return
[sse, fitdata] = expfun(estimates);


% exp accepts curve parameters as inputs, and outputs
% the sum of the square error for ydata - C*exp(a*x) (sse)
% it also outputs the FittedCurve data points so you can plot it
    function [sse, FittedCurve] = expfun(params)
        C = params(1);
        a = params(2);
        FittedCurve = C * exp(a*x);
        ErrorVector = log(FittedCurve) - log(ydata);
        sse = sum(ErrorVector .^ 2);
    end
end
