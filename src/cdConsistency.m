clc
clear
close all

gamma = 1.4;

% downstream
pd = 10.333333333333334;
rhod = 3.857142857142857;
vd = 0.920279077371051;
ud = 0.354964786985977;
RTd = 2.679012345679012;
Ed = 7.183987654320989;

wd = [rhod
      rhod*ud
	  rhod*vd
	  rhod*Ed];
barwd = wd;

% perturbations
wdprime = [0.0
           0.0
           0.0
           0.0];
wd = wd+wdprime;



cdexact = sqrt( gamma*(gamma-1)/wd(1) * (wd(4) - ((wd(2))^2 + (wd(3))^2)/(2*(wd(1)))) )

barcd = sqrt( gamma*(gamma-1)/barwd(1) * (barwd(4) - ((barwd(2))^2 + (barwd(3))^2)/(2*(barwd(1)))) );
w1coeff = wdprime(1) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1)^2 * ( (barwd(2)^2+barwd(3)^2)/barwd(1) - barwd(4) ) );
w2coeff = wdprime(2) * gamma*(gamma-1)/(2*barcd) * ( -barwd(2)/barwd(1)^2 );
w3coeff = wdprime(3) * gamma*(gamma-1)/(2*barcd) * ( -barwd(3)/barwd(1)^2 );
w4coeff = wdprime(4) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1) );
cd = barcd + w1coeff + w2coeff + w3coeff + w4coeff

ratio = cd/cdexact



