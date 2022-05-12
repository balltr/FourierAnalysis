function [dwdt,dxdt] = get_dwdt_dxdt(wd,x)
Nx = size(wd,2);
dwdt = zeros(3,Nx);
dxdt = get_dxdt(wd);
F = flux(wd,dxdt);
for i = 1:Nx
    dwdt(:,i) = 1/(x(i+1)-x(i)) * ( -wd(:,i)*(dxdt(i+1)-dxdt(i)) + F(:,i)-F(:,i+1) );
end
end