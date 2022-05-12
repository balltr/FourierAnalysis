function [wd, x] = RK4(wd,x,dt)
[K1w,K1x] = get_dwdt_dxdt(wd,x);
[K2w,K2x] = get_dwdt_dxdt(wd + dt*K1w/2,x + dt*K1x/2);
[K3w,K3x] = get_dwdt_dxdt(wd + dt*K2w/2,x + dt*K2x/2);
[K4w,K4x] = get_dwdt_dxdt(wd + dt*K3w,x + dt*K3x);

wd = wd + 1/6*dt*(K1w + 2*K2w + 2*K3w + K4w);
x = x + 1/6*dt*(K1x + 2*K2x + 2*K3x + K4x);
end