% This is for a 6x6 B0 grid

close all
clear
clc

addpath('/Users/Tristan/Packages/share/petsc/matlab')

%  [U,S,V] = svds(J,1,'smallest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% central dtinv = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Jcen = PetscBinaryRead('Jacobian_cen');

Jcenfull = full(Jcen);
Jcenrow = sparse(Jcen(42,:))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% upwind = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Jup = PetscBinaryRead('Jacobian_up');

Jupfull = full(Jup);
Juprow = sparse(Jup(42,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M central %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JMcene22 = PetscBinaryRead('Jacobian_Mcen_e22');
JMcene15 = PetscBinaryRead('Jacobian_Mcen_e15');
JMcene10 = PetscBinaryRead('Jacobian_Mcen_e10');
JMcene9 = PetscBinaryRead('Jacobian_Mcen_e9');
JMcene8 = PetscBinaryRead('Jacobian_Mcen_e8');
JMcene7 = PetscBinaryRead('Jacobian_Mcen_e7');
JMcene6 = PetscBinaryRead('Jacobian_Mcen_e6');
JMcene5 = PetscBinaryRead('Jacobian_Mcen_e5');

% JMcene22full = full(JMcene22);
% JMcene22row = sparse(JMcene22(42,:))
% 
% JMcene15full = full(JMcene15);
% JMcene15row = sparse(JMcene15(42,:))
% 
% JMcene10full = full(JMcene10);
% JMcene10row = sparse(JMcene10(42,:))
% 
% JMcene9full = full(JMcene9);
% JMcene9row = sparse(JMcene9(42,:))
% 
% JMcene8full = full(JMcene8);
% JMcene8row = sparse(JMcene8(42,:))
% 
JMcene7full = full(JMcene7);
JMcene7row = sparse(JMcene7(42,:))
% 
% JMcene6full = full(JMcene6);
% JMcene6row = sparse(JMcene6(42,:))
% 
% JMcene5full = full(JMcene5);
% JMcene5row = sparse(JMcene5(42,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%% M upwind %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JMupe22 = PetscBinaryRead('Jacobian_Mup_e22');
JMupe15 = PetscBinaryRead('Jacobian_Mup_e15');
JMupe10 = PetscBinaryRead('Jacobian_Mup_e10');
JMupe9 = PetscBinaryRead('Jacobian_Mup_e9');
JMupe8 = PetscBinaryRead('Jacobian_Mup_e8');
JMupe7 = PetscBinaryRead('Jacobian_Mup_e7');
JMupe6 = PetscBinaryRead('Jacobian_Mup_e6');
JMupe5 = PetscBinaryRead('Jacobian_Mup_e5');

% JMupe22full = full(JMupe22);
% JMupe22row = sparse(JMupe22(42,:))
% 
% JMupe15full = full(JMupe15);
% JMupe15row = sparse(JMupe15(42,:))
% 
% JMupe10full = full(JMupe10);
% JMupe10row = sparse(JMupe10(42,:))
% 
% JMupe9full = full(JMupe9);
% JMupe9row = sparse(JMupe9(42,:))
% 
% JMupe8full = full(JMupe8);
% JMupe8row = sparse(JMupe8(42,:))
% 
JMupe7full = full(JMupe7);
JMupe7row = sparse(JMupe7(42,:))
% 
% JMupe6full = full(JMupe6);
% JMupe6row = sparse(JMupe6(42,:))
% 
% JMupe5full = full(JMupe5);
% JMupe5row = sparse(JMupe5(42,:))

%%%%%%%%%%%%%%%%%%%%%% eig analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mat_up = inv(JMupe7full/1e7)*-Jupfull;
% [V,E_mat] = eig(mat_up);
[V,E_mat] = eig(-Jupfull,JMupe7full/1e7);
E_vec = diag(E_mat);
figure
plot(real(E_vec),imag(E_vec),'marker','x','LineStyle','none')
title('dtinv = 1e7, adis = 1.0')

% mat_cen = inv(JMcene7full/1e7)*Jcenfull;
% [V,E_mat] = eig(mat_cen);
[V,E_mat] = eig(-Jcenfull,JMcene7full/1e7);
E_vec = diag(E_mat);
figure
plot(real(E_vec),imag(E_vec),'marker','x','LineStyle','none')
title('dtinv = 1e7, adis = 0.0')

