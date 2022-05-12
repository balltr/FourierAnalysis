function [] = tExport(t)
filename = '../euler1D_FVM_moving_mesh/Results/time.txt';
fid = fopen(filename,'w');
fprintf(fid,'%f\n',t);
fclose(fid);
end