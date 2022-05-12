function [] = dataExport(itr,x,wd,varType)
% convert variables if needed
if varType == 1
    wd = varConvert(wd);
end

% file 1: x positions
filename = fullfile('../euler1D_FVM_moving_mesh/Results',strcat('x',num2str(itr),'.txt'));
fid = fopen(filename,'w');
fprintf(fid,'%e\n',x);
fclose(fid);

% file 2: wd1, wd2, wd3
filename = fullfile('../euler1D_FVM_moving_mesh/Results',strcat('rstrt',num2str(itr),'.txt'));
fid = fopen(filename,'w');
for i = 1:size(wd,2)
    fprintf(fid,'%.16e\t%.16e\t%.16e\t\n',wd(1,i),wd(2,i),wd(3,i));
end
fclose(fid);
end