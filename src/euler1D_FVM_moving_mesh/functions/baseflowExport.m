function [] = baseflowExport(wdinit,varType)
% convert variables if needed
if varType == 1
    wdinit = varConvert(wdinit);
end

filename = fullfile('../euler1D_FVM_moving_mesh/Results',strcat('baseflow','.txt'));
fid = fopen(filename,'w');
fprintf(fid,'%.16e\n',wdinit);
fclose(fid);
end

