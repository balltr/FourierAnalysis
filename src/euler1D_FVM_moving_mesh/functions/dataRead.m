function [data] = dataRead(path,prefix)
filename = strcat(path,prefix,'*');
files = dir(filename);
numfiles = length(files);
data = cell(1, numfiles);
for k = 1:numfiles
  filename = sprintf('%s%d.txt', prefix, k-1);
  data{k} = importdata(filename);
end
end