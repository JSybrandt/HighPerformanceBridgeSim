function [ele, nodes]=element(n,l,L)
s = [1:(n+1)]; % global node indices
st = transpose(s); %transpose into rows
se = [1:n];
set = transpose(se);
se1 = [1:2:(n*2-1)];
set1 = transpose(se1);
se2 = [2:2:n*2];
set2 = transpose(se2);
se3 = [3:2:(n*2+1)];
set3 = transpose(se3);
se4 = [4:2:(n*2+2)];
set4 = transpose(se4);
x = [0:l:L];
xt = transpose(x);
nodes = [st,xt]; %nodal matrix
ele = [set,set1,set2,set3,set4]; %element matrix
end