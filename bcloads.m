function [PC,QC,RC]=bcloads(n,PC1,PC,QC1,QC,RC1,RC)
if nargin < 6
RC1=zeros(2*(n+1),1);
RC(1:2*n-1,1) = RC1(2:2*n,:);
RC(2*n,:) = RC1(2*(n+1),:);
end
PC(1:2*n-1,1) = PC1(2:2*n,:);
PC(2*n,:) = PC1(2*(n+1),:);
QC(1:2*n-1,1) = QC1(2:2*n,:);
QC(2*n,:) = QC1(2*(n+1),:);
end