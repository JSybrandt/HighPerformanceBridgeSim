function [PC,QC,RC]=bcloads(NumberElements,PC,QC,RC)
if nargin < 4
RC=zeros(2*(NumberElements+1),1);
end
PC(2*NumberElements+1,:)=[];
PC(1,:) = [];
QC(2*NumberElements+1,:)=[];
QC(1,:) = [];
RC(2*NumberElements+1,:)=[];
RC(1,:) = [];
end