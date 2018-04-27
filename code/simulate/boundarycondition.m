function [M,K]=boundarycondition(M,K,NumberElements)
% if nargin < 6
% CB1=zeros(2*(n+1),2*(n+1));
% C=zeros(2*n,2*n);
% end  
M(2*NumberElements+1,:)=[];
M(:,2*NumberElements+1)=[];
M(1,:)=[];
M(:,1)=[];

K(2*NumberElements+1,:)=[];
K(:,2*NumberElements+1)=[];
K(1,:)=[];
K(:,1)=[];

end
