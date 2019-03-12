function [M,K,C]=boundarycondition(M,K,NumberElements,C)
if nargin < 4
C=zeros(2*(NumberElements+1),2*(NumberElements+1));
end  
M(2*NumberElements+1,:)=[];
M(:,2*NumberElements+1)=[];
M(1,:)=[];
M(:,1)=[];

K(2*NumberElements+1,:)=[];
K(:,2*NumberElements+1)=[];
K(1,:)=[];
K(:,1)=[];

C(2*NumberElements+1,:)=[];
C(:,2*NumberElements+1)=[];
C(1,:)=[];
C(:,1)=[];

end
