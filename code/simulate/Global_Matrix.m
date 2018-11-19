function [KB,MB,cor1,cor2]=Global_Matrix(NumberElements,number_vehicles,kb,mb,ele1,ele2)
KB=zeros(2*(NumberElements+1),2*(NumberElements+1));
MB=zeros(2*(NumberElements+1),2*(NumberElements+1));
kk=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
mm=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
cor1=zeros(NumberElements,4);
cor2=zeros(NumberElements,4);
for i=1:NumberElements
    if number_vehicles==1
        %Monitoring Vehicle
   cor1(i,:)=ele1(i,2:5);
cor2=zeros(NumberElements,4);
   kk(:,:,i)=KInsert(kb,cor1(i,:),2*(NumberElements+1));
   KB=KB+kk(:,:,i); % Beam stiffness matrix
   mm(:,:,i)=KInsert(mb,cor1(i,:),2*(NumberElements+1));
   MB=MB+mm(:,:,i); % Beam mass matrix
    else
   cor1(i,:)=ele1(i,2:5);
   cor2(i,:)=ele2(i,2:5);
   kk(:,:,i)=KInsert(kb,cor1(i,:),2*(NumberElements+1));
   KB=KB+kk(:,:,i); % Beam stiffness matrix
   mm(:,:,i)=KInsert(mb,cor1(i,:),2*(NumberElements+1));
   MB=MB+mm(:,:,i); % Beam mass matrix
    end
end
end