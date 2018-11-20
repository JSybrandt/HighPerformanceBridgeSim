temperature(1,end)=40;
for Day=1:50
for ii=1:n
   if hour==23
TopHourTemp=(temperature(Day+1,2)-45)*5/9;
BottomHourTemp=(temperature(Day,hour+1)-45)*5/9;
Tact(Day,ii)=BottomHourTemp+(Td(ii)-hour)*(TopHourTemp-BottomHourTemp);
    else
TopHourTemp=(temperature(Day,hour+2)-45)*5/9;
BottomHourTemp=(temperature(Day,hour+1)-45)*5/9;
Tact(Day,ii)=BottomHourTemp+(Td(ii)-hour)*(TopHourTemp-BottomHourTemp);
    end
% Variables for modulus modification factor
Q=normrnd(1.0129,.003);
S=normrnd(-.0048,.0001);
R=normrnd(.1977,.0027);
tu=normrnd(3.1466,.0861);
lam=normrnd(-1.1012,.0513);
% Modified Modulus
u0=Q+S*Tact(Day,ii)+R*(1-erf((Tact(Day,ii)-lam)/tu)); % Modification factor
E0=u0*E;

% Update the hour
if Td(ii)>=(hour+1)
    hour=hour+1;
end

% Healthy elemental matrices
kb=KBeam(E0,I,l); % Stiffness matrix for bridge
mb=MBeam(mu(Day),l); % Consistent mass matrix for bridge

% Global Beam Matricies
KB=zeros(2*(NumberElements+1),2*(NumberElements+1));
MB=zeros(2*(NumberElements+1),2*(NumberElements+1));
cor_Mon=zeros(NumberElements,4);
kk=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
mm=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
for i=1:NumberElements
   cor_Mon(i,:)=ele(i,2:5); 
   kk(:,:,i)=KInsert(kb,cor_Mon(i,:),2*(NumberElements+1)); 
 KB=KB+kk(:,:,i); % Beam stiffness matrix
 mm(:,:,i)=KInsert(mb,cor_Mon(i,:),2*(NumberElements+1));
 MB=MB+mm(:,:,i); % Beam mass matrix
end

[M,K]=boundarycondition(MB,KB,NumberElements);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s] 
wn_FEA(Day,ii)=ef(1,1)/(2*pi); % sorted natural angular frequencies




end
figure(1)
set(gcf,'color','white')
scatter(Tact(Day,:),wn_FEA(Day,:),'MarkerEdgeColor',[0 0 0]); hold on
xlabel('Temperature ({\circ}C)')
ylabel('Frequency (Hz)')
plotformat
box('on')
end

