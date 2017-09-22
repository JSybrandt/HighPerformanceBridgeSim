function exitCode = SimulateDay(inPath, Day, outPath)
load(inPath);

for ii=1:n % This is the loop for changing vehicle variables over time

% Temperature bounds
if Temp==1;
    if hour==23
TopHourTemp=(temperature(Day+1,2)-32)*5/9;
BottomHourTemp=(temperature(Day,hour+1)-32)*5/9;
Tact(Day,ii)=BottomHourTemp+(Td(ii)-hour)*(TopHourTemp-BottomHourTemp);
    else
TopHourTemp=(temperature(Day,hour+2)-32)*5/9;
BottomHourTemp=(temperature(Day,hour+1)-32)*5/9;
Tact(Day,ii)=BottomHourTemp+(Td(ii)-hour)*(TopHourTemp-BottomHourTemp);
    end
% Variables for modulus modification factor
Q=1.0099+rand(1)*(1.0159-1.0099);
S=-.0049+rand(1)*(-.0047+.0049);
R=.195+rand(1)*(.2004-.195);
tu=3.0605+rand(1)*(3.2327-3.0605);
lam=-1.1525+rand(1)*(-1.0499+-1.1525);
% Modified Modulus
u0=Q+S*Tact(Day,ii)+R*(1-erf((Tact(Day,ii)-lam)/tu)); % Modification factor
E0=u0*E;
else
u0=1;
E0=E;
end

% Wind loads
if Wind==1
ce=1.3; % Exposure factor (1.3 for normal structures)
cfz=-.9; % Lift coefficient (constant value)
CW=ce*cfz; % Wind load factor
gas_constant=287.05; %J/(kg*k)
air_pressure=101325; %Pa (Assuming constant air pressure through out the year)
air_density=air_pressure/((Tact(Day,ii)+273)*gas_constant); %kg/m^3
    if hour==23
TopHourWind=wind(Day+1,2)*0.44704;
BottomHourWind=wind(Day,hour+1)*0.44704;
WindVelocity(Day,ii)=BottomHourWind+(Td(ii)-hour)*(TopHourWind-BottomHourWind);
ForceWind(Day,ii)=.5*air_density*WindVelocity(Day,ii)^2*CW*W*l; %(This load will be applied at every node. It is assumed the wind load is equally distribuuted across the surface of the stucture)
    else
TopHourWind=wind(Day,hour+2)*0.44704;
BottomHourWind=wind(Day,hour+1)*0.44704;
WindVelocity(Day,ii)=BottomHourWind+(Td(ii)-hour)*(TopHourWind-BottomHourWind);
ForceWind(Day,ii)=.5*air_density*WindVelocity(Day,ii)^2*CW*W*l; %(This load will be applied at every node. It is assumed the wind load is equally distribuuted across the surface of the stucture)
    end
if Td(ii)>=(hour+1)
    hour=hour+1;
end
end

% Applying damaged modulus to elemental matrix
       if Day>=DayDamage1; 
           if sum(ED1(1,1:(Day-DayDamage1+1)))/E >= .4
               E_damaged1=E*.6;
           else
E_damaged1=(E-sum(ED1(1,1:(Day-DayDamage1+1)))); % Overall Damaged Modulus 1
           end
           kb_damaged1=KBeam(E_damaged1*u0,I,l); 
       end

       if Day>=DayDamage2;
           if sum(ED2(1,1:(Day-DayDamage2+1)))/E >= .4
               E_damaged2=E*.6;
           else
E_damaged2=(E-sum(ED2(1,1:(Day-DayDamage2+1)))); % Overall Damaged Modulus 1
           end
           kb_damaged2=KBeam(E_damaged2*u0,I,l);
       end

% Healthy elemental matrices
kb=KBeam(E0,I,l); % Stiffness matrix for bridge
mb=MBeam(mu(Day),l); % Consistent mass matrix for bridge

% Vehicle parameters
mv=VehicleVariables(row(Day,ii),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 3 vehicles)
if row(Day,ii)==1;
mw=VehicleVariables(row(Day,ii),2)+round(-.2+randi(1,1)*(2+2)); % wheel mass of vehicle kg
else
mw=VehicleVariables(row(Day,ii),2)+randi([-3 3],1,1);    
end
kv=VehicleVariables(row(Day,ii),3); %Stiffness of vehicle spring N/m
kw=VehicleVariables(row(Day,ii),6); %Stiffness of vehicle tire N/m
cs=VehicleVariables(row(Day,ii),4); % Damping of vehicle spring N*s/m
cw=VehicleVariables(row(Day,ii),5); % Damping of vehicle tire N*s/m
fv(Day,ii)=sqrt(kv/mv)/(2*pi); % Natural frequency of sprung mass

% 
VehicleMass(Day,ii)=VehicleVariables(row(Day,ii),1);
WheelMass(Day,ii)=VehicleVariables(row(Day,ii),2);
SuspensionStiffness(Day,ii)=VehicleVariables(row(Day,ii),3);
WheelStiffness(Day,ii)=VehicleVariables(row(Day,ii),6);
SuspensionDamping(Day,ii)=VehicleVariables(row(Day,ii),4);
WheelDamping(Day,ii)=VehicleVariables(row(Day,ii),5);

% Time and position
dT=.001; % Time Step
T=0:dT:L/V(Day,ii); % Total time to cross bridge
Kt=length(T);
xg=0; % Initial global position
j=1; % Initial row for elemental matrix
J=2; % Initial row for nodal matrix


%% Matricies and Initial Conditions
% Beam Matricies
KB=zeros(2*(NumberElements+1),2*(NumberElements+1));
MB=zeros(2*(NumberElements+1),2*(NumberElements+1));
cor=zeros(NumberElements,4);
kk=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
mm=zeros(2*(NumberElements+1),2*(NumberElements+1),NumberElements);
for i=1:NumberElements
   cor(i,:)=ele(i,2:5); 
   kk(:,:,i)=KInsert(kb,cor(i,:),2*(NumberElements+1));
   if Damage==1
   if i==DamageLocation(1,1); % Insert damage state 1
       if Day>=DayDamage1;       
   kk(:,:,i)=KInsert(kb_damaged1,cor(i,:),2*(NumberElements+1));
       end
   end
   
   if i==DamageLocation(2,1); % Insert damage state 2
       if Day>=DayDamage2;
   kk(:,:,i)=KInsert(kb_damaged2,cor(i,:),2*(NumberElements+1));
       end
   end
   end
 KB=KB+kk(:,:,i); % Beam stiffness matrix
 mm(:,:,i)=KInsert(mb,cor(i,:),2*(NumberElements+1));
 MB=MB+mm(:,:,i); % Beam mass matrix
end

% Load Matricies
pc=zeros(2*(NumberElements+1),1);
qc=zeros(2*(NumberElements+1),1);

% Vehicle matricies
muu=mv; muw=0;  mwu=0;  mww=mw;
kuu=kv; kuw=-kv;    kwu=-kv;    kww=kv+kw;
cuu=cs; cuw=-cs;    cwu=-cs;    cww=cs+cw;
fue_t_dt=0; fwe_t_dt=-9.81*mv-9.81*mw;
lu=0;
lw=1;

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt+1);
AB1=zeros(2*(NumberElements+1),Kt+1);
VB1=zeros(2*(NumberElements+1),Kt+1);

% Initial Condition Vehicle
zu(:,1)=[0;0]; % Initial displacement vehicle
zv(:,1)=[0;0]; % Initial velocity of vehicle
za(:,1)=[0;0]; % Initial acceleration of vehicle

%% Newmark Method
%Integration Parameters
gamma=.5;
phi=.25;
% Integration Constants
a0 = 1/(phi*dT^2);     a4 = gamma*dT;
a1 = 1/(phi*dT);       a5 = gamma/(phi*dT);
a2 = 1/(2*phi)-1;      a6 =  (gamma/phi)-1;
a3 = dT*(1-gamma);     a7 =  (dT/2)*((gamma/phi)-2);

% Contact Matricies
PSIuu=a0*muu+a5*cuu+kuu;
PSIwu=a0*mwu+a5*cwu+kwu;
mc=lw\(mww-PSIwu*PSIuu\muw); %Mass contact matrix
cc=lw\(cww-PSIwu*PSIuu\cuw); %Damping contact matrix
kc=lw\(kww-PSIwu*PSIuu\kuw); %Stiffness contact matrix

%% Pre-allocate matrix sizes
dz=zeros(1,Kt-1); % Change in upper vehicle displacement
quc_tdt=zeros(1,Kt);% Uppdated contact Force
ac=zeros(1,Kt-1);  % Contact acceleration
dc=zeros(1,Kt-1);  % Contact displacement
vc=zeros(1,Kt-1);  % Contact velocity
ab=zeros(4,Kt-1);  % Local element acceleration
db=zeros(4,Kt-1);  % Local element displacement
vb=zeros(4,Kt-1);  % Local element velocity
du=zeros(NumberElements*2,Kt-1); % Change in displacement of global matrices
Nc=zeros(Kt-1,4);   % Shape function
Ncd=zeros(Kt-1,4);  % Derivative of shape function
Ncdd=zeros(Kt-1,4); % Double derivative of shape function
M=zeros(2*NumberElements,2*NumberElements);
K=zeros(2*NumberElements,2*NumberElements);
C=zeros(2*NumberElements,2*NumberElements);
wn_FEA=zeros(2*NumberElements,Kt-1);


%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
xe=nodes(J,2); % elment end location
cor(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1; % This loop is for calculating the accleration response for each vehicle crossing
xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
t=xc/V(Day,ii); % Local time
Nc(i,:)=[1-3*xb^2+2*xb^3, xc*(1-2*xb+xb^2), 3*xb^2-2*xb^3, xc*(xb^2-xb)]; % Shape function Row Vector
Ncd(i,:)=[-6*V(Day,ii)^2*t/l^2+6*V(Day,ii)^3*t^2/l^3, V(Day,ii)-4*V(Day,ii)^2*t/l+3*V(Day,ii)^3*t^2/l^2, 6*V(Day,ii)^2*t/l^2-6*V(Day,ii)^3*t^2/l^3, 3*V(Day,ii)^3*t^2/l^2-2*V(Day,ii)^2*t/l];
Ncdd(i,:)=[-6*V(Day,ii)^2/l^2+12*V(Day,ii)^3*t/l^3, -4*V(Day,ii)^2/l+6*V(Day,ii)^3*t/l^2, 6*V(Day,ii)^2/l^2-12*V(Day,ii)^3*t/l^3, 6*V(Day,ii)^3*t/l^2-2*V(Day,ii)^2/l];
Nct=transpose(Nc); % Column Vector

% Load Vectors    
qu_t=muu*(a1*zv(1,i)+a2*za(1,i))+cuu*(a6*zv(1,i)+a7*za(1,i))-kuu*zu(1,i);
qw_t=mwu*(a1*zv(1,i)+a2*za(1,i))+cwu*(a6*zv(1,i)+a7*za(1,i))-kwu*zu(1,i);
pc_tdt=lw\(PSIwu*PSIuu\fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu*PSIuu\qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct(:,i)*mc*Nc(i,:); %Mass
cc_st=Nct(:,i)*(2*V(Day,ii)*mc*Ncd(i,:)+cc*Nc(i,:)); %Damping
kc_st=Nct(:,i)*(V(Day,ii)^2*mc*Ncdd(i,:)+V(Day,ii)*cc*Ncd(i,:)+kc*Nc(i,:)); %Stiffness

% Equivalent nodal loads
pc_tdt_st=Nct(:,i)*pc_tdt;
qc_t_st=Nct(:,i)*qc_t;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh=KInsert(kc_st,cor(j,:),2*(NumberElements+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor(j,:),2*(NumberElements+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

% Apply boundary conditions to calculate damping matirx 
[M,K]=boundarycondition(MB1,KB1,M,K,NumberElements);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s] 
wn_FEA(:,i)=ef/(2*pi); % sorted natural angular frequencies [Hz]

% Beam damping matrix
al=2*beta*wn_FEA(1,i)*wn_FEA(2,i)/(wn_FEA(1,i)+wn_FEA(2,i)); % Alpha for Rayleigh Damping
be=2*beta/(wn_FEA(1,i)+wn_FEA(2,i)); % Beta for Rayleigh Damping
CB=al*MB1+be*KB1; % Damping Matrix for beam

% Global Damping
cch=KInsert(cc_st,cor(j,:),2*(NumberElements+1));
CB=CB+cch; % Updated beam damping matrix
% Apply boundary condition for damping matrix
C(1:2*NumberElements-1,1:2*NumberElements-1)=CB(2:2*NumberElements,2:2*NumberElements);
C(2*NumberElements,1:2*NumberElements-1)=CB(2*(NumberElements+1),2:2*NumberElements);
C(1:2*NumberElements-1,2*NumberElements)=CB(2:2*NumberElements,2*(NumberElements+1));
C(2*NumberElements,2*NumberElements)=CB(2*(NumberElements+1),2*(NumberElements+1));

% Global Contact Loads
PC1=pc; 
pch=KInsert2(pc_tdt_st,cor(j,:),2*(NumberElements+1));
PC1=PC1+pch;

QC1=qc;
qch=KInsert2(qc_t_st,cor(j,:),2*(NumberElements+1));
QC1=QC1+qch;

% Apply Boundary Conditions for Load Matrices
% Delete first and next to last row
PC=zeros(2*NumberElements,1);
QC=zeros(2*NumberElements,1);
[PC,QC]=bcloads(NumberElements,PC1,PC,QC1,QC);

if Wind==1
    FW(:)=ForceWind(Day,ii);
    % Apply Numark Method
    RS=FW-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
    LS=(a0*M+a4*a0*C+K);
else 
    % Apply Numark Method
    RS=-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
    LS=(a0*M+a4*a0*C+K);
end

% Finding unknown displacements
du(:,i)=LS\RS;

% Future Disp, Vel, Acc
UB(:,i+1)=UB(:,i)+du(:,i); % Future displacement
AB(:,i+1)=a0*du(:,i)-a1*VB(:,i)-a2*AB(:,i); % Future Acceleration
VB(:,i+1)=VB(:,i)+a3*AB(:,i)+a4*AB(:,i+1); % Future Velocity

% Adding noise to the acceleration of bridge
AB(:,i+1)=AB(:,i+1)+(AB(:,i+1)*(-.025))+randn(2*NumberElements,1).*(AB(:,i+1)*(.025)-AB(:,i+1)*(-.025)); % Future Acceleration

% Add back constrained dof
UB1(2:2*NumberElements,i)=UB(1:2*NumberElements-1,i);
UB1(2*(NumberElements+1),i)=UB(2*NumberElements,i);
AB1(2:2*NumberElements,i)=AB(1:2*NumberElements-1,i);
AB1(2*(NumberElements+1),i)=AB(2*NumberElements,i);
VB1(2:2*NumberElements,i)=VB(1:2*NumberElements-1,i);
VB1(2*(NumberElements+1),i)=VB(2*NumberElements,i);

% Element displacement, velocity, acceleration
db(:,i)=UB1(cor(j,:),i); % Local Displacement
vb(:,i)=VB1(cor(j,:),i); % Local Velocity
ab(:,i)=AB1(cor(j,:),i); % Local Acceleration

% Contact points
dc(:,i)=Nc(i,:)*db(:,i); % Contact displacement
vc(:,i)=V(Day,ii)*Ncd(i,:)*db(:,i)+Nc(i,:)*vb(:,i); % Contact velocity
ac(:,i)=Nc(i,:)*ab(:,i)+2*V(Day,ii)*Ncd(i,:)*vb(:,i)+(V(Day,ii)^2)*Ncdd(i,:)*db(:,i); % Contact acceleration

% Future displacement, velocity and acceleration in lower vehicle
zu(2,i+1)=dc(:,i); % Future displacement
za(2,i+1)=vc(:,i); % Future Acceleration
zv(2,i+1)=ac(:,i); % Future Velocity

% Contact Force
% Vi_Tdt(:,i+1)=pc_tdt+qc_t+(mc*ac(:,i)+cc*vc(:,i)+kc*dc(:,i));
quc_tdt(:,i+1)=muw*ac(:,i)+cuw*vc(:,i)+kuw*dc(:,i);

% Change in upper vehicle displacement
dz(:,i)=PSIuu\(0-quc_tdt(:,i+1)+qu_t);

% Future displacement, vnocity and accneration in upper vehicle
zu(1,i+1)=zu(1,i)+dz(:,i); % Future displacement
za(1,i+1)=a0*dz(:,i)-a1*zv(1,i)-a2*za(1,i); % Future Acceleration
zv(1,i+1)=zv(1,i)+a3*za(1,i)+a4*za(1,i+1); % Future Velocity

if xc>=l
j=j+1;
J=J+1;
xo=nodes(J-1,2); % Updated for new element
xe=nodes(J,2); % Updated for new element
cor(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V(Day,ii);
end
AllFrequencyData{Day,ii}=wn_FEA(1:2,:);
AccelerationVehicle{Day,ii}=za;
Time{Day,ii}=T;
end
save(outPath)
exitCode = 0;
end
