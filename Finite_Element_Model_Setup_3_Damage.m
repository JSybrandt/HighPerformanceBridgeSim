clc; clear all;
% Pre-Processing (Change values manually)
Surface=0; % 1 if surface is considered, 0 otherwise
Wind=1; % 1 if wind effects are considered, 0 otherwise
Temp=1; % 1 if temp effects are considered, 0 otherwise
RainEffects=1; % 1 if rain effects are considered, 0 otherwise
Damage=1; % 1 if damage effects are considered, 0 otherwise
lim=362; % Number of days monitoring a single bridge
NumberBridges=1; %How many bridges you want to monitor 
NumberElements=10; % Number of elements bridge is divided into

% Load variable arrays
BridgeVariables=load('BridgeVariables.dat');
VehicleVariables=load('VehicleData.dat');
load('USW00003870.hourly.mat'); 


%% Pre Analysis Calculations (Sets up environmental and vehicle parameters)
Td=0:.016667:24; %Time of day record was taken (Assumes vehicles cross bridge every min of a day)
n=length(Td); % Number of cases looking at a single car crossing bridge
Day=1; % Number of days spent monitoring a single bridge (Starts at 1)

% Bridge Parameters
Bridge=1; % Indicates which bridge is being tested
L=BridgeVariables(Bridge,1); % Length m
W=BridgeVariables(Bridge,6);% Width m
l=L/NumberElements; % Length of each individual element
CN=(NumberElements+2); % Central node of bridge
E=BridgeVariables(Bridge,3); % Modulus of elasticity of bridge N/m^2
I=BridgeVariables(Bridge,4); % Moment of Inertia m^4
mu=BridgeVariables(Bridge,2); % mass per unit length kg/m
wn_ref=sqrt((pi^4/L^4)*(E*I/mu)); % First natural frequency of bridge without temp effects
beta=BridgeVariables(Bridge,5); % Damping of bridge
% mb=(mu*l/2)*[1 0 0 0;0 1/12*(l)^2 0 0;0 0 1 0; 0 0 0 1/12*(l)^2]; % Lumped mass matrix
[ele, nodes]=element(NumberElements,l,L);
% Pangle=2*pi*rand(size(10*L)); % Phase angle to be used for bridge surface profile

% Damage Variables
if Damage==1;
DamageLocation=randsample(ele(:,1),2); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
DayDamage2=round(lim*.5)+round(rand(1)*(lim*.67-lim*.5)); % The day second damage is iniciated on bridge
ED1=[((.1+rand(1)*(.2-.1))*E),.005*rand(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1
ED2=[((.1+rand(1)*(.2-.1))*E),.005*rand(1,(lim-DayDamage2+1))*E]; % Damaged Modulus 2
else
    DayDamage1=1000000;
    DayDamage2=1000000;
end

% Environment matrices (Rain, Temp, Wind)
if RainEffects==1;
Rain=[0,0,randi([0 5],1,NumberBridges*lim+1)];
Rain(Rain<4)=0;
Rain(Rain>0)=1;
end
if Temp==1
    Tact=zeros(lim,n); %Actual temperature degrees Celsius
    hour=0;
end
if Wind==1
    WindVelocity=zeros(NumberBridges*lim,n);
    ForceWind=zeros(NumberBridges*lim,n);
    FW=zeros(2*NumberElements,1);
end

% Randomized vehicle selection and speed
row=randi([1 3],NumberBridges*lim,n); % Randomly selects which vehicle is crossing bridge
V=randi([10 25],NumberBridges*lim,n); % Speed of vehicle m/s

% Storage matrices and cell arrays (used to store information for machine
% learning)
Time=cell(lim*NumberBridges,n);
VehicleMass=zeros(NumberBridges*lim,n);
WheelMass=zeros(NumberBridges*lim,n);
SuspensionStiffness=zeros(NumberBridges*lim,n);
WheelStiffness=zeros(NumberBridges*lim,n);
SuspensionDamping=zeros(NumberBridges*lim,n);
WheelDamping=zeros(NumberBridges*lim,n);
AccelerationVehicle=cell(lim*NumberBridges,n);
AllFrequencyData=cell(lim*NumberBridges,n);
fv=zeros(lim*NumberBridges,n);



%% Analysis
% (Parallel Loop)

for jjj=1:(NumberBridges*lim); % This is the loop for the number of days

if Day>lim;
Day=1;
Bridge=Bridge+1; % Start monitoring another bridge
L=BridgeVariables(Bridge,1); % Length m
W=BridgeVariables(Bridge,6);% Width m
l=L/NumberElements; % Length of each individual element
CN=(NumberElements+2); % Central node of bridge
E=BridgeVariables(Bridge,3); % Modulus of elasticity of bridge N/m^2
I=BridgeVariables(Bridge,4); % Moment of Inertia m^4
mu=BridgeVariables(Bridge,2); % mass per unit length kg/m
beta=BridgeVariables(Bridge,5); % Damping of bridge
wn_ref=sqrt((pi^4/L^4)*(E*I/mu)); % First natural frequency of bridge without temp effects
% mb=(mu*l/2)*[1 0 0 0;0 1/12*(l)^2 0 0;0 0 1 0; 0 0 0 1/12*(l)^2]; % Lumped mass matrix
[ele, nodes]=element(NumberElements,l,L);
% Pangle=2*pi*rand(size(10*L)); % Phase angle to be used for bridge surface profile

if Damage==1;
DamageLocation=randsample(ele(:,1),2); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
DayDamage2=round(lim*.5)+round(rand(1)*(lim*.67-lim*.5)); % The day second damage is iniciated on bridge
ED1=[((.1+rand(1)*(.2-.1))*E),.005*rand(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1
ED2=[((.1+rand(1)*(.2-.1))*E),.005*rand(1,(lim-DayDamage2+1))*E]; % Damaged Modulus 2
else
    DayDamage1=1000000;
    DayDamage2=1000000;
end
end

% Rain bounds
if RainEffects==1;
if jjj>=3
if Rain(jjj)==1 && mu>=1.01*BridgeVariables(Bridge,2);
       mu=1.01*BridgeVariables(Bridge,2);
else if Rain(jjj)==1
       mu=mu+.001*rand(1,1)*BridgeVariables(Bridge,2);
    else if Rain(jjj)==0 && Rain(jjj-1)==1
   mu=mu-.001*rand(1,1)*BridgeVariables(Bridge,2);
        else if Rain(jjj)==0 && mu>BridgeVariables(Bridge,2)
           mu=mu-.001*rand(1,1)*BridgeVariables(Bridge,2);
            else if Rain(jjj)==0 && mu<=.99*BridgeVariables(Bridge,2)
                    mu=.99*BridgeVariables(Bridge,2);
                end
            end
        end
    end
end
end
end
 
for ii=1:n % This is the loop for changing vehicle variables over time

% Temperature bounds
if Temp==1;
TopHourTemp=(temperature(jjj,hour+2)-32)*5/9;
BottomHourTemp=(temperature(jjj,hour+1)-32)*5/9;
Tact(jjj,ii)=BottomHourTemp+(Td(ii)-hour)*(TopHourTemp-BottomHourTemp);
if Td(ii+1)>=(hour+1)
    hour=hour+1;
end
% Variables for modulus modification factor
Q=1.0099+rand(1)*(1.0159-1.0099);
S=-.0049+rand(1)*(-.0047+.0049);
R=.195+rand(1)*(.2004-.195);
tu=3.0605+rand(1)*(3.2327-3.0605);
lam=-1.1525+rand(1)*(-1.0499+-1.1525);
% Modified Modulus
u0=Q+S*Tact(jjj,ii)+R*(1-erf((Tact(jjj,ii)-lam)/tu)); % Modification factor
E0=u0*E;
else
u0=1;
E0=E;
end

% wn=sqrt((pi^4/L^4)*(E0*I/mu)); % First natural frequency modified by temperature
% f1=wn/(2*pi); % First natural frequency in Hz
% w2=sqrt((2^4*pi^4/L^4)*(E0*I/mu)); % Second natural frequency

% Wind loads
if Wind==1
ce=1.3; % Exposure factor (1.3 for normal structures)
cfz=-.9; % Lift coefficient (constant value)
CW=ce*cfz; % Wind load factor
gas_constant=287.05; %J/(kg*k)
air_pressure=101325; %Pa (Assuming constant air pressure through out the year)
air_density=air_pressure/((Tact(jjj,ii)+273)*gas_constant); %kg/m^3
TopHourWind=wind(jjj,hour+2)*0.44704;
BottomHourWind=wind(jjj,hour+1)*0.44704;
WindVelocity(jjj,ii)=BottomHourWind+(Td(ii)-hour)*(TopHourWind-BottomHourWind);
ForceWind(jjj,ii)=.5*air_density*WindVelocity(jjj,ii)^2*CW*W*l; %(This load will be applied at every node. It is assumed the wind load is equally distribuuted across the surface of the stucture)
% (add a guest factor that can turn on and off randomly, will only do this
% for higher than normal wind speeds)  
end

% Applying damaged modulus to elemental matrix
       if Day>=DayDamage1; 
           if sum(ED1(1,1:(jjj-DayDamage1+1)))/E >= .4
               E_damaged1=E*.6;
           else
E_damaged1=(E-sum(ED1(1,1:(jjj-DayDamage1+1)))); % Overall Damaged Modulus 1
           end
           kb_damaged1=KBeam(E_damaged1*u0,I,l); 
       end

       if Day>=DayDamage2;
           if sum(ED2(1,1:(jjj-DayDamage2+1)))/E >= .4
               E_damaged2=E*.6;
           else
E_damaged2=(E-sum(ED2(1,1:(jjj-DayDamage2+1)))); % Overall Damaged Modulus 1
           end
           kb_damaged2=KBeam(E_damaged2*u0,I,l);
       end

% Healthy elemental matrices
kb=KBeam(E0,I,l); % Stiffness matrix for bridge
mb=MBeam(mu,l); % Consistent mass matrix for bridge

% Vehicle parameters
mv=VehicleVariables(row(jjj,ii),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 3 vehicles)
if row(jjj,ii)==1;
mw=VehicleVariables(row(jjj,ii),2)+round(-.2+randi(1,1)*(2+2)); % wheel mass of vehicle kg
else
mw=VehicleVariables(row(jjj,ii),2)+randi([-3 3],1,1);    
end
kv=VehicleVariables(row(jjj,ii),3); %Stiffness of vehicle spring N/m
kw=VehicleVariables(row(jjj,ii),6); %Stiffness of vehicle tire N/m
cs=VehicleVariables(row(jjj,ii),4); % Damping of vehicle spring N*s/m
cw=VehicleVariables(row(jjj,ii),5); % Damping of vehicle tire N*s/m
fv(jjj,ii)=sqrt(kv/mv)/(2*pi); % Natural frequency of sprung mass

% 
VehicleMass(jjj,ii)=VehicleVariables(row(jjj,ii),1);
WheelMass(jjj,ii)=VehicleVariables(row(jjj,ii),2);
SuspensionStiffness(jjj,ii)=VehicleVariables(row(jjj,ii),3);
WheelStiffness(jjj,ii)=VehicleVariables(row(jjj,ii),6);
SuspensionDamping(jjj,ii)=VehicleVariables(row(jjj,ii),4);
WheelDamping(jjj,ii)=VehicleVariables(row(jjj,ii),5);

% Time and position
dT=.001; % Time Step
T=0:dT:L/V(jjj,ii); % Total time to cross bridge
Kt=length(T);
xg=0; % Initial global position
j=1; % Initial row for elemental matrix
J=2; % Initial row for nodal matrix

% Road Surface Roughness
% k=3; % A-B Roughness value
% B=L/Kt ; % Sampling Interval (m)
% dN=1/L;  % Frequency Band
% n_o=0.1; % Spatial Frequency (cycles/m)
% n=dN:dN:Kt*dN % Spatial Frequency Band
% x = 0:B:L-B; % Abscissa Variable from 0 to L
% hx = zeros(size(x));
%     hx(i) = sum(sqrt(dN)*(2^k)*(1e-3)*(n_o./n)*cos(2*pi*n*x(i)+ Pangle(1:Kt)));



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
% Mv=[muu,muw;mwu,mww]; % Vehicle mass matrix
kuu=kv; kuw=-kv;    kwu=-kv;    kww=kv+kw;
% Kv=[kuu,kuw;kwu,kww]; % Vehicle stiffness matrix
cuu=cs; cuw=-cs;    cwu=-cs;    cww=cs+cw;
% Cv=[cuu,cuw;cwu,cww]; % Vehicle damping matrix
fue_t_dt=0; fwe_t_dt=-9.81*mv-9.81*mw;
% Fe=[fue_t_dt;fwe_t_dt]; % External forces
lu=0;
lw=1;
% ll=[lu;lw]; % Transformation matrix
% um=[1]; % Unit matrix

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
t=xc/V(jjj,ii); % Local time
Nc(i,:)=[1-3*xb^2+2*xb^3, xc*(1-2*xb+xb^2), 3*xb^2-2*xb^3, xc*(xb^2-xb)]; % Shape function Row Vector
Ncd(i,:)=[-6*V(jjj,ii)^2*t/l^2+6*V(jjj,ii)^3*t^2/l^3, V(jjj,ii)-4*V(jjj,ii)^2*t/l+3*V(jjj,ii)^3*t^2/l^2, 6*V(jjj,ii)^2*t/l^2-6*V(jjj,ii)^3*t^2/l^3, 3*V(jjj,ii)^3*t^2/l^2-2*V(jjj,ii)^2*t/l];
Ncdd(i,:)=[-6*V(jjj,ii)^2/l^2+12*V(jjj,ii)^3*t/l^3, -4*V(jjj,ii)^2/l+6*V(jjj,ii)^3*t/l^2, 6*V(jjj,ii)^2/l^2-12*V(jjj,ii)^3*t/l^3, 6*V(jjj,ii)^3*t/l^2-2*V(jjj,ii)^2/l];
Nct=transpose(Nc); % Column Vector

% Load Vectors    
qu_t=muu*(a1*zv(1,i)+a2*za(1,i))+cuu*(a6*zv(1,i)+a7*za(1,i))-kuu*zu(1,i);
qw_t=mwu*(a1*zv(1,i)+a2*za(1,i))+cwu*(a6*zv(1,i)+a7*za(1,i))-kwu*zu(1,i);
pc_tdt=lw\(PSIwu*PSIuu\fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu*PSIuu\qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct(:,i)*mc*Nc(i,:); %Mass
cc_st=Nct(:,i)*(2*V(jjj,ii)*mc*Ncd(i,:)+cc*Nc(i,:)); %Damping
kc_st=Nct(:,i)*(V(jjj,ii)^2*mc*Ncdd(i,:)+V(jjj,ii)*cc*Ncd(i,:)+kc*Nc(i,:)); %Stiffness

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
    FW(:)=ForceWind(jjj,ii);
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
vc(:,i)=V(jjj,ii)*Ncd(i,:)*db(:,i)+Nc(i,:)*vb(:,i); % Contact velocity
ac(:,i)=Nc(i,:)*ab(:,i)+2*V(jjj,ii)*Ncd(i,:)*vb(:,i)+(V(jjj,ii)^2)*Ncdd(i,:)*db(:,i); % Contact acceleration

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
xg=xg+dT*V(jjj,ii);
end
AllFrequencyData{jjj,ii}=wn_FEA(1:2,:);
AccelerationVehicle{jjj,ii}=za;
Time{jjj,ii}=T;
end
Day=Day+1;
end

% %% Plots
% figure(1)
% subplot(2,1,1);
% plot(V*T/L,dispB1,'color','b'); hold on
% plot(V*T/L,UB(CN,:),'color','g');
% xlabel('V*T/L ');
% ylabel('Displacement (m)');
% legend('Analytical','FEM')
% title('Displacement of Bridge')
% subplot(2,1,2);
% plot(V*T/L,dispV,'color','r'); hold on
% plot(V*T/L,zu(1,:),'color','b');
% xlabel('V*T/L ');
% ylabel('Displacement (m)');
% legend('Analytical','FEM')
% title('Displacement of Vehicle')
% 
% figure(2)
% subplot(2,1,1)
% plot(V*T/L,accB1,'color','r'); hold on
% plot(V*T/L,AB(CN,:),'color','b');
% xlabel('V*T/L ');
% ylabel('Acceleration (m/s^2)');
% legend('Analytical','FEM')
% title('Acceleration of Bridge')
% subplot(2,1,2);
% plot(V*T/L, accV1,'color','r'); hold on
% plot(V*T/L,za(1,:),'color','b');
% xlabel('V*T/L ');
% ylabel('Acceleration (m/s^2)');
% legend('Analytical','FEM')
% title('Acceleration of Vehicle')

%% Fast Fourier Transform
% Fs = 1/dT; % Sampling frequency                          
% t = (0:Kt-1)*dT; % Time vector
% f = Fs*(0:(Kt/2))/Kt; % Frequency domain
% 
% % Executing FFT for Vehicle
% fftV_FE=abs(fft(za(1,:)));
% fftV_Ana=abs(fft(accV1));
% 
% Twosided_FE = fftV_FE/Kt; % two-sided spectrum
% Twosided_Ana= fftV_Ana/Kt; 
% onesided_FE = Twosided_FE(1:Kt/2+1); % Single-sided spectrum
% onesided_Ana=Twosided_Ana(1:Kt/2+1);
% onesided_FE(2:end-1) = 2*onesided_FE(2:end-1);
% onesided_Ana(2:end-1) = 2*onesided_Ana(2:end-1);
% 
% figure(3)
% [pks_Ana,locs_Ana] = findpeaks(onesided_Ana);
% semilogx(f(1:100),onesided_Ana(1:100),f(locs_Ana),pks_Ana,'or','color','r'); hold on
% [pks_FE,locs_FE] = findpeaks(onesided_FE);
% semilogx(f(1:100),onesided_FE(1:100),f(locs_FE),pks_FE,'or','color','b');
% axis([0 100 0 .04]);
% title('Power Spectrum of Vehicle Acceleration')
% legend('Analytical','Ana Peaks','FEM','FEM Peaks')
% xlabel('f (Hz)')
% ylabel('Spectrum Amplitude')
% 
% % Executing FFt for Bridge
% fftB_FE=fft(AB(CN,:));
% fftB_Ana=fft(accB1);
% 
% Twosided_FE = abs(fftB_FE)/Kt; % two-sided spectrum
% Twosided_Ana=abs(fftB_Ana)/Kt; 
% onesided_FE = Twosided_FE(1:Kt/2+1); % Single-sided spectrum
% onesided_Ana=Twosided_Ana(1:Kt/2+1);
% onesided_FE(2:end-1) = 2*onesided_FE(2:end-1);
% onesided_Ana(2:end-1) = 2*onesided_Ana(2:end-1);
% 
% figure(4)
% semilogx(f(1:100),onesided_Ana(1:100),'color','r'); hold on
% semilogx(f(1:100),onesided_FE(1:100),'color','b');
% axis([0 100 0 .02]);
% title('Power Spectrum of Bridge Acceleration')
% legend('Analytical','FEM')
% xlabel('f (Hz)')
% ylabel('Spectrum Amplitude')


% % Obtaining Bridge Frequency
% [PeakValue_FE, maxIndex_FE] = max(onesided_FE);
% [PeakValue_Ana, maxIndex_Ana] = max(onesided_Ana);
% bridgefrequency_FE=f(maxIndex_FE);
% bridgefrequency_Ana=f(maxIndex_Ana);
% line([bridgefrequency_FE,bridgefrequency_FE],[0,PeakValue_FE],'Linestyle','--','color','k');
% 
% %Adjust Bridge Frequency
% FE_wn=bridgefrequency_FE+(w/(2*pi));
% Ana_en=bridgefrequency_Ana+(w/(2*pi));
