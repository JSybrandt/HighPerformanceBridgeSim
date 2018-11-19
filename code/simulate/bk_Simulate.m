function exitCode = SimulateDay(inPath, DaySTR, outPath)
Day = str2num(DaySTR);
load(inPath);

for ii=1:n % This is the loop for changing vehicle variables over time

% Temperature bounds
if Temp==1
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
if Damage==1 && Damage_Case==1
  if Day>=DayDamage1
    coef = sum(ED1(1,1:(Day-DayDamage1+1)))/E;
    if coef >= .5
      E_damaged1=E*.5;
      DamageClass(Day, ii) = 5;
    else
      % Overall Damaged Modulus 1
      E_damaged1=(E-sum(ED1(1,1:(Day-DayDamage1+1))));
      % Damage class between 0 and 5
      DamageClass(Day, ii) = floor(coef*10);
    end
    kb_damaged1=KBeam(E_damaged1*u0,I,l);
  else
    DamageClass(Day, ii) = 0;
  end

%        if Day>=DayDamage2
%            if sum(ED2(1,1:(Day-DayDamage2+1)))/E >= .5
%                E_damaged2=E*.5;
%            else
% E_damaged2=(E-sum(ED2(1,1:(Day-DayDamage2+1)))); % Overall Damaged Modulus 1
%            end
%            kb_damaged2=KBeam(E_damaged2*u0,I,l);
%        end
elseif Damage==1 && Damage_Case==2

    if Day>=DayDamage5
E_damaged1=E-ED1-ED2-ED3-ED4-ED5; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
      DamageClass(Day, ii) = 5;
    elseif Day>=DayDamage4
E_damaged1=E-ED1-ED2-ED3-ED4; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
      DamageClass(Day, ii) = 4;
    elseif Day>=DayDamage3
E_damaged1=E-ED1-ED2-ED3; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
      DamageClass(Day, ii) = 3;
    elseif Day>=DayDamage2
E_damaged1=E-ED1-ED2; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
      DamageClass(Day, ii) = 2;
    elseif Day>=DayDamage1
E_damaged1=E-ED1; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
      DamageClass(Day, ii) = 1;
    else
      DamageClass(Day, ii) = 0;
    end
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
   if Damage==1 && Damage_Case==1
   if i==DamageLocation % Insert damage state 1
       if Day>=DayDamage1
   kk(:,:,i)=KInsert(kb_damaged1,cor_Mon(i,:),2*(NumberElements+1));
       end
   end

%    if i==DamageLocation(2,1) % Insert damage state 2
%        if Day>=DayDamage2
%    kk(:,:,i)=KInsert(kb_damaged2,cor_Mon(i,:),2*(NumberElements+1));
%        end
%    end

   elseif Damage==1 && Damage_Case==2
       if i==DamageLocation
           if Day>=DayDamage5 || Day>=DayDamage4 || Day>=DayDamage3 || Day>=DayDamage2 ||Day>=DayDamage1
    kk(:,:,i)=KInsert(kb_damaged1,cor_Mon(i,:),2*(NumberElements+1));
           end
       end
   end

 KB=KB+kk(:,:,i); % Beam stiffness matrix
 mm(:,:,i)=KInsert(mb,cor_Mon(i,:),2*(NumberElements+1));
 MB=MB+mm(:,:,i); % Beam mass matrix
end


%% If statement for number of vehicles

if Multiple_Vehicles==1 && number_vehicles(Day,ii)==1
    % Vehicle parameters
if row(Day,ii,1)<=3
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-50 50],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2); % wheel mass of vehicle kg
elseif row(Day,ii,1)<=7
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2)+randi([-3 3],1,1);
else
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-1000 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2)+randi([-10 10],1,1);
end
kv_Mon=VehicleVariables(row(Day,ii,1),3); %Stiffness of vehicle spring N/m
kw_Mon=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cs_Mon=VehicleVariables(row(Day,ii,1),4); % Damping of vehicle spring N*s/m
cw_Mon=VehicleVariables(row(Day,ii,1),5); % Damping of vehicle tire N*s/m

K_Mon=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
M_Mon=[mv_Mon, 0; 0, mw_Mon];
ei_Mon=eig(K_Mon,M_Mon); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

% Data Storage Matrices
MonitorVehicleMass(Day,ii)=VehicleVariables(row(Day,ii,1),1);
MonitorWheelMass(Day,ii)=VehicleVariables(row(Day,ii,1),2);
MonitorSuspensionStiffness(Day,ii)=VehicleVariables(row(Day,ii,1),3);
MonitorWheelStiffness(Day,ii)=VehicleVariables(row(Day,ii,1),6);
MonitorSuspensionDamping(Day,ii)=VehicleVariables(row(Day,ii,1),4);
MonitorWheelDamping(Day,ii)=VehicleVariables(row(Day,ii,1),5);

% Time and position
dT=.001; % Time Step
T=0:dT:L/V(Day,ii,1); % Total time to cross bridge
Kt=length(T);
xg=0; % Initial global position
j=1; % Initial row for elemental matrix
J=2; % Initial row for nodal matrix

Start_Time_Following_Vehicle(Day,ii)=1000000;
Order_of_Vehicles{Day,ii}=1;
%% Matricies and Initial Conditions

% Load Matricies
pc=zeros(2*(NumberElements+1),1);
qc=zeros(2*(NumberElements+1),1);

% Vehicle matricies
muu=mv_Mon; muw=0;          mwu=0;          mww=mw_Mon;
kuu=kv_Mon; kuw=-kv_Mon;    kwu=-kv_Mon;    kww=kv_Mon+kw_Mon;
cuu=cs_Mon; cuw=-cs_Mon;    cwu=-cs_Mon;    cww=cs_Mon+cw_Mon;
fue_t_dt=0; fwe_t_dt=-9.81*mv_Mon-9.81*mw_Mon;
lw=1;

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt);
AB1=zeros(2*(NumberElements+1),Kt);
VB1=zeros(2*(NumberElements+1),Kt);

% Initial Condition Vehicle
zu=zeros(2,Kt-1); % Initial displacement Monitoring Vehicle
zv=zeros(2,Kt-1);  % Initial velocity of Monitoring Vehicle
za=zeros(2,Kt-1);  % Initial acceleration of Monitoring Vehicle

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
wn_FEA=zeros(2*NumberElements,Kt-1);
rc=zeros(2*(NumberElements+1),1);
rx=zeros(1,Kt-1);
drx=zeros(1,Kt-1);
ddrx=zeros(1,Kt-1);

%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
xe=nodes(J,2); % elment end location
cor_Mon(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1 % This loop is for calculating the accleration response for each vehicle crossing

xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
t=xc/V(Day,ii,1); % Local time
Nc(i,:)=[1-3*xb^2+2*xb^3, xc*(1-2*xb+xb^2), 3*xb^2-2*xb^3, xc*(xb^2-xb)]; % Shape function Row Vector
Ncd(i,:)=[-6*V(Day,ii,1)^2*t/l^2+6*V(Day,ii,1)^3*t^2/l^3, V(Day,ii,1)-4*V(Day,ii,1)^2*t/l+3*V(Day,ii,1)^3*t^2/l^2, 6*V(Day,ii,1)^2*t/l^2-6*V(Day,ii,1)^3*t^2/l^3, 3*V(Day,ii,1)^3*t^2/l^2-2*V(Day,ii,1)^2*t/l];
Ncdd(i,:)=[-6*V(Day,ii,1)^2/l^2+12*V(Day,ii,1)^3*t/l^3, -4*V(Day,ii,1)^2/l+6*V(Day,ii,1)^3*t/l^2, 6*V(Day,ii,1)^2/l^2-12*V(Day,ii,1)^3*t/l^3, 6*V(Day,ii,1)^3*t/l^2-2*V(Day,ii,1)^2/l];
Nct=transpose(Nc); % Column Vector

% Surface Rougness
if Surface==1
Column=find(RoadMatrix(1,:)==round(xg,3));
rx(i)=RoadMatrix(2,Column(1,1));
drx(i)=RoadMatrix(3,Column(1,1));
ddrx(i)=RoadMatrix(4,Column(1,1));
rc_tdt_st=Nct(:,i)*(V(Day,ii,1)^2*mc*ddrx(i)+V(Day,ii,1)*cc*drx(i)+kc*rx(i)); %Contact force from road irregularities
RC=rc;
rch=KInsert2(rc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
RC=RC+rch;
end

% Load Vectors
qu_t=muu*(a1*zv(1,i)+a2*za(1,i))+cuu*(a6*zv(1,i)+a7*za(1,i))-kuu*zu(1,i);
qw_t=mwu*(a1*zv(1,i)+a2*za(1,i))+cwu*(a6*zv(1,i)+a7*za(1,i))-kwu*zu(1,i);
pc_tdt=lw\(PSIwu*PSIuu\fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu*PSIuu\qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct(:,i)*mc*Nc(i,:); %Mass
cc_st=Nct(:,i)*(2*V(Day,ii,1)*mc*Ncd(i,:)+cc*Nc(i,:)); %Damping
kc_st=Nct(:,i)*(V(Day,ii,1)^2*mc*Ncdd(i,:)+V(Day,ii,1)*cc*Ncd(i,:)+kc*Nc(i,:)); %Stiffness

% Equivalent nodal loads
pc_tdt_st=Nct(:,i)*pc_tdt;
qc_t_st=Nct(:,i)*qc_t;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh=KInsert(kc_st,cor_Mon(j,:),2*(NumberElements+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor_Mon(j,:),2*(NumberElements+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

% Apply boundary conditions to calculate damping matirx
[M,K]=boundarycondition(MB1,KB1,NumberElements);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s]
wn_FEA(:,i)=ef/(2*pi); % sorted natural angular frequencies [Hz]

% Beam damping matrix
al=2*bbeta*wn_FEA(1,i)*wn_FEA(2,i)/(wn_FEA(1,i)+wn_FEA(2,i)); % Alpha for Rayleigh Damping
be=2*bbeta/(wn_FEA(1,i)+wn_FEA(2,i)); % Beta for Rayleigh Damping
C=al*MB1+be*KB1; % Damping Matrix for beam

% Global Damping
cch=KInsert(cc_st,cor_Mon(j,:),2*(NumberElements+1));
C=C+cch; % Updated beam damping matrix

% Apply boundary condition for damping matrix
C(2*NumberElements+1,:)=[];
C(:,2*NumberElements+1)=[];
C(1,:)=[];
C(:,1)=[];


% Global Contact Loads
PC=pc;
pch=KInsert2(pc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
PC=PC+pch;

QC=qc;
qch=KInsert2(qc_t_st,cor_Mon(j,:),2*(NumberElements+1));
QC=QC+qch;

% Apply Boundary Conditions for Load Matrices
% Delete first and next to last row
if Surface==1
[PC,QC,RC]=bcloads(NumberElements,PC,QC,RC);
else
[PC,QC]=bcloads(NumberElements,PC,QC);
end

if Wind==1
    FW(:)=ForceWind(Day,ii);
end

if Wind==1 && Surface==1
RS=FW-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==1 && Surface==0
RS=FW-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==0 && Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
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
db(:,i)=UB1(cor_Mon(j,:),i); % Local Displacement
vb(:,i)=VB1(cor_Mon(j,:),i); % Local Velocity
ab(:,i)=AB1(cor_Mon(j,:),i); % Local Acceleration

% Contact points
dc(:,i)=Nc(i,:)*db(:,i); % Contact displacement
vc(:,i)=V(Day,ii,1)*Ncd(i,:)*db(:,i)+Nc(i,:)*vb(:,i); % Contact velocity
ac(:,i)=Nc(i,:)*ab(:,i)+2*V(Day,ii,1)*Ncd(i,:)*vb(:,i)+(V(Day,ii,1)^2)*Ncdd(i,:)*db(:,i); % Contact acceleration

if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zu(2,i+1)=dc(:,i)+rx(i); % Future displacement
zv(2,i+1)=vc(:,i)+V(Day,ii,1)*drx(i); % Future Velocity
za(2,i+1)=ac(:,i)+V(Day,ii,1)^2*ddrx(i); % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle (without
% surface roughness)
zu(2,i+1)=dc(:,i); % Future displacement
zv(2,i+1)=vc(:,i); % Future Velocity
za(2,i+1)=ac(:,i); % Future Acceleration
end

% Contact Force
% Vi_Tdt(:,i+1)=pc_tdt+qc_t+(mc*ac(:,i)+cc*vc(:,i)+kc*dc(:,i));
quc_tdt(:,i+1)=muw*za(2,i+1)+cuw*zv(2,i+1)+kuw*zu(2,i+1);

% Change in upper vehicle displacement
dz(:,i)=PSIuu\(-quc_tdt(:,i+1)+qu_t);

% Future displacement, vnocity and accneration in upper vehicle
zu(1,i+1)=zu(1,i)+dz(:,i); % Future displacement
za(1,i+1)=a0*dz(:,i)-a1*zv(1,i)-a2*za(1,i); % Future Acceleration
zv(1,i+1)=zv(1,i)+a3*za(1,i)+a4*za(1,i+1); % Future Velocity

if xc>=l
j=j+1;
J=J+1;
xo=nodes(J-1,2); % Updated for new element
xe=nodes(J,2); % Updated for new element
cor_Mon(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V(Day,ii,1);
end % end single vehicle loop

% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
t = (0:Kt-1)*dT; % Time vector
f = Fs*(0:(Kt))/(Kt*2); % Frequency domain

% Executing FFT for Vehicle
fftV_FE=abs(fft(za(1,:),2*Kt));

Twosided_FE = fftV_FE/Kt; % two-sided spectrum
onesided_FE = Twosided_FE(1:Kt+1); % Single-sided spectrum
onesided_FE(2:end-1) = onesided_FE(2:end-1);

% % Storage matrices and cell arrays (used to store information for machine
% % learning)
Monitor_Vehicle_Time{Day,ii}=T;
Monitor_Vehicle_Acceleration{Day,ii}=za;
Monitor_Vehicle_Frequency_Amp_Data{Day,ii}=onesided_FE;
Monitor_Vehicle_Frequency_Data{Day,ii}=f;
Monitor_Vehicle_Road_Profile{Day,ii}=rx;
Monitor_Vehicle_Derivative_Road_Profile{Day,ii}=drx;
Monitor_Vehicle_Other_Derivative_Road_Profile{Day,ii}=ddrx;
%% End of first section of "If Multiple_Vehicle" loop

%% Beginning of Second section of "If Multiple_Vehicle" loop
elseif Multiple_Vehicles==1 && number_vehicles(Day,ii)==2

     % Vehicle parameters
if row(Day,ii,1)<=3
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-50 50],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2); % wheel mass of vehicle kg
elseif row(Day,ii,1)<=7
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2)+randi([-3 3],1,1);
else
mv_Mon=VehicleVariables(row(Day,ii,1),1)+randi([-1000 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Mon=VehicleVariables(row(Day,ii,1),2)+randi([-10 10],1,1);
end
kv_Mon=VehicleVariables(row(Day,ii,1),3); %Stiffness of vehicle spring N/m
kw_Mon=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cs_Mon=VehicleVariables(row(Day,ii,1),4); % Damping of vehicle spring N*s/m
cw_Mon=VehicleVariables(row(Day,ii,1),5); % Damping of vehicle tire N*s/m

K_Mon=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
M_Mon=[mv_Mon, 0; 0, mw_Mon];
ei_Mon=eig(K_Mon,M_Mon); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

if row(Day,ii,2)<=3
mv_Sec=VehicleVariables(row(Day,ii,2),1)+randi([-50 50],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Sec=VehicleVariables(row(Day,ii,2),2); % wheel mass of vehicle kg
elseif row(Day,ii,2)<=7
mv_Sec=VehicleVariables(row(Day,ii,2),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Sec=VehicleVariables(row(Day,ii,2),2)+randi([-3 3],1,1);
else
mv_Sec=VehicleVariables(row(Day,ii,2),1)+randi([-1000 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw_Sec=VehicleVariables(row(Day,ii,2),2)+randi([-10 10],1,1);
end
kv_Sec=VehicleVariables(row(Day,ii,2),3); %Stiffness of vehicle spring N/m
kw_Sec=VehicleVariables(row(Day,ii,2),6); %Stiffness of vehicle tire N/m
cs_Sec=VehicleVariables(row(Day,ii,2),4); % Damping of vehicle spring N*s/m
cw_Sec=VehicleVariables(row(Day,ii,2),5); % Damping of vehicle tire N*s/m

K_Sec=[kv_Sec, -kv_Sec; -kv_Sec, kv_Sec+kw_Sec];
M_Sec=[mv_Sec, 0; 0, mw_Sec];
ei_Sec=eig(K_Sec,M_Sec); % eigenvalues
ef_Sec=sort(real(sqrt(ei_Sec))); % sorted natural angular frequencies [rad/s]
Secondfv{Day,ii}=ef_Sec/(2*pi); % 1st and second natural frequencies of sprung mass

% Data Storage Matrices
MonitorVehicleMass(Day,ii)=VehicleVariables(row(Day,ii,1),1);
MonitorWheelMass(Day,ii)=VehicleVariables(row(Day,ii,1),2);
MonitorSuspensionStiffness(Day,ii)=VehicleVariables(row(Day,ii,1),3);
MonitorWheelStiffness(Day,ii)=VehicleVariables(row(Day,ii,1),6);
MonitorSuspensionDamping(Day,ii)=VehicleVariables(row(Day,ii,1),4);
MonitorWheelDamping(Day,ii)=VehicleVariables(row(Day,ii,1),5);

SecondVehicleMass(Day,ii)=VehicleVariables(row(Day,ii,2),1);
SecondWheelMass(Day,ii)=VehicleVariables(row(Day,ii,2),2);
SecondSuspensionStiffness(Day,ii)=VehicleVariables(row(Day,ii,2),3);
SecondWheelStiffness(Day,ii)=VehicleVariables(row(Day,ii,2),6);
SecondSuspensionDamping(Day,ii)=VehicleVariables(row(Day,ii,2),4);
SecondWheelDamping(Day,ii)=VehicleVariables(row(Day,ii,2),5);

% Time and position
dT=.001; % Time Step
T_Mon=0:dT:L/V(Day,ii,1); % Total time to cross bridge
Kt_Mon=length(T_Mon);
T_Sec=0:dT:L/V(Day,ii,2); % Total time to cross bridge
Kt_Sec=length(T_Sec);
Vehicle_order=randperm(2); % gives the order the vehicles arrive to the bridge randperm(2)

if Vehicle_order(1)==1
Start_Time_Last_Vehicle=randi([1 (Kt_Mon-200)],1,1);

    if (Start_Time_Last_Vehicle+Kt_Sec)<=Kt_Mon
Kt=Kt_Mon; % Total time for all vehicles to cross bridge
    else
Kt=Start_Time_Last_Vehicle+Kt_Sec;% Total time for all vehicles to cross bridge
    end

elseif Vehicle_order(1)==2
Start_Time_Last_Vehicle=randi([1 (Kt_Sec-200)],1,1);

    if (Start_Time_Last_Vehicle+Kt_Mon)<=Kt_Sec
Kt=Kt_Sec;% Total time for all vehicles to cross bridge
    else
Kt=Start_Time_Last_Vehicle+Kt_Mon;% Total time for all vehicles to cross bridge
    end
end
% Storage matrix to indicate what timestep the second vehicle entered the bridge
Start_Time_Following_Vehicle(Day,ii)=Start_Time_Last_Vehicle;
Order_of_Vehicles{Day,ii}=Vehicle_order;
%% Matricies and Initial Conditions

% Load Matricies
pc=zeros(2*(NumberElements+1),1);
qc=zeros(2*(NumberElements+1),1);

% Vehicle matricies
muu_Mon=mv_Mon;     muw_Mon=0;          mwu_Mon=0;          mww_Mon=mw_Mon;
kuu_Mon=kv_Mon;     kuw_Mon=-kv_Mon;    kwu_Mon=-kv_Mon;    kww_Mon=kv_Mon+kw_Mon;
cuu_Mon=cs_Mon;     cuw_Mon=-cs_Mon;    cwu_Mon=-cs_Mon;    cww_Mon=cs_Mon+cw_Mon;
fue_t_dt_Mon=0;     fwe_t_dt_Mon=-9.81*mv_Mon-9.81*mw_Mon;

muu_Sec=mv_Sec;     muw_Sec=0;          mwu_Sec=0;          mww_Sec=mw_Sec;
kuu_Sec=kv_Sec;     kuw_Sec=-kv_Sec;    kwu_Sec=-kv_Sec;    kww_Sec=kv_Sec+kw_Sec;
cuu_Sec=cs_Sec;     cuw_Sec=-cs_Sec;    cwu_Sec=-cs_Sec;    cww_Sec=cs_Sec+cw_Sec;
fue_t_dt_Sec=0;     fwe_t_dt_Sec=-9.81*mv_Sec-9.81*mw_Sec;
lw=1;

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt-1);
AB1=zeros(2*(NumberElements+1),Kt-1);
VB1=zeros(2*(NumberElements+1),Kt-1);

% Initial Condition Vehicles
zu_Mon=zeros(2,Kt-1); % Initial displacement Monitoring Vehicle
zv_Mon=zeros(2,Kt-1);  % Initial velocity of Monitoring Vehicle
za_Mon=zeros(2,Kt-1);  % Initial acceleration of Monitoring Vehicle

zu_Sec=zeros(2,Kt-1);  % Initial displacement Second Vehicle
zv_Sec=zeros(2,Kt-1);  % Initial velocity of Second Vehicle
za_Sec=zeros(2,Kt-1);  % Initial acceleration of Second Vehicle

%% Newmark Method Variables
%Integration Parameters
gamma=.5;
phi=.25;
% Integration Constants
a0 = 1/(phi*dT^2);     a4 = gamma*dT;
a1 = 1/(phi*dT);       a5 = gamma/(phi*dT);
a2 = 1/(2*phi)-1;      a6 =  (gamma/phi)-1;
a3 = dT*(1-gamma);     a7 =  (dT/2)*((gamma/phi)-2);

% Contact Matricies
PSIuu_Mon=a0*muu_Mon+a5*cuu_Mon+kuu_Mon;
PSIwu_Mon=a0*mwu_Mon+a5*cwu_Mon+kwu_Mon;
mc_Mon=lw\(mww_Mon-PSIwu_Mon*PSIuu_Mon\muw_Mon); %Monitoring Mass contact matrix
cc_Mon=lw\(cww_Mon-PSIwu_Mon*PSIuu_Mon\cuw_Mon); %Monitoring Damping contact matrix
kc_Mon=lw\(kww_Mon-PSIwu_Mon*PSIuu_Mon\kuw_Mon); %Monitoring Stiffness contact matrix

PSIuu_Sec=a0*muu_Sec+a5*cuu_Sec+kuu_Sec;
PSIwu_Sec=a0*mwu_Sec+a5*cwu_Sec+kwu_Sec;
mc_Sec=lw\(mww_Sec-PSIwu_Sec*PSIuu_Sec\muw_Sec); %Monitoring Mass contact matrix
cc_Sec=lw\(cww_Sec-PSIwu_Sec*PSIuu_Sec\cuw_Sec); %Monitoring Damping contact matrix
kc_Sec=lw\(kww_Sec-PSIwu_Sec*PSIuu_Sec\kuw_Sec); %Monitoring Stiffness contact matrix

%% Pre-allocate matrix sizes
dz_Mon=zeros(1,Kt-1); % Change in Monitoring Vehicle upper displacement
quc_tdt_Mon=zeros(1,Kt);% Uppdated Monitoring Vehicle contact force
ac_Mon=zeros(1,Kt-1);  % Monitoring Vehicle contact acceleration
dc_Mon=zeros(1,Kt-1);  % Monitoring Vehicle contact displacement
vc_Mon=zeros(1,Kt-1);  % Monitoring Vehicle contact velocity
ab_Mon=zeros(4,Kt-1);  % Monitoring Vehicle Local element acceleration
db_Mon=zeros(4,Kt-1);  % Monitoring Vehicle Local element displacement
vb_Mon=zeros(4,Kt-1);  % Monitoring Vehicle Local element velocity
Nc_Mon=zeros(Kt-1,4);   % Monitoring Vehicle Shape function
Ncd_Mon=zeros(Kt-1,4);  % Monitoring Vehicle Derivative of shape function
Ncdd_Mon=zeros(Kt-1,4); % Monitoring Vehicle Double derivative of shape function
rx_Mon=zeros(1,Kt-1);
drx_Mon=zeros(1,Kt-1);
ddrx_Mon=zeros(1,Kt-1);

dz_Sec=zeros(1,Kt-1); % Change in Second Vehicle upper displacement
quc_tdt_Sec=zeros(1,Kt);% Uppdated Second Vehicle contact force
ac_Sec=zeros(1,Kt-1);  % Second Vehicle contact acceleration
dc_Sec=zeros(1,Kt-1);  % Second Vehicle contact displacement
vc_Sec=zeros(1,Kt-1);  % Second Vehicle contact velocity
ab_Sec=zeros(4,Kt-1);  % Second Vehicle Local element acceleration
db_Sec=zeros(4,Kt-1);  % Second Vehicle Local element displacement
vb_Sec=zeros(4,Kt-1);  % Second Vehicle Local element velocity
Nc_Sec=zeros(Kt-1,4);   % Second Vehicle Shape function
Ncd_Sec=zeros(Kt-1,4);  % Second Vehicle Derivative of shape function
Ncdd_Sec=zeros(Kt-1,4); % Second Vehicle Double derivative of shape function
rx_Sec=zeros(1,Kt-1);
drx_Sec=zeros(1,Kt-1);
ddrx_Sec=zeros(1,Kt-1);

du=zeros(NumberElements*2,Kt-1); % Change in displacement of global matrices
wn_FEA=zeros(2*NumberElements,Kt-1);
rc=zeros(2*(NumberElements+1),1);

%% Calc global position and shape functions
xg_Mon=0; % Initial global position of Monitoring Vehicle
j_Mon=1; % Initial row for Monitoring Vehicle elemental matrix
J_Mon=2; % Initial row for monitoring Vehicle nodal matrix
xo_Mon=nodes(J_Mon-1,2); % Monitoring Vehicle current element start location
cor_Mon(j_Mon,:)=ele(j_Mon,2:5); % Monitoring Vehicle current element coordinate

xg_Sec=L; % Initial global position of Second Vehicle
j_Sec=1; % Initial row for Second Vehicle elemental matrix
J_Sec=2; % Initial row for Second Vehicle nodal matrix
xo_Sec=nodes2(J_Sec-1,2); % Monitoring Vehicle current element start location
cor_Sec=zeros(NumberElements,4);
cor_Sec(j_Sec,:)=ele2(j_Sec,2:5); % Monitoring Vehicle current element coordinate

for i=1:Kt-1 % This loop is for calculating the accleration response for each vehicle crossing

if i<Start_Time_Last_Vehicle

    if Vehicle_order(1)==1
xc_Mon=(xg_Mon-xo_Mon); % Local position of vehicle on bridge
xb_Mon=xc_Mon/l; % Local coordinate
t_Mon=xc_Mon/V(Day,ii,1); % Local time
Nc_Mon(i,:)=[1-3*xb_Mon^2+2*xb_Mon^3, xc_Mon*(1-2*xb_Mon+xb_Mon^2), 3*xb_Mon^2-2*xb_Mon^3, xc_Mon*(xb_Mon^2-xb_Mon)]; % Shape function Row Vector
Ncd_Mon(i,:)=[-6*V(Day,ii,1)^2*t_Mon/l^2+6*V(Day,ii,1)^3*t_Mon^2/l^3, V(Day,ii,1)-4*V(Day,ii,1)^2*t_Mon/l+3*V(Day,ii,1)^3*t_Mon^2/l^2, 6*V(Day,ii,1)^2*t_Mon/l^2-6*V(Day,ii,1)^3*t_Mon^2/l^3, 3*V(Day,ii,1)^3*t_Mon^2/l^2-2*V(Day,ii,1)^2*t_Mon/l];
Ncdd_Mon(i,:)=[-6*V(Day,ii,1)^2/l^2+12*V(Day,ii,1)^3*t_Mon/l^3, -4*V(Day,ii,1)^2/l+6*V(Day,ii,1)^3*t_Mon/l^2, 6*V(Day,ii,1)^2/l^2-12*V(Day,ii,1)^3*t_Mon/l^3, 6*V(Day,ii,1)^3*t_Mon/l^2-2*V(Day,ii,1)^2/l];
Nct_Mon=transpose(Nc_Mon); % Column Vector

xc_Sec=0; % Local position of vehicle on bridge
Nc_Sec(i,:)=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Sec(i,:)=[0, 0, 0, 0];
Ncdd_Sec(i,:)=[0, 0, 0, 0];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    else
xc_Sec=abs((xg_Sec-xo_Sec)); % Local position of vehicle on bridge
xb_Sec=xc_Sec/l; % Local coordiante
t_Sec=xc_Sec/V(Day,ii,2); % Local time
Nc_Sec(i,:)=[1-3*xb_Sec^2+2*xb_Sec^3, xc_Sec*(1-2*xb_Sec+xb_Sec^2), 3*xb_Sec^2-2*xb_Sec^3, xc_Sec*(xb_Sec^2-xb_Sec)]; % Shape function Row Vector
Ncd_Sec(i,:)=[-6*V(Day,ii,2)^2*t_Sec/l^2+6*V(Day,ii,2)^3*t_Sec^2/l^3, V(Day,ii,2)-4*V(Day,ii,2)^2*t_Sec/l+3*V(Day,ii,2)^3*t_Sec^2/l^2, 6*V(Day,ii,2)^2*t_Sec/l^2-6*V(Day,ii,2)^3*t_Sec^2/l^3, 3*V(Day,ii,2)^3*t_Sec^2/l^2-2*V(Day,ii,2)^2*t_Sec/l];
Ncdd_Sec(i,:)=[-6*V(Day,ii,2)^2/l^2+12*V(Day,ii,2)^3*t_Sec/l^3, -4*V(Day,ii,2)^2/l+6*V(Day,ii,2)^3*t_Sec/l^2, 6*V(Day,ii,2)^2/l^2-12*V(Day,ii,2)^3*t_Sec/l^3, 6*V(Day,ii,2)^3*t_Sec/l^2-2*V(Day,ii,2)^2/l];
Nct_Sec=transpose(Nc_Sec); % Column Vector

xc_Mon=0; % Local position of vehicle on bridge
Nc_Mon(i,:)=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Mon(i,:)=[0, 0, 0, 0];
Ncdd_Mon(i,:)=[0, 0, 0, 0];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    end

else

    if xg_Mon<L
xc_Mon=(xg_Mon-xo_Mon); % Local position of vehicle on bridge
xb_Mon=xc_Mon/l; % Local coordinate
t_Mon=xc_Mon/V(Day,ii,1); % Local time
Nc_Mon(i,:)=[1-3*xb_Mon^2+2*xb_Mon^3, xc_Mon*(1-2*xb_Mon+xb_Mon^2), 3*xb_Mon^2-2*xb_Mon^3, xc_Mon*(xb_Mon^2-xb_Mon)]; % Shape function Row Vector
Ncd_Mon(i,:)=[-6*V(Day,ii,1)^2*t_Mon/l^2+6*V(Day,ii,1)^3*t_Mon^2/l^3, V(Day,ii,1)-4*V(Day,ii,1)^2*t_Mon/l+3*V(Day,ii,1)^3*t_Mon^2/l^2, 6*V(Day,ii,1)^2*t_Mon/l^2-6*V(Day,ii,1)^3*t_Mon^2/l^3, 3*V(Day,ii,1)^3*t_Mon^2/l^2-2*V(Day,ii,1)^2*t_Mon/l];
Ncdd_Mon(i,:)=[-6*V(Day,ii,1)^2/l^2+12*V(Day,ii,1)^3*t_Mon/l^3, -4*V(Day,ii,1)^2/l+6*V(Day,ii,1)^3*t_Mon/l^2, 6*V(Day,ii,1)^2/l^2-12*V(Day,ii,1)^3*t_Mon/l^3, 6*V(Day,ii,1)^3*t_Mon/l^2-2*V(Day,ii,1)^2/l];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    else
c_Mon=0; % Local position of vehicle on bridge
Nc_Mon(i,:)=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Mon(i,:)=[0, 0, 0, 0];
Ncdd_Mon(i,:)=[0, 0, 0, 0];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    end

    if xg_Sec > 0
xc_Sec=abs((xg_Sec-xo_Sec)); % Local position of vehicle on bridge
xb_Sec=xc_Sec/l; % Local coordiante
t_Sec=xc_Sec/V(Day,ii,2); % Local time
Nc_Sec(i,:)=[1-3*xb_Sec^2+2*xb_Sec^3, xc_Sec*(1-2*xb_Sec+xb_Sec^2), 3*xb_Sec^2-2*xb_Sec^3, xc_Sec*(xb_Sec^2-xb_Sec)]; % Shape function Row Vector
Ncd_Sec(i,:)=[-6*V(Day,ii,2)^2*t_Sec/l^2+6*V(Day,ii,2)^3*t_Sec^2/l^3, V(Day,ii,2)-4*V(Day,ii,2)^2*t_Sec/l+3*V(Day,ii,2)^3*t_Sec^2/l^2, 6*V(Day,ii,2)^2*t_Sec/l^2-6*V(Day,ii,2)^3*t_Sec^2/l^3, 3*V(Day,ii,2)^3*t_Sec^2/l^2-2*V(Day,ii,2)^2*t_Sec/l];
Ncdd_Sec(i,:)=[-6*V(Day,ii,2)^2/l^2+12*V(Day,ii,2)^3*t_Sec/l^3, -4*V(Day,ii,2)^2/l+6*V(Day,ii,2)^3*t_Sec/l^2, 6*V(Day,ii,2)^2/l^2-12*V(Day,ii,2)^3*t_Sec/l^3, 6*V(Day,ii,2)^3*t_Sec/l^2-2*V(Day,ii,2)^2/l];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    else
xc_Sec=0; % Local position of vehicle on bridge
Nc_Sec(i,:)=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Sec(i,:)=[0, 0, 0, 0];
Ncdd_Sec(i,:)=[0, 0, 0, 0];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    end
end


if Surface==1
RC=rc; % Resets Surface Roughness matrix each loop
    if i<Start_Time_Last_Vehicle && Vehicle_order(1)==2
rx_Mon(i)=0;
drx_Mon(i)=0;
ddrx_Mon(i)=0;
    else
        if xg_Mon<=L
Column_Mon=find(RoadMatrix(1,:)==round(xg_Mon,3));
rx_Mon(i)=RoadMatrix(2,Column_Mon(1,1));
drx_Mon(i)=RoadMatrix(3,Column_Mon(1,1));
ddrx_Mon(i)=RoadMatrix(4,Column_Mon(1,1));
rc_tdt_st_Mon=Nct_Mon(:,i)*(V(Day,ii,1)^2*mc_Mon*ddrx_Mon(i)+V(Day,ii,1)*cc_Mon*drx_Mon(i)+kc_Mon*rx_Mon(i)); %Contact force from road irregularities
rch=KInsert2(rc_tdt_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
RC=RC+rch;
        else
rx_Mon(i)=RoadMatrix(2,end);
drx_Mon(i)=RoadMatrix(3,end);
ddrx_Mon(i)=RoadMatrix(4,end);
        end
    end


    if i<Start_Time_Last_Vehicle && Vehicle_order(1)==1
rx_Sec(i)=0;
drx_Sec(i)=0;
ddrx_Sec(i)=0;
    else
        if xg_Sec>=0
Column_Sec=find(RoadMatrix(1,:)==round(xg_Sec,3));
rx_Sec(i)=RoadMatrix(2,Column_Sec(1,1));
drx_Sec(i)=RoadMatrix(3,Column_Sec(1,1));
ddrx_Sec(i)=RoadMatrix(4,Column_Sec(1,1));
rc_tdt_st_Sec=Nct_Sec(:,i)*(V(Day,ii,2)^2*mc_Sec*ddrx_Sec(i)+V(Day,ii,2)*cc_Sec*drx_Sec(i)+kc_Sec*rx_Sec(i)); %Contact force from road irregularities
rch=KInsert2(rc_tdt_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
RC=RC+rch;
        else
rx_Sec(i)=RoadMatrix(2,1);
drx_Sec(i)=RoadMatrix(3,1);
ddrx_Sec(i)=RoadMatrix(4,1);
        end
    end
end

% Load Vectors
qu_t_Mon=muu_Mon*(a1*zv_Mon(1,i)+a2*za_Mon(1,i))+cuu_Mon*(a6*zv_Mon(1,i)+a7*za_Mon(1,i))-kuu_Mon*zu_Mon(1,i);
qw_t_Mon=mwu_Mon*(a1*zv_Mon(1,i)+a2*za_Mon(1,i))+cwu_Mon*(a6*zv_Mon(1,i)+a7*za_Mon(1,i))-kwu_Mon*zu_Mon(1,i);
pc_tdt_Mon=lw\(PSIwu_Mon*PSIuu_Mon\fue_t_dt_Mon-fwe_t_dt_Mon);
qc_t_Mon=lw\(PSIwu_Mon*PSIuu_Mon\qu_t_Mon-qw_t_Mon);

qu_t_Sec=muu_Sec*(a1*zv_Sec(1,i)+a2*za_Sec(1,i))+cuu_Sec*(a6*zv_Sec(1,i)+a7*za_Sec(1,i))-kuu_Sec*zu_Sec(1,i);
qw_t_Sec=mwu_Sec*(a1*zv_Sec(1,i)+a2*za_Sec(1,i))+cwu_Sec*(a6*zv_Sec(1,i)+a7*za_Sec(1,i))-kwu_Sec*zu_Sec(1,i);
pc_tdt_Sec=lw\(PSIwu_Sec*PSIuu_Sec\fue_t_dt_Sec-fwe_t_dt_Sec);
qc_t_Sec=lw\(PSIwu_Sec*PSIuu_Sec\qu_t_Sec-qw_t_Sec);

% Interaction contact matrices
mc_st_Mon=Nct_Mon(:,i)*mc_Mon*Nc_Mon(i,:); %Mass
cc_st_Mon=Nct_Mon(:,i)*(2*V(Day,ii,1)*mc_Mon*Ncd_Mon(i,:)+cc_Mon*Nc_Mon(i,:)); %Damping
kc_st_Mon=Nct_Mon(:,i)*(V(Day,ii,1)^2*mc_Mon*Ncdd_Mon(i,:)+V(Day,ii,1)*cc_Mon*Ncd_Mon(i,:)+kc_Mon*Nc_Mon(i,:)); %Stiffness

mc_st_Sec=Nct_Sec(:,i)*mc_Sec*Nc_Sec(i,:); %Mass
cc_st_Sec=Nct_Sec(:,i)*(2*V(Day,ii,2)*mc_Sec*Ncd_Sec(i,:)+cc_Sec*Nc_Sec(i,:)); %Damping
kc_st_Sec=Nct_Sec(:,i)*(V(Day,ii,2)^2*mc_Sec*Ncdd_Sec(i,:)+V(Day,ii,2)*cc_Sec*Ncd_Sec(i,:)+kc_Sec*Nc_Sec(i,:)); %Stiffness

% Equivalent nodal loads
pc_tdt_st_Mon=Nct_Mon(:,i)*pc_tdt_Mon;
qc_t_st_Mon=Nct_Mon(:,i)*qc_t_Mon;

pc_tdt_st_Sec=Nct_Sec(:,i)*pc_tdt_Sec;
qc_t_st_Sec=Nct_Sec(:,i)*qc_t_Sec;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh_Mon=KInsert(kc_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
kkh_Sec=KInsert(kc_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
KB1=KB1+kkh_Mon+kkh_Sec; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh_Mon=KInsert(mc_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
mmh_Sec=KInsert(mc_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
MB1=MB1+mmh_Mon+mmh_Sec; % Updated beam stiffness matrix

% Apply boundary conditions to calculate damping matirx
[M,K]=boundarycondition(MB1,KB1,NumberElements);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s]
wn_FEA(:,i)=ef/(2*pi); % sorted natural angular frequencies [Hz]

% Beam damping matrix
al=2*bbeta*wn_FEA(1,i)*wn_FEA(2,i)/(wn_FEA(1,i)+wn_FEA(2,i)); % Alpha for Rayleigh Damping
be=2*bbeta/(wn_FEA(1,i)+wn_FEA(2,i)); % Beta for Rayleigh Damping
C=al*MB1+be*KB1; % Damping Matrix for beam

% Global Damping
cch_Mon=KInsert(cc_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
cch_Sec=KInsert(cc_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
C=C+cch_Mon+cch_Sec; % Updated beam damping matrix

% Apply boundary condition for damping matrix
C(2*NumberElements+1,:)=[];
C(:,2*NumberElements+1)=[];
C(1,:)=[];
C(:,1)=[];

% Global Contact Loads
PC=pc;
pch_Mon=KInsert2(pc_tdt_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
pch_Sec=KInsert2(pc_tdt_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
PC=PC+pch_Mon+pch_Sec;

QC=qc;
qch_Mon=KInsert2(qc_t_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
qch_Sec=KInsert2(qc_t_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
QC=QC+qch_Mon+qch_Sec;

% Apply Boundary Conditions for Load Matrices
% Delete first and next to last row
if Surface==1
[PC,QC,RC]=bcloads(NumberElements,PC,QC,RC);
else
[PC,QC]=bcloads(NumberElements,PC,QC);
end

if Wind==1
    FW(:)=ForceWind(Day,ii);
end

if Wind==1 && Surface==1
RS=FW-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==1 && Surface==0
RS=FW-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==0 && Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
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
% AB(:,i+1)=AB(:,i+1)+(AB(:,i+1)*(-.025))+randn(2*NumberElements,1).*(AB(:,i+1)*(.025)-AB(:,i+1)*(-.025)); % Future Acceleration

% Add back constrained dof
UB1(2:2*NumberElements,i)=UB(1:2*NumberElements-1,i);
UB1(2*(NumberElements+1),i)=UB(2*NumberElements,i);
AB1(2:2*NumberElements,i)=AB(1:2*NumberElements-1,i);
AB1(2*(NumberElements+1),i)=AB(2*NumberElements,i);
VB1(2:2*NumberElements,i)=VB(1:2*NumberElements-1,i);
VB1(2*(NumberElements+1),i)=VB(2*NumberElements,i);

% Element displacement, velocity, acceleration
db_Mon(:,i)=UB1(cor_Mon(j_Mon,:),i); % Local Displacement
vb_Mon(:,i)=VB1(cor_Mon(j_Mon,:),i); % Local Velocity
ab_Mon(:,i)=AB1(cor_Mon(j_Mon,:),i); % Local Acceleration

db_Sec(:,i)=UB1(cor_Sec(j_Sec,:),i); % Local Displacement
vb_Sec(:,i)=VB1(cor_Sec(j_Sec,:),i); % Local Velocity
ab_Sec(:,i)=AB1(cor_Sec(j_Sec,:),i); % Local Acceleration

% Contact points
dc_Mon(:,i)=Nc_Mon(i,:)*db_Mon(:,i); % Contact displacement
vc_Mon(:,i)=V(Day,ii,1)*Ncd_Mon(i,:)*db_Mon(:,i)+Nc_Mon(i,:)*vb_Mon(:,i); % Contact velocity
ac_Mon(:,i)=Nc_Mon(i,:)*ab_Mon(:,i)+2*V(Day,ii,1)*Ncd_Mon(i,:)*vb_Mon(:,i)+(V(Day,ii,1)^2)*Ncdd_Mon(i,:)*db_Mon(:,i); % Contact acceleration

dc_Sec(:,i)=Nc_Sec(i,:)*db_Sec(:,i); % Contact displacement
vc_Sec(:,i)=V(Day,ii,2)*Ncd_Sec(i,:)*db_Sec(:,i)+Nc_Sec(i,:)*vb_Sec(:,i); % Contact velocity
ac_Sec(:,i)=Nc_Sec(i,:)*ab_Sec(:,i)+2*V(Day,ii,2)*Ncd_Sec(i,:)*vb_Sec(:,i)+(V(Day,ii,2)^2)*Ncdd_Sec(i,:)*db_Sec(:,i); % Contact acceleration


if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zu_Mon(2,i+1)=dc_Mon(:,i)+rx_Mon(i); % Future displacement
zv_Mon(2,i+1)=vc_Mon(:,i)+V(Day,ii,1)*drx_Mon(i); % Future Velocity
za_Mon(2,i+1)=ac_Mon(:,i)+V(Day,ii,1)^2*ddrx_Mon(i); % Future Acceleration

zu_Sec(2,i+1)=dc_Sec(:,i)+rx_Sec(i); % Future displacement
zv_Sec(2,i+1)=vc_Sec(:,i)+V(Day,ii,2)*drx_Sec(i); % Future Velocity
za_Sec(2,i+1)=ac_Sec(:,i)+V(Day,ii,2)^2*ddrx_Sec(i); % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle (without
% surface roughness)
zu_Mon(2,i+1)=dc_Mon(:,i); % Future displacement
zv_Mon(2,i+1)=vc_Mon(:,i); % Future Velocity
za_Mon(2,i+1)=ac_Mon(:,i); % Future Acceleration

zu_Sec(2,i+1)=dc_Sec(:,i); % Future displacement
zv_Sec(2,i+1)=vc_Sec(:,i); % Future Velocity
za_Sec(2,i+1)=ac_Sec(:,i); % Future Acceleration
end

% Contact Force
% Vi_Tdt(:,i+1)=pc_tdt+qc_t+(mc*ac(:,i)+cc*vc(:,i)+kc*dc(:,i));
quc_tdt_Mon(:,i+1)=muw_Mon*za_Mon(2,i+1)+cuw_Mon*zv_Mon(2,i+1)+kuw_Mon*zu_Mon(2,i+1);
quc_tdt_Sec(:,i+1)=muw_Sec*za_Sec(2,i+1)+cuw_Sec*zv_Sec(2,i+1)+kuw_Sec*zu_Sec(2,i+1);

% Change in upper vehicle displacement
dz_Mon(:,i)=PSIuu_Mon\(-quc_tdt_Mon(:,i+1)+qu_t_Mon);
dz_Sec(:,i)=PSIuu_Sec\(-quc_tdt_Sec(:,i+1)+qu_t_Sec);

% Future displacement, vnocity and accneration in upper vehicle
zu_Mon(1,i+1)=zu_Mon(1,i)+dz_Mon(:,i); % Future displacement
za_Mon(1,i+1)=a0*dz_Mon(:,i)-a1*zv_Mon(1,i)-a2*za_Mon(1,i); % Future Acceleration
zv_Mon(1,i+1)=zv_Mon(1,i)+a3*za_Mon(1,i)+a4*za_Mon(1,i+1); % Future Velocity

zu_Sec(1,i+1)=zu_Sec(1,i)+dz_Sec(:,i); % Future displacement
za_Sec(1,i+1)=a0*dz_Sec(:,i)-a1*zv_Sec(1,i)-a2*za_Sec(1,i); % Future Acceleration
zv_Sec(1,i+1)=zv_Sec(1,i)+a3*za_Sec(1,i)+a4*za_Sec(1,i+1); % Future Velocity

if xc_Mon>=l
j_Mon=j_Mon+1;
J_Mon=J_Mon+1;
xo_Mon=nodes(J_Mon-1,2); % Updated for new element
cor_Mon(j_Mon,:)=ele(j_Mon,2:5); % New element coordinates
end

if xc_Sec>=l
j_Sec=j_Sec+1;
J_Sec=J_Sec+1;
xo_Sec=nodes2(J_Sec-1,2); % Updated for new element
cor_Sec(j_Sec,:)=ele2(j_Sec,2:5); % New element coordinates
end

if i<Start_Time_Last_Vehicle && Vehicle_order(1)==1
% Update global x position for next loop
xg_Mon=xg_Mon+dT*V(Day,ii,1);
elseif i<Start_Time_Last_Vehicle && Vehicle_order(1)==2
% Update global x position for next loop
xg_Sec=xg_Sec-dT*V(Day,ii,2);
elseif i>=Start_Time_Last_Vehicle

    if xg_Mon<L
% Update global x position for next loop
xg_Mon=xg_Mon+dT*V(Day,ii,1);
    else
zu_Mon(:,i+1)=[0;0]; % Future displacement
za_Mon(:,i+1)=[0;0]; % Future Acceleration
zv_Mon(:,i+1)=[0;0]; % Future Velocity

      xg_Mon=L;
    end

    if xg_Sec>0
% Update global x position for next loop
xg_Sec=xg_Sec-dT*V(Day,ii,2);
    else
zu_Sec(:,i+1)=[0;0]; % Future displacement
za_Sec(:,i+1)=[0;0]; % Future Acceleration
zv_Sec(:,i+1)=[0;0]; % Future Velocity
      xg_Sec=0;
    end
end

end % end vehicle loop

% Acceleration plots before excess data is trimmed off the ends
% figure(1)
% set(gcf,'color','white')
% plot(1:Kt,za_Sec(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat
%
% figure(2)
% set(gcf,'color','white')
% plot(1:Kt,za_Mon(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat

%Test matrix for debugging
% Matrix1=[Kt_Mon, Kt_Sec; length(rx_Mon), length(rx_Sec); length(za_Mon), length(za_Sec)]
%  Matrix2=[1:Kt-1;rx_Mon; za_Mon(1,1:end-1); rx_Sec; za_Sec(1,1:end-1)]

if Vehicle_order(1)==1
    if (Start_Time_Last_Vehicle+Kt_Sec)<Kt_Mon
za_Sec(:,Start_Time_Last_Vehicle+Kt_Sec:end)=[];
rx_Sec(:,Start_Time_Last_Vehicle+Kt_Sec-1:end)=[];
drx_Sec(:,Start_Time_Last_Vehicle+Kt_Sec-1:end)=[];
ddrx_Sec(:,Start_Time_Last_Vehicle+Kt_Sec-1:end)=[];
za_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
rx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
drx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
ddrx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
    end

    if (Start_Time_Last_Vehicle+Kt_Sec)>Kt_Mon
za_Sec(:,end)=[];
za_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
rx_Sec(:,end)=[];
rx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
drx_Sec(:,end)=[];
drx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
ddrx_Sec(:,end)=[];
ddrx_Sec(:,1:Start_Time_Last_Vehicle-1)=[];

za_Mon(:,Kt_Mon+1:end)=[];
rx_Mon(:,Kt_Mon:end)=[];
drx_Mon(:,Kt_Mon:end)=[];
ddrx_Mon(:,Kt_Mon:end)=[];
    end

elseif Vehicle_order(1)==2
        if (Start_Time_Last_Vehicle+Kt_Mon)<Kt_Sec
za_Mon(:,Start_Time_Last_Vehicle+Kt_Mon:end)=[];
rx_Mon(:,Start_Time_Last_Vehicle+Kt_Mon-1:end)=[];
drx_Mon(:,Start_Time_Last_Vehicle+Kt_Mon-1:end)=[];
ddrx_Mon(:,Start_Time_Last_Vehicle+Kt_Mon-1:end)=[];
za_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
rx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
drx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
ddrx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
        end

    if (Start_Time_Last_Vehicle+Kt_Mon)>Kt_Sec
za_Mon(:,end)=[];
za_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
rx_Mon(:,end)=[];
rx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
drx_Mon(:,end)=[];
drx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
ddrx_Mon(:,end)=[];
ddrx_Mon(:,1:Start_Time_Last_Vehicle-1)=[];

za_Sec(:,Kt_Sec+1:end)=[];
rx_Sec(:,Kt_Sec:end)=[];
drx_Sec(:,Kt_Sec:end)=[];
ddrx_Sec(:,Kt_Sec:end)=[];
    end

end

% Test matrix for debugging
% Matrix1=[Kt_Mon, Kt_Sec; length(rx_Mon), length(rx_Sec); length(za_Mon), length(za_Sec)]
%  Matrix2=[rx_Mon(1), rx_Mon(end); rx_Sec(1), rx_Sec(end)]

% Acceleration plots after excess data has been trimmed off of the ends
% figure(3)
% set(gcf,'color','white')
% plot(0:dT:L/V(Day,ii,2),za_Sec(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat
%
% figure(4)
% set(gcf,'color','white')
% plot(0:dT:L/V(Day,ii,1),za_Mon(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat

% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
f_Mon = Fs*(0:(Kt_Mon))/(Kt_Mon*2); % Frequency domain
f_Sec = Fs*(0:(Kt_Sec))/(Kt_Sec*2); % Frequency domain

% Executing FFT for Vehicle
fftV_FE_Mon=abs(fft(za_Mon(1,:),2*Kt_Mon));
fftV_FE_Sec=abs(fft(za_Sec(1,:),2*Kt_Sec));

Twosided_FE_Mon = fftV_FE_Mon/Kt_Mon; % two-sided spectrum
onesided_FE_Mon = Twosided_FE_Mon(1:Kt_Mon+1); % Single-sided spectrum
onesided_FE_Mon(2:end-1) = onesided_FE_Mon(2:end-1);

Twosided_FE_Sec = fftV_FE_Sec/Kt_Sec; % two-sided spectrum
onesided_FE_Sec = Twosided_FE_Sec(1:Kt_Sec+1); % Single-sided spectrum
onesided_FE_Sec(2:end-1) = onesided_FE_Sec(2:end-1);

% FFT plot to check is frequency domain looks normal
% Matrix=[rx_Mon(1), rx_Mon(end); rx_Sec(1), rx_Sec(end)]
% figure(5)
% set(gcf,'color','white')
%  semilogy(f_Mon(1:400),onesided_FE_Mon(1:400),'color','b','LineWidth',2);hold on
%   semilogy(f_Sec(1:400),onesided_FE_Sec(1:400),'color','r','LineWidth',2);hold on
% set(gca,'fontsize',16);
% title('Power Spectrum of Vehicle Acceleration','Fontname','Timesnewroman')
% xlabel('Frequency (Hz)','Fontname','Timesnewroman')
% ylabel('Spectrum Amplitude','Fontname','Timesnewroman')
% plotformat

% Frequency information for vehicle
% FV_Mon=ef_Mon/(2*pi)
% FV_Sec=ef_Sec/(2*pi)

% Storage matrices and cell arrays (used to store information for machine learning)
Monitor_Vehicle_Time{Day,ii}=T_Mon;
Monitor_Vehicle_Acceleration{Day,ii}=za_Mon;
Monitor_Vehicle_Frequency_Amp_Data{Day,ii}=onesided_FE_Mon;
Monitor_Vehicle_Frequency_Data{Day,ii}=f_Mon;
Monitor_Vehicle_Road_Profile{Day,ii}=rx_Mon;
Monitor_Vehicle_Derivative_Road_Profile{Day,ii}=drx_Mon;
Monitor_Vehicle_Other_Derivative_Road_Profile{Day,ii}=ddrx_Mon;

Other_Vehicle_Time{Day,ii}=T_Sec;
Other_Vehicle_Acceleration{Day,ii}=za_Sec;
Other_Vehicle_Frequency_Amp_Data{Day,ii}=onesided_FE_Sec;
Other_Vehicle_Frequency_Data{Day,ii}=f_Sec;
Other_Vehicle_Road_Profile{Day,ii}=rx_Sec;
Other_Vehicle_Derivative_Road_Profile{Day,ii}=drx_Sec;
Other_Vehicle_Other_Derivative_Road_Profile{Day,ii}=ddrx_Sec;

%% End of second section of "If Multiple_Vehicle" loop

% elseif Multiple_Vehicles==1 && number_vehicles(Day,ii)==3
%
%% Beginning of last section of "If Multiple_Vehicle" loop (This section is for if there is only one vehicle ever being considered)
else
    % Vehicle parameters
if row(Day,ii)<=3
mv=VehicleVariables(row(Day,ii),1)+randi([-50 50],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw=VehicleVariables(row(Day,ii),2); % wheel mass of vehicle kg
elseif row(Day,ii)<=7
mv=VehicleVariables(row(Day,ii),1)+randi([-500 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw=VehicleVariables(row(Day,ii),2)+randi([-3 3],1,1);
else
mv=VehicleVariables(row(Day,ii),1)+randi([-1000 500],1,1);% sprung mass of vehicle kg (Randomly selects 1 of 10 vehicles)
mw=VehicleVariables(row(Day,ii),2)+randi([-10 10],1,1);
end
kv=VehicleVariables(row(Day,ii),3); %Stiffness of vehicle spring N/m
kw=VehicleVariables(row(Day,ii),6); %Stiffness of vehicle tire N/m
cs=VehicleVariables(row(Day,ii),4); % Damping of vehicle spring N*s/m
cw=VehicleVariables(row(Day,ii),5); % Damping of vehicle tire N*s/m

K_Mon=[kv, -kv; -kv, kv+kw];
M_Mon=[mv, 0; 0, mw];
ei_Mon=eig(K_Mon,M_Mon); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

% Data Storage Matrices
MonitorVehicleMass(Day,ii)=VehicleVariables(row(Day,ii),1);
MonitorWheelMass(Day,ii)=VehicleVariables(row(Day,ii),2);
MonitorSuspensionStiffness(Day,ii)=VehicleVariables(row(Day,ii),3);
MonitorWheelStiffness(Day,ii)=VehicleVariables(row(Day,ii),6);
MonitorSuspensionDamping(Day,ii)=VehicleVariables(row(Day,ii),4);
MonitorWheelDamping(Day,ii)=VehicleVariables(row(Day,ii),5);

% Time and position
dT=.001; % Time Step
T=0:dT:L/V(Day,ii); % Total time to cross bridge
Kt=length(T);
xg=0; % Initial global position
j=1; % Initial row for elemental matrix
J=2; % Initial row for nodal matrix


%% Matricies and Initial Conditions

% Load Matricies
pc=zeros(2*(NumberElements+1),1);
qc=zeros(2*(NumberElements+1),1);

% Vehicle matricies
muu=mv; muw=0;      mwu=0;      mww=mw;
kuu=kv; kuw=-kv;    kwu=-kv;    kww=kv+kw;
cuu=cs; cuw=-cs;    cwu=-cs;    cww=cs+cw;
fue_t_dt=0; fwe_t_dt=-9.81*mv-9.81*mw;
lw=1;

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt);
AB1=zeros(2*(NumberElements+1),Kt);
VB1=zeros(2*(NumberElements+1),Kt);

% Initial Condition Vehicle
zu=zeros(2,Kt-1); % Initial displacement Vehicle
zv=zeros(2,Kt-1);  % Initial velocity of Vehicle
za=zeros(2,Kt-1);  % Initial acceleration of Vehicle

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
wn_FEA=zeros(2*NumberElements,Kt-1);
rc=zeros(2*(NumberElements+1),1);
rx=zeros(1,Kt-1);
drx=zeros(1,Kt-1);
ddrx=zeros(1,Kt-1);

%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
xe=nodes(J,2); % elment end location
cor_Mon(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1 % This loop is for calculating the accleration response for each vehicle crossing

xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
t=xc/V(Day,ii); % Local time
Nc(i,:)=[1-3*xb^2+2*xb^3, xc*(1-2*xb+xb^2), 3*xb^2-2*xb^3, xc*(xb^2-xb)]; % Shape function Row Vector
Ncd(i,:)=[-6*V(Day,ii)^2*t/l^2+6*V(Day,ii)^3*t^2/l^3, V(Day,ii)-4*V(Day,ii)^2*t/l+3*V(Day,ii)^3*t^2/l^2, 6*V(Day,ii)^2*t/l^2-6*V(Day,ii)^3*t^2/l^3, 3*V(Day,ii)^3*t^2/l^2-2*V(Day,ii)^2*t/l];
Ncdd(i,:)=[-6*V(Day,ii)^2/l^2+12*V(Day,ii)^3*t/l^3, -4*V(Day,ii)^2/l+6*V(Day,ii)^3*t/l^2, 6*V(Day,ii)^2/l^2-12*V(Day,ii)^3*t/l^3, 6*V(Day,ii)^3*t/l^2-2*V(Day,ii)^2/l];
Nct=transpose(Nc); % Column Vector

% Surface Rougness
if Surface==1
Column=find(RoadMatrix(1,:)==round(xg,3));
rx(i)=RoadMatrix(2,Column(1,1));
drx(i)=RoadMatrix(3,Column(1,1));
ddrx(i)=RoadMatrix(4,Column(1,1));
rc_tdt_st=Nct(:,i)*(V(Day,ii)^2*mc*ddrx(i)+V(Day,ii)*cc*drx(i)+kc*rx(i)); %Contact force from road irregularities
RC=rc;
rch=KInsert2(rc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
RC=RC+rch;
end

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
kkh=KInsert(kc_st,cor_Mon(j,:),2*(NumberElements+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor_Mon(j,:),2*(NumberElements+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

% Apply boundary conditions to calculate damping matirx
[M,K]=boundarycondition(MB1,KB1,NumberElements);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s]
wn_FEA(:,i)=ef/(2*pi); % sorted natural angular frequencies [Hz]

% Beam damping matrix
al=2*bbeta*wn_FEA(1,i)*wn_FEA(2,i)/(wn_FEA(1,i)+wn_FEA(2,i)); % Alpha for Rayleigh Damping
be=2*bbeta/(wn_FEA(1,i)+wn_FEA(2,i)); % Beta for Rayleigh Damping
C=al*MB1+be*KB1; % Damping Matrix for beam

% Global Damping
cch=KInsert(cc_st,cor_Mon(j,:),2*(NumberElements+1));
C=C+cch; % Updated beam damping matrix
% Apply boundary condition for damping matrix
C(2*NumberElements+1,:)=[];
C(:,2*NumberElements+1)=[];
C(1,:)=[];
C(:,1)=[];

% Global Contact Loads
PC=pc;
pch=KInsert2(pc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
PC=PC+pch;

QC=qc;
qch=KInsert2(qc_t_st,cor_Mon(j,:),2*(NumberElements+1));
QC=QC+qch;

% Apply Boundary Conditions for Load Matrices
% Delete first and next to last row
if Surface==1
[PC,QC,RC]=bcloads(NumberElements,PC,QC,RC);
else
[PC,QC]=bcloads(NumberElements,PC,QC);
end

if Wind==1
    FW(:)=ForceWind(Day,ii);
end

if Wind==1 && Surface==1
RS=FW-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==1 && Surface==0
RS=FW-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
elseif Wind==0 && Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
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
db(:,i)=UB1(cor_Mon(j,:),i); % Local Displacement
vb(:,i)=VB1(cor_Mon(j,:),i); % Local Velocity
ab(:,i)=AB1(cor_Mon(j,:),i); % Local Acceleration

% Contact points
dc(:,i)=Nc(i,:)*db(:,i); % Contact displacement
vc(:,i)=V(Day,ii)*Ncd(i,:)*db(:,i)+Nc(i,:)*vb(:,i); % Contact velocity
ac(:,i)=Nc(i,:)*ab(:,i)+2*V(Day,ii)*Ncd(i,:)*vb(:,i)+(V(Day,ii)^2)*Ncdd(i,:)*db(:,i); % Contact acceleration

if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zu(2,i+1)=dc(:,i)+rx(i); % Future displacement
zv(2,i+1)=vc(:,i)+V(Day,ii)*drx(i); % Future Velocity
za(2,i+1)=ac(:,i)+V(Day,ii)^2*ddrx(i); % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle (without
% surface roughness)
zu(2,i+1)=dc(:,i); % Future displacement
zv(2,i+1)=vc(:,i); % Future Velocity
za(2,i+1)=ac(:,i); % Future Acceleration
end

% Contact Force
% Vi_Tdt(:,i+1)=pc_tdt+qc_t+(mc*ac(:,i)+cc*vc(:,i)+kc*dc(:,i));
quc_tdt(:,i+1)=muw*za(2,i+1)+cuw*zv(2,i+1)+kuw*zu(2,i+1);

% Change in upper vehicle displacement
dz(:,i)=PSIuu\(-quc_tdt(:,i+1)+qu_t);

% Future displacement, vnocity and accneration in upper vehicle
zu(1,i+1)=zu(1,i)+dz(:,i); % Future displacement
za(1,i+1)=a0*dz(:,i)-a1*zv(1,i)-a2*za(1,i); % Future Acceleration
zv(1,i+1)=zv(1,i)+a3*za(1,i)+a4*za(1,i+1); % Future Velocity

if xc>=l
j=j+1;
J=J+1;
xo=nodes(J-1,2); % Updated for new element
xe=nodes(J,2); % Updated for new element
cor_Mon(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V(Day,ii);
end % end single vehicle loop
% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
t = (0:Kt-1)*dT; % Time vector
f = Fs*(0:(Kt))/(Kt*2); % Frequency domain

% Executing FFT for Vehicle
fftV_FE=abs(fft(za(1,:),2*Kt));

Twosided_FE = fftV_FE/Kt; % two-sided spectrum
onesided_FE = Twosided_FE(1:Kt+1); % Single-sided spectrum
onesided_FE(2:end-1) = onesided_FE(2:end-1);

Monitor_Vehicle_Time{Day,ii}=T;
Monitor_Vehicle_Acceleration{Day,ii}=za;
Monitor_Vehicle_Frequency_Amp_Data{Day,ii}=onesided_FE;
Monitor_Vehicle_Frequency_Data{Day,ii}=f;
Monitor_Vehicle_Road_Profile{Day,ii}=rx;
Monitor_Vehicle_Derivative_Road_Profile{Day,ii}=drx;
Monitor_Vehicle_Other_Derivative_Road_Profile{Day,ii}=ddrx;

end


end % end day loop

save(outPath)
exitCode = 0;
end
