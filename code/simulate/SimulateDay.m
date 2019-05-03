function exitCode = SimulateDay(inPath, DaySTR, outPath)
Day = str2num(DaySTR);
load(inPath);

% Record 4 - 9
FFT_STEP = 0.05
FFT_SELECT_START=floor(3/FFT_STEP)+1
FFT_SELECT_END=floor(10/FFT_STEP)+1

hour = 0;
for ii=1:n

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
Q=normrnd(1.0129,.003);
S=normrnd(-.0048,.0001);
R=normrnd(.1977,.0027);
tu=normrnd(3.1466,.0861);
lam=normrnd(-1.1012,.0513);
% Modified Modulus
u0=Q+S*Tact(Day,ii)+R*(1-erf((Tact(Day,ii)-lam)/tu)); % Modification factor
E0=u0*E;
else
u0=1;
E0=E;
end

% Update the hour
if Temp==1
if Td(ii)>=(hour+1)
    hour=hour+1;
end
end

% Applying damaged modulus to elemental matrix
if Damage==1 && Damage_Case==1
  if Day>=DayDamage1
    if sum(ED1(1,1:(Day-DayDamage1+1)))/E >= .25
      E_damaged1=E*.75;
    else
      E_damaged1=(E-sum(ED1(1,1:(Day-DayDamage1+1)))); % Overall Damaged Modulus 1
    end
    kb_damaged1=KBeam(E_damaged1*u0,I,l);
  end
elseif Damage==1 && Damage_Case==2

    if Day>=DayDamage5
E_damaged1=E-ED1-ED2-ED3-ED4-ED5; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
    elseif Day>=DayDamage4
E_damaged1=E-ED1-ED2-ED3-ED4; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
    elseif Day>=DayDamage3
E_damaged1=E-ED1-ED2-ED3; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
    elseif Day>=DayDamage2
E_damaged1=E-ED1-ED2; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
    elseif Day>=DayDamage1
E_damaged1=E-ED1; % Overall Damaged Modulus 1
kb_damaged1=KBeam(E_damaged1*u0,I,l);
    end
end

% Healthy elemental matrices
kb=KBeam(E0,I,l); % Stiffness matrix for bridge
mb=MBeam(mu(1),l); % Consistent mass matrix for bridge

% Global Beam Matricies
KB=zeros(2*(NumberElements+1),2*(NumberElements+1));
MB=zeros(2*(NumberElements+1),2*(NumberElements+1));
cor_Mon=zeros(NumberElements,4);
for i=1:NumberElements
   cor_Mon(i,:)=ele(i,2:5);
   kk=KInsert(kb,cor_Mon(i,:),2*(NumberElements+1));
   if Damage==1 && Damage_Case==1
   if i==DamageLocation % Insert damage state 1
       if Day>=DayDamage1
   kk=KInsert(kb_damaged1,cor_Mon(i,:),2*(NumberElements+1));
       end
   end

   elseif Damage==1 && Damage_Case==2
       if i==DamageLocation
           if Day>=DayDamage5 || Day>=DayDamage4 || Day>=DayDamage3 || Day>=DayDamage2 ||Day>=DayDamage1
    kk=KInsert(kb_damaged1,cor_Mon(i,:),2*(NumberElements+1));
           end
       end
   end

 KB=KB+kk; % Beam stiffness matrix
 mm=KInsert(mb,cor_Mon(i,:),2*(NumberElements+1));
 MB=MB+mm; % Beam mass matrix
end

% Remove boundary conditions only so we can calculate damping matirx

[M,K]=boundarycondition(MB,KB,NumberElements);
% Encountered errors regarding bad eig values
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s]
wn_FEA=ef; % sorted natural angular frequencies

% Beam damping matrix
al=2*bbeta*wn_FEA(1,1)*wn_FEA(2,1)/(wn_FEA(1,1)+wn_FEA(2,1)); % Alpha for Rayleigh Damping
be=2*bbeta/(wn_FEA(1,1)+wn_FEA(2,1)); % Beta for Rayleigh Damping
CB=al*MB+be*KB; % Damping Matrix for beam


%% If statement for number of vehicles

if Multiple_Vehicles==1 && number_vehicles(Day,ii)==1
    % Vehicle parameters
mv_Mon=normrnd(VehicleVariables(row(Day,ii,1),1),VehicleVariables(row(Day,ii,1),2));% sprung mass of vehicle kg (Randomly selects 1 of 8 vehicles)
mw_Mon=normrnd(VehicleVariables(row(Day,ii,1),3),VehicleVariables(row(Day,ii,1),4)); % wheel mass of vehicle kg
kv_Mon=VehicleVariables(row(Day,ii,1),5); %Stiffness of vehicle spring N/m
kw_Mon=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cv_Mon=VehicleVariables(row(Day,ii,1),7); % Damping of vehicle spring N*s/m
cw_Mon=VehicleVariables(row(Day,ii,1),8); % Damping of vehicle tire N*s/m


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
muu=diag([mv_Mon mw_Mon]); muw=[0;0];  mwu=[0,0];  mww=0;
M_Mon=[muu,muw;mwu,mww]; % Vehicle mass matrix
kuu=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
kuw=[0;-kw_Mon]; kwu=[0,-kw_Mon]; kww=kw_Mon;
K_Mon=[kuu,kuw;kwu,kww]; % Vehicle stiffness matrix
cuu=[cv_Mon, -cv_Mon; -cv_Mon, cv_Mon+cw_Mon];
cuw=[0;-cw_Mon]; cwu=[0,-cw_Mon]; cww=cw_Mon;
fue_t_dt=[0;0]; fwe_t_dt=-9.81*mv_Mon-9.81*mw_Mon;
lw=1;

ei_Mon=eig(kuu,muu); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt);
AB1=zeros(2*(NumberElements+1),Kt);
VB1=zeros(2*(NumberElements+1),Kt);

% Initial Condition Vehicle
zu=zeros(2,Kt-1); % Initial displacement upper Vehicle (vertical, rotational, wheel 1, wheel 2)
zw=zeros(1,Kt-1); % Initial displacement lower Vehicle (veritcal at both tires)
zv_u=zeros(2,Kt-1);  % Initial velocity of upper Vehicle
zv_w=zeros(1,Kt-1);  % Initial velocity of lower Vehicle
za_u=zeros(2,Kt-1);  % Initial acceleration of upper Vehicle
za_w=zeros(1,Kt-1);  % Initial acceleration of lower Vehicle

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
mc=lw\(mww-PSIwu/PSIuu*muw); %Mass contact matrix
cc=lw\(cww-PSIwu/PSIuu*cuw); %Damping contact matrix
kc=lw\(kww-PSIwu/PSIuu*kuw); %Stiffness contact matrix

%% Pre-allocate matrix sizes
rc=zeros(2*(NumberElements+1),1);

%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
cor_Mon(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1 % This loop is for calculating the accleration response for each vehicle crossing

xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
Nc=[1-3*xb.^2+2*xb.^3, xc.*(1-2*xb+xb.^2), 3*xb.^2-2*xb.^3, xc.*(xb.^2-xb)];
Ncd=[-6*(xc/l^2)+6*(xc.^2/l^3), 1-4*xb+3*xb.^2, 6*(xc/l^2)-6*(xc.^2/l^3), 3*xb.^2-2*xb];
Ncdd=[-6/l^2+12*xc/l^3, -4/l+6*xc/l^2, 6/l^2-12*xc/l^3, 6*xc/l^2-2/l];
Nct=transpose(Nc); % Column Vector

% Surface Rougness
if Surface==1
difference = abs(RoadMatrix(1,:)-xg);
[~, Column] = min(difference);
rx=RoadMatrix(2,Column(1,1));
drx=RoadMatrix(3,Column(1,1));
ddrx=RoadMatrix(4,Column(1,1));
rc_tdt_st=Nct*(V(Day,ii,1)^2*mc*ddrx+V(Day,ii,1)*cc*drx+kc*rx); %Contact force from road irregularities
RC=rc;
rch=KInsert2(rc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
RC=RC+rch;
end

% Load Vectors
qu_t=muu*(a1*zv_u(:,i)+a2*za_u(:,i))+cuu*(a6*zv_u(:,i)+a7*za_u(:,i))-kuu*zu(:,i);
qw_t=mwu*(a1*zv_u(:,i)+a2*za_u(:,i))+cwu*(a6*zv_u(:,i)+a7*za_u(:,i))-kwu*zu(:,i);
pc_tdt=lw\(PSIwu/PSIuu*fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu/PSIuu*qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct*mc*Nc; %Mass
cc_st=Nct*(2*V(Day,ii,1)*mc*Ncd+cc*Nc); %Damping
kc_st=Nct*(V(Day,ii,1)^2*mc*Ncdd+V(Day,ii,1)*cc*Ncd+kc*Nc); %Stiffness

% Equivalent nodal loads
pc_tdt_st=Nct*pc_tdt;
qc_t_st=Nct*qc_t;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh=KInsert(kc_st,cor_Mon(j,:),2*(NumberElements+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor_Mon(j,:),2*(NumberElements+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

CB1=CB; % resets global matrix each time
cch=KInsert(cc_st,cor_Mon(j,:),2*(NumberElements+1));
CB1=CB1+cch; % Updated beam stiffness matrix

% Apply boundary conditions
[M,K,C]=boundarycondition(MB1,KB1,NumberElements,CB1);

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

if Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
RS=-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
end

% Finding unknown displacements
du=LS\RS;

% Future Disp, Vel, Acc
UB(:,i+1)=UB(:,i)+du; % Future displacement
AB(:,i+1)=a0*du-a1*VB(:,i)-a2*AB(:,i); % Future Acceleration
VB(:,i+1)=VB(:,i)+a3*AB(:,i)+a4*AB(:,i+1); % Future Velocity

% Add back constrained dof
UB1(2:2*NumberElements,i)=UB(1:2*NumberElements-1,i+1);
UB1(2*(NumberElements+1),i)=UB(2*NumberElements,i+1);
AB1(2:2*NumberElements,i)=AB(1:2*NumberElements-1,i+1);
AB1(2*(NumberElements+1),i)=AB(2*NumberElements,i+1);
VB1(2:2*NumberElements,i)=VB(1:2*NumberElements-1,i+1);
VB1(2*(NumberElements+1),i)=VB(2*NumberElements,i+1);

% Element displacement, velocity, acceleration
db=UB1(cor_Mon(j,:),i); % Local Displacement
vb=VB1(cor_Mon(j,:),i); % Local Velocity
ab=AB1(cor_Mon(j,:),i); % Local Acceleration

% Contact points
dc=Nc*db; % Contact displacement
vc=V(Day,ii,1)*Ncd*db+Nc*vb; % Contact velocity
ac=Nc*ab+2*V(Day,ii,1)*Ncd*vb+(V(Day,ii,1)^2)*Ncdd*db; % Contact acceleration

if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zw(:,i+1)=dc+rx; % Future displacement
zv_w(:,i+1)=vc+V(Day,ii,1)*drx; % Future Velocity
za_w(:,i+1)=ac+V(Day,ii,1)^2*ddrx; % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle
zw(:,i+1)=dc; % Future displacement
za_w(:,i+1)=ac; % Future Acceleration
zv_w(:,i+1)=vc; % Future Velocity
end

% Contact Force
quc_tdt=muw*za_w(:,i+1)+cuw*zv_w(:,i+1)+kuw*zw(:,i+1);

% Change in upper vehicle displacement
dz=PSIuu\(fue_t_dt-quc_tdt+qu_t);

% Future displacement, vnocity and accneration in upper vehicle
zu(:,i+1)=zu(:,i)+dz; % Future displacement
za_u(:,i+1)=a0*dz-a1*zv_u(:,i)-a2*za_u(:,i); % Future Acceleration
zv_u(:,i+1)=zv_u(:,i)+a3*za_u(:,i)+a4*za_u(:,i+1); % Future Velocity

if xc>=l
j=j+1;
J=J+1;
xo=nodes(J-1,2); % Updated for new element
cor_Mon(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V(Day,ii,1);
end % end single vehicle loop

% Adding noise to the acceleration of upper vehicle
% za_u=za_u+(za_u*(-.025))+randn(2,1).*(za_u*(.025)-za_u*(-.025));

% figure(4)
% set(gcf,'color','white')
% plot(T,za_u(1,:),'linewidth',3);hold on
% title('Acceleration of Upper Vehicle')
% set(gca,'fontsize',16);
% xlabel('Time (s) ');
% ylabel('Displacement (m)');
% % legend('11.11','22.22','33.33','44.44','location','northwest')
% plotformat
%
% figure(5)
% set(gcf,'color','white')
% plot(T,za_u(2,:),'linewidth',3);hold on
% title('Acceleration of Upper Vehicle')
% set(gca,'fontsize',16);
% xlabel('Time (s) ');
% ylabel('Displacement (m)');
% % legend('11.11','22.22','33.33','44.44','location','northwest')
% plotformat


% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
nFFT = Fs/FFT_STEP;
if rem(nFFT,2)>0
    nFFT = nFFT+1;
end
f = Fs*(0:(nFFT/2))/(nFFT); % Frequency domain

za_u = lowpass(za_u',100,Fs,'Steepness',.99)';

% Executing FFt for Upper Vehicle
fft_V=abs(fft(za_u,nFFT,2)/nFFT);
onesided_FE_Original = fft_V(:,1:nFFT/2+1); % Single-sided spectrum
onesided_FE_Original(:,2:end-1) = 2*onesided_FE_Original(:,2:end-1); % Scale Power by 2

% Apply Filters
% Bandpass
%fpass=[4 9];
%Filt_BF = bandpass(za_u',fpass,Fs,'Steepness',.99);
%fft_V = abs(fft(Filt_BF',nFFT,2)/nFFT);
%onesided_FE_BF = fft_V(:,1:nFFT/2+1); % Single-sided spectrum
%onesided_FE_BF(:,2:end-1) = 2*onesided_FE_BF(:,2:end-1); % Scale Power by 2

% % Storage matrices and cell arrays (used to store information for machine
% % learning)
Monitor_Vehicle_Time{Day,ii}=T;
Monitor_Vehicle_Acceleration{Day,ii}=za_u;
Monitor_Vehicle_Frequency_Amp_Data.Original{Day,ii}=onesided_FE_Original(:, FFT_SELECT_START:FFT_SELECT_END);
%Monitor_Vehicle_Frequency_Amp_Data.Filtered{Day,ii}=onesided_FE_BF;
Monitor_Vehicle_Frequency_Data{Day,ii}=f(FFT_SELECT_START:FFT_SELECT_END); % only store band

%% End of first section of "If Multiple_Vehicle" loop

%% Beginning of Second section of "If Multiple_Vehicle" loop
elseif Multiple_Vehicles==1 && number_vehicles(Day,ii)==2

     % Vehicle parameters
mv_Mon=normrnd(VehicleVariables(row(Day,ii,1),1),VehicleVariables(row(Day,ii,1),2));% sprung mass of vehicle kg (Randomly selects 1 of 8 vehicles)
mw_Mon=normrnd(VehicleVariables(row(Day,ii,1),3),VehicleVariables(row(Day,ii,1),4)); % wheel mass of vehicle kg
kv_Mon=VehicleVariables(row(Day,ii,1),5); %Stiffness of vehicle spring N/m
kw_Mon=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cv_Mon=VehicleVariables(row(Day,ii,1),7); % Damping of vehicle spring N*s/m
cw_Mon=VehicleVariables(row(Day,ii,1),8); % Damping of vehicle tire N*s/m

mv_Sec=normrnd(VehicleVariables(row(Day,ii,1),1),VehicleVariables(row(Day,ii,1),2));% sprung mass of vehicle kg (Randomly selects 1 of 8 vehicles)
mw_Sec=normrnd(VehicleVariables(row(Day,ii,1),3),VehicleVariables(row(Day,ii,1),4)); % wheel mass of vehicle kg
kv_Sec=VehicleVariables(row(Day,ii,1),5); %Stiffness of vehicle spring N/m
kw_Sec=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cv_Sec=VehicleVariables(row(Day,ii,1),7); % Damping of vehicle spring N*s/m
cw_Sec=VehicleVariables(row(Day,ii,1),8); % Damping of vehicle tire N*s/m

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
muu_Mon=diag([mv_Mon mw_Mon]); muw_Mon=[0;0];  mwu_Mon=[0,0];  mww_Mon=0;
kuu_Mon=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
kuw_Mon=[0;-kw_Mon]; kwu_Mon=[0,-kw_Mon]; kww_Mon=kw_Mon;
cuu_Mon=[cv_Mon, -cv_Mon; -cv_Mon, cv_Mon+cw_Mon];
cuw_Mon=[0;-cw_Mon]; cwu_Mon=[0,-cw_Mon]; cww_Mon=cw_Mon;
fue_t_dt_Mon=[0;0]; fwe_t_dt_Mon=-9.81*mv_Mon-9.81*mw_Mon;

ei_Mon=eig(kuu_Mon,muu_Mon); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

muu_Sec=diag([mv_Sec mw_Sec]); muw_Sec=[0;0];  mwu_Sec=[0,0];  mww_Sec=0;
kuu_Sec=[kv_Sec, -kv_Sec; -kv_Sec, kv_Sec+kw_Sec];
kuw_Sec=[0;-kw_Sec]; kwu_Sec=[0,-kw_Sec]; kww_Sec=kw_Sec;
cuu_Sec=[cv_Sec, -cv_Sec; -cv_Sec, cv_Sec+cw_Sec];
cuw_Sec=[0;-cw_Sec]; cwu_Sec=[0,-cw_Sec]; cww_Sec=cw_Sec;
fue_t_dt_Sec=[0;0]; fwe_t_dt_Sec=-9.81*mv_Sec-9.81*mw_Sec;
lw=1;

ei_Sec=eig(kuu_Sec,muu_Sec); % eigenvalues
ef_Sec=sort(real(sqrt(ei_Sec))); % sorted natural angular frequencies [rad/s]
Secondfv{Day,ii}=ef_Sec/(2*pi); % 1st and second natural frequencies of sprung mass

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt-1);
AB1=zeros(2*(NumberElements+1),Kt-1);
VB1=zeros(2*(NumberElements+1),Kt-1);

% Initial Condition Vehicles
zu_Mon=zeros(2,Kt-1); % Initial displacement upper Vehicle
zw_Mon=zeros(1,Kt-1); % Initial displacement lower Vehicle
zv_u_Mon=zeros(2,Kt-1);  % Initial velocity of upper Vehicle
zv_w_Mon=zeros(1,Kt-1);  % Initial velocity of lower Vehicle
za_u_Mon=zeros(2,Kt-1);  % Initial acceleration of upper Vehicle
za_w_Mon=zeros(1,Kt-1);  % Initial acceleration of lower Vehicle

zu_Sec=zeros(2,Kt-1); % Initial displacement upper Vehicle
zw_Sec=zeros(1,Kt-1); % Initial displacement lower Vehicle
zv_u_Sec=zeros(2,Kt-1);  % Initial velocity of upper Vehicle
zv_w_Sec=zeros(1,Kt-1);  % Initial velocity of lower Vehicle
za_u_Sec=zeros(2,Kt-1);  % Initial acceleration of upper Vehicle
za_w_Sec=zeros(1,Kt-1);  % Initial acceleration of lower Vehicle

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
mc_Mon=lw\(mww_Mon-PSIwu_Mon/PSIuu_Mon*muw_Mon); %Monitoring Mass contact matrix
cc_Mon=lw\(cww_Mon-PSIwu_Mon/PSIuu_Mon*cuw_Mon); %Monitoring Damping contact matrix
kc_Mon=lw\(kww_Mon-PSIwu_Mon/PSIuu_Mon*kuw_Mon); %Monitoring Stiffness contact matrix

PSIuu_Sec=a0*muu_Sec+a5*cuu_Sec+kuu_Sec;
PSIwu_Sec=a0*mwu_Sec+a5*cwu_Sec+kwu_Sec;
mc_Sec=lw\(mww_Sec-PSIwu_Sec/PSIuu_Sec*muw_Sec); %Monitoring Mass contact matrix
cc_Sec=lw\(cww_Sec-PSIwu_Sec/PSIuu_Sec*cuw_Sec); %Monitoring Damping contact matrix
kc_Sec=lw\(kww_Sec-PSIwu_Sec/PSIuu_Sec*kuw_Sec); %Monitoring Stiffness contact matrix

%% Pre-allocate matrix sizes
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

Nc_Mon=[1-3*xb_Mon.^2+2*xb_Mon.^3, xc_Mon.*(1-2*xb_Mon+xb_Mon.^2), 3*xb_Mon.^2-2*xb_Mon.^3, xc_Mon.*(xb_Mon.^2-xb_Mon)];
Ncd_Mon=[-6*(xc_Mon/l^2)+6*(xc_Mon.^2/l^3), 1-4*xb_Mon+3*xb_Mon.^2, 6*(xc_Mon/l^2)-6*(xc_Mon.^2/l^3), 3*xb_Mon.^2-2*xb_Mon];
Ncdd_Mon=[-6/l^2+12*xc_Mon/l^3, -4/l+6*xc_Mon/l^2, 6/l^2-12*xc_Mon/l^3, 6*xc_Mon/l^2-2/l];
Nct_Mon=transpose(Nc_Mon); % Column Vector

xc_Sec=0; % Local position of vehicle on bridge
Nc_Sec=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Sec=[0, 0, 0, 0];
Ncdd_Sec=[0, 0, 0, 0];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    else
xc_Sec=abs((xg_Sec-xo_Sec)); % Local position of vehicle on bridge
xb_Sec=xc_Sec/l; % Local coordiante

Nc_Sec=[1-3*xb_Sec.^2+2*xb_Sec.^3, xc_Sec.*(1-2*xb_Sec+xb_Sec.^2), 3*xb_Sec.^2-2*xb_Sec.^3, xc_Sec.*(xb_Sec.^2-xb_Sec)];
Ncd_Sec=[-6*(xc_Sec/l^2)+6*(xc_Sec.^2/l^3), 1-4*xb_Sec+3*xb_Sec.^2, 6*(xc_Sec/l^2)-6*(xc_Sec.^2/l^3), 3*xb_Sec.^2-2*xb_Sec];
Ncdd_Sec=[-6/l^2+12*xc_Sec/l^3, -4/l+6*xc_Sec/l^2, 6/l^2-12*xc_Sec/l^3, 6*xc_Sec/l^2-2/l];
Nct_Sec=transpose(Nc_Sec); % Column Vector

xc_Mon=0; % Local position of vehicle on bridge
Nc_Mon=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Mon=[0, 0, 0, 0];
Ncdd_Mon=[0, 0, 0, 0];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    end

else

    if xg_Mon<L
xc_Mon=(xg_Mon-xo_Mon); % Local position of vehicle on bridge
xb_Mon=xc_Mon/l; % Local coordinate
Nc_Mon=[1-3*xb_Mon.^2+2*xb_Mon.^3, xc_Mon.*(1-2*xb_Mon+xb_Mon.^2), 3*xb_Mon.^2-2*xb_Mon.^3, xc_Mon.*(xb_Mon.^2-xb_Mon)];
Ncd_Mon=[-6*(xc_Mon/l^2)+6*(xc_Mon.^2/l^3), 1-4*xb_Mon+3*xb_Mon.^2, 6*(xc_Mon/l^2)-6*(xc_Mon.^2/l^3), 3*xb_Mon.^2-2*xb_Mon];
Ncdd_Mon=[-6/l^2+12*xc_Mon/l^3, -4/l+6*xc_Mon/l^2, 6/l^2-12*xc_Mon/l^3, 6*xc_Mon/l^2-2/l];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    else
c_Mon=0; % Local position of vehicle on bridge
Nc_Mon=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Mon=[0, 0, 0, 0];
Ncdd_Mon=[0, 0, 0, 0];
Nct_Mon=transpose(Nc_Mon); % Column Vector
    end

    if xg_Sec > 0
xc_Sec=abs((xg_Sec-xo_Sec)); % Local position of vehicle on bridge
xb_Sec=xc_Sec/l; % Local coordiante
Nc_Sec=[1-3*xb_Sec.^2+2*xb_Sec.^3, xc_Sec.*(1-2*xb_Sec+xb_Sec.^2), 3*xb_Sec.^2-2*xb_Sec.^3, xc_Sec.*(xb_Sec.^2-xb_Sec)];
Ncd_Sec=[-6*(xc_Sec/l^2)+6*(xc_Sec.^2/l^3), 1-4*xb_Sec+3*xb_Sec.^2, 6*(xc_Sec/l^2)-6*(xc_Sec.^2/l^3), 3*xb_Sec.^2-2*xb_Sec];
Ncdd_Sec=[-6/l^2+12*xc_Sec/l^3, -4/l+6*xc_Sec/l^2, 6/l^2-12*xc_Sec/l^3, 6*xc_Sec/l^2-2/l];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    else
xc_Sec=0; % Local position of vehicle on bridge
Nc_Sec=[0, 0, 0, 0]; % Shape function Row Vector
Ncd_Sec=[0, 0, 0, 0];
Ncdd_Sec=[0, 0, 0, 0];
Nct_Sec=transpose(Nc_Sec); % Column Vector
    end
end


if Surface==1
RC=rc; % Resets Surface Roughness matrix each loop
    if i<Start_Time_Last_Vehicle && Vehicle_order(1)==2
rx_Mon=0;
drx_Mon=0;
ddrx_Mon=0;
    else
        if xg_Mon<=L
difference = abs(RoadMatrix(1,:)-xg_Mon);
[~, Column_Mon] = min(difference);
rx_Mon=RoadMatrix(2,Column_Mon(1,1));
drx_Mon=RoadMatrix(3,Column_Mon(1,1));
ddrx_Mon=RoadMatrix(4,Column_Mon(1,1));
rc_tdt_st_Mon=Nct_Mon*(V(Day,ii,1)^2*mc_Mon*ddrx_Mon+V(Day,ii,1)*cc_Mon*drx_Mon+kc_Mon*rx_Mon); %Contact force from road irregularities
rch=KInsert2(rc_tdt_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
RC=RC+rch;
        else
rx_Mon=RoadMatrix(2,end);
drx_Mon=RoadMatrix(3,end);
ddrx_Mon=RoadMatrix(4,end);
        end
    end


    if i<Start_Time_Last_Vehicle && Vehicle_order(1)==1
rx_Sec=0;
drx_Sec=0;
ddrx_Sec=0;
    else
        if xg_Sec>=0
difference = abs(RoadMatrix(1,:)-xg_Sec);
[~, Column_Sec] = min(difference);
rx_Sec=RoadMatrix(2,Column_Sec(1,1));
drx_Sec=RoadMatrix(3,Column_Sec(1,1));
ddrx_Sec=RoadMatrix(4,Column_Sec(1,1));
rc_tdt_st_Sec=Nct_Sec*(V(Day,ii,2)^2*mc_Sec*ddrx_Sec+V(Day,ii,2)*cc_Sec*drx_Sec+kc_Sec*rx_Sec); %Contact force from road irregularities
rch=KInsert2(rc_tdt_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
RC=RC+rch;
        else
rx_Sec=RoadMatrix(2,1);
drx_Sec=RoadMatrix(3,1);
ddrx_Sec=RoadMatrix(4,1);
        end
    end
end

% Load Vectors
qu_t_Mon=muu_Mon*(a1*zv_u_Mon(:,i)+a2*za_u_Mon(:,i))+cuu_Mon*(a6*zv_u_Mon(:,i)+a7*za_u_Mon(:,i))-kuu_Mon*zu_Mon(:,i);
qw_t_Mon=mwu_Mon*(a1*zv_u_Mon(:,i)+a2*za_u_Mon(:,i))+cwu_Mon*(a6*zv_u_Mon(:,i)+a7*za_u_Mon(:,i))-kwu_Mon*zu_Mon(:,i);
pc_tdt_Mon=lw\(PSIwu_Mon/PSIuu_Mon*fue_t_dt_Mon-fwe_t_dt_Mon);
qc_t_Mon=lw\(PSIwu_Mon/PSIuu_Mon*qu_t_Mon-qw_t_Mon);

qu_t_Sec=muu_Sec*(a1*zv_u_Sec(:,i)+a2*za_u_Sec(:,i))+cuu_Sec*(a6*zv_u_Sec(:,i)+a7*za_u_Sec(:,i))-kuu_Sec*zu_Sec(:,i);
qw_t_Sec=mwu_Sec*(a1*zv_u_Sec(:,i)+a2*za_u_Sec(:,i))+cwu_Sec*(a6*zv_u_Sec(:,i)+a7*za_u_Sec(:,i))-kwu_Sec*zu_Sec(:,i);
pc_tdt_Sec=lw\(PSIwu_Sec/PSIuu_Sec*fue_t_dt_Sec-fwe_t_dt_Sec);
qc_t_Sec=lw\(PSIwu_Sec/PSIuu_Sec*qu_t_Sec-qw_t_Sec);

% Interaction contact matrices
mc_st_Mon=Nct_Mon*mc_Mon*Nc_Mon; %Mass
cc_st_Mon=Nct_Mon*(2*V(Day,ii,1)*mc_Mon*Ncd_Mon+cc_Mon*Nc_Mon); %Damping
kc_st_Mon=Nct_Mon*(V(Day,ii,1)^2*mc_Mon*Ncdd_Mon+V(Day,ii,1)*cc_Mon*Ncd_Mon+kc_Mon*Nc_Mon); %Stiffness

mc_st_Sec=Nct_Sec*mc_Sec*Nc_Sec; %Mass
cc_st_Sec=Nct_Sec*(2*V(Day,ii,2)*mc_Sec*Ncd_Sec+cc_Sec*Nc_Sec); %Damping
kc_st_Sec=Nct_Sec*(V(Day,ii,2)^2*mc_Sec*Ncdd_Sec+V(Day,ii,2)*cc_Sec*Ncd_Sec+kc_Sec*Nc_Sec); %Stiffness

% Equivalent nodal loads
pc_tdt_st_Mon=Nct_Mon*pc_tdt_Mon;
qc_t_st_Mon=Nct_Mon*qc_t_Mon;

pc_tdt_st_Sec=Nct_Sec*pc_tdt_Sec;
qc_t_st_Sec=Nct_Sec*qc_t_Sec;

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

% Global Damping
CB1=CB; % resets global matrix each time
cch_Mon=KInsert(cc_st_Mon,cor_Mon(j_Mon,:),2*(NumberElements+1));
cch_Sec=KInsert(cc_st_Sec,cor_Sec(j_Sec,:),2*(NumberElements+1));
CB1=CB1+cch_Mon+cch_Sec; % Updated beam damping matrix

% Apply boundary conditions
[M,K,C]=boundarycondition(MB1,KB1,NumberElements,CB1);

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

if Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
RS=-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
end

% Finding unknown displacements
du=LS\RS;

% Future Disp, Vel, Acc
UB(:,i+1)=UB(:,i)+du; % Future displacement
AB(:,i+1)=a0*du-a1*VB(:,i)-a2*AB(:,i); % Future Acceleration
VB(:,i+1)=VB(:,i)+a3*AB(:,i)+a4*AB(:,i+1); % Future Velocity

% Add back constrained dof
UB1(2:2*NumberElements,i)=UB(1:2*NumberElements-1,i);
UB1(2*(NumberElements+1),i)=UB(2*NumberElements,i);
AB1(2:2*NumberElements,i)=AB(1:2*NumberElements-1,i);
AB1(2*(NumberElements+1),i)=AB(2*NumberElements,i);
VB1(2:2*NumberElements,i)=VB(1:2*NumberElements-1,i);
VB1(2*(NumberElements+1),i)=VB(2*NumberElements,i);

% Element displacement, velocity, acceleration
db_Mon=UB1(cor_Mon(j_Mon,:),i); % Local Displacement
vb_Mon=VB1(cor_Mon(j_Mon,:),i); % Local Velocity
ab_Mon=AB1(cor_Mon(j_Mon,:),i); % Local Acceleration

db_Sec=UB1(cor_Sec(j_Sec,:),i); % Local Displacement
vb_Sec=VB1(cor_Sec(j_Sec,:),i); % Local Velocity
ab_Sec=AB1(cor_Sec(j_Sec,:),i); % Local Acceleration

% Contact points
dc_Mon=Nc_Mon*db_Mon; % Contact displacement
vc_Mon=V(Day,ii,1)*Ncd_Mon*db_Mon+Nc_Mon*vb_Mon; % Contact velocity
ac_Mon=Nc_Mon*ab_Mon+2*V(Day,ii,1)*Ncd_Mon*vb_Mon+(V(Day,ii,1)^2)*Ncdd_Mon*db_Mon; % Contact acceleration

dc_Sec=Nc_Sec*db_Sec; % Contact displacement
vc_Sec=V(Day,ii,2)*Ncd_Sec*db_Sec+Nc_Sec*vb_Sec; % Contact velocity
ac_Sec=Nc_Sec*ab_Sec+2*V(Day,ii,2)*Ncd_Sec*vb_Sec+(V(Day,ii,2)^2)*Ncdd_Sec*db_Sec; % Contact acceleration


if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zw_Mon(:,i+1)=dc_Mon+rx_Mon; % Future displacement
zv_w_Mon(:,i+1)=vc_Mon+V(Day,ii,2)*drx_Mon; % Future Velocity
za_w_Mon(:,i+1)=ac_Mon+V(Day,ii,2)^2*ddrx_Mon; % Future Acceleration

% Future displacement, velocity and acceleration in lower vehicle
zw_Sec(:,i+1)=dc_Sec+rx_Sec; % Future displacement
zv_w_Sec(:,i+1)=vc_Sec+V(Day,ii,2)*drx_Sec; % Future Velocity
za_w_Sec(:,i+1)=ac_Sec+V(Day,ii,2)^2*ddrx_Sec; % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle
zw_Mon(:,i+1)=dc_Mon; % Future displacement
zv_w_Mon(:,i+1)=vc_Mon; % Future Velocity
za_w_Mon(:,i+1)=ac_Mon; % Future Acceleration

% Future displacement, velocity and acceleration in lower vehicle
zw_Sec(:,i+1)=dc_Sec; % Future displacement
zv_w_Sec(:,i+1)=vc_Sec; % Future Velocity
za_w_Sec(:,i+1)=ac_Sec; % Future Acceleration
end

% Contact Force
quc_tdt_Mon=muw_Mon*za_w_Mon(:,i+1)+cuw_Mon*zv_w_Mon(:,i+1)+kuw_Mon*zw_Mon(:,i+1);
quc_tdt_Sec=muw_Sec*za_w_Sec(:,i+1)+cuw_Sec*zv_w_Sec(:,i+1)+kuw_Sec*zw_Sec(:,i+1);

% Change in upper vehicle displacement
dz_Mon=PSIuu_Mon\(fue_t_dt_Mon-quc_tdt_Mon+qu_t_Mon);
dz_Sec=PSIuu_Sec\(fue_t_dt_Sec-quc_tdt_Sec+qu_t_Sec);

% Future displacement, vnocity and accneration in upper vehicle
zu_Mon(:,i+1)=zu_Mon(:,i)+dz_Mon; % Future displacement
za_u_Mon(:,i+1)=a0*dz_Mon-a1*zv_u_Mon(:,i)-a2*za_u_Mon(:,i); % Future Acceleration
zv_u_Mon(:,i+1)=zv_u_Mon(:,i)+a3*za_u_Mon(:,i)+a4*za_u_Mon(:,i+1); % Future Velocity

zu_Sec(:,i+1)=zu_Sec(:,i)+dz_Sec; % Future displacement
za_u_Sec(:,i+1)=a0*dz_Sec-a1*zv_u_Sec(:,i)-a2*za_u_Sec(:,i); % Future Acceleration
zv_u_Sec(:,i+1)=zv_u_Sec(:,i)+a3*za_u_Sec(:,i)+a4*za_u_Sec(:,i+1); % Future Velocity

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

zu_Mon(:,i+1)=zeros(2,1); % Future displacement upper Vehicle
zw_Mon(:,i+1)=zeros(1,1); % Future displacement lower Vehicle
zv_u_Mon(:,i+1)=zeros(2,1);  % Future Velocity of upper Vehicle
zv_w_Mon(:,i+1)=zeros(1,1);  % Future Velocity of lower Vehicle
za_u_Mon(:,i+1)=zeros(2,1);  % Future Acceleration of upper Vehicle
za_w_Mon(:,i+1)=zeros(1,1);  % Future Acceleration of lower Vehicle

      xg_Mon=L;
    end

    if xg_Sec>0
% Update global x position for next loop
xg_Sec=xg_Sec-dT*V(Day,ii,2);
    else
zu_Sec(:,i+1)=zeros(2,1); % Future displacement upper Vehicle
zw_Sec(:,i+1)=zeros(1,1); % Future displacement lower Vehicle
zv_u_Sec(:,i+1)=zeros(2,1);  % Future Velocity of upper Vehicle
zv_w_Sec(:,i+1)=zeros(1,1);  % Future Velocity of lower Vehicle
za_u_Sec(:,i+1)=zeros(2,1);  % Future Acceleration of upper Vehicle
za_w_Sec(:,i+1)=zeros(1,1);  % Future Acceleration of lower Vehicle
      xg_Sec=0;
    end
end

end % end vehicle loop

% Adding noise to the acceleration of vehicles
% za_u_Mon=za_u_Mon+(za_u_Mon*(-.025))+randn(2,1).*(za_u_Mon*(.025)-za_u_Mon*(-.025));
% za_u_Sec=za_u_Sec+(za_u_Sec*(-.025))+randn(2,1).*(za_u_Sec*(.025)-za_u_Sec*(-.025)); % Future Acceleration

% Acceleration plots before excess data is trimmed off the ends
% figure(1)
% set(gcf,'color','white')
% plot(1:Kt,za_Sec(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat
% %
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
za_u_Sec(:,Start_Time_Last_Vehicle+Kt_Sec:end)=[];
za_u_Sec(:,1:Start_Time_Last_Vehicle-1)=[];
    end

    if (Start_Time_Last_Vehicle+Kt_Sec)>Kt_Mon
za_u_Sec(:,end)=[];
za_u_Sec(:,1:Start_Time_Last_Vehicle-1)=[];

za_u_Mon(:,Kt_Mon+1:end)=[];
    end

elseif Vehicle_order(1)==2
        if (Start_Time_Last_Vehicle+Kt_Mon)<Kt_Sec
za_u_Mon(:,Start_Time_Last_Vehicle+Kt_Mon:end)=[];
za_u_Mon(:,1:Start_Time_Last_Vehicle-1)=[];
        end

    if (Start_Time_Last_Vehicle+Kt_Mon)>Kt_Sec
za_u_Mon(:,end)=[];
za_u_Mon(:,1:Start_Time_Last_Vehicle-1)=[];

za_u_Sec(:,Kt_Sec+1:end)=[];
    end

end

% Test matrix for debugging
% Matrix1=[Kt_Mon, Kt_Sec; length(rx_Mon), length(rx_Sec); length(za_Mon), length(za_Sec)]
%  Matrix2=[rx_Mon(1), rx_Mon(end); rx_Sec(1), rx_Sec(end)]

% Acceleration plots after excess data has been trimmed off of the ends
% figure(3)
% set(gcf,'color','white')
% plot(0:dT:L/V(Day,ii,2),za_u_Sec(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat
% %
% figure(4)
% set(gcf,'color','white')
% plot(0:dT:L/V(Day,ii,1),za_u_Mon(1,:),'b','linewidth',3);hold on
% title('Monitoring Vehicle Acceleration')
% xlabel('Time (s) ');
% ylabel('Acceleration (m)');
% plotformat

% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
nFFT = Fs/FFT_STEP;
if rem(nFFT,2)>0
    nFFT = nFFT+1;
end
f = Fs*(0:(nFFT/2))/(nFFT); % Frequency domain

za_u_Mon = lowpass(za_u_Mon',100,Fs,'Steepness',.99)';
za_u_Sec = lowpass(za_u_Sec',100,Fs,'Steepness',.99)';

% Executing FFt for Upper Vehicle
fft_V=abs(fft(za_u_Mon,nFFT,2)/nFFT);
onesided_FE_Mon_Original = fft_V(:,1:nFFT/2+1); % Single-sided spectrum

fft_V=abs(fft(za_u_Sec,nFFT,2)/nFFT);
onesided_FE_Sec_Original = fft_V(:,1:nFFT/2+1); % Single-sided spectrum

% Apply Filters
% Bandpass
%fpass=[4 9];
%Filt_BF = bandpass(za_u_Mon',fpass,Fs,'Steepness',.99);
%fft_V = abs(fft(Filt_BF',nFFT,2)/nFFT);
%onesided_FE_Mon_BF = fft_V(:,1:nFFT/2+1); % Single-sided spectrum

%Filt_BF = bandpass(za_u_Sec',fpass,Fs,'Steepness',.99);
%fft_V = abs(fft(Filt_BF',nFFT,2)/nFFT);
%onesided_FE_Sec_BF = fft_V(:,1:nFFT/2+1); % Single-sided spectrum

% figure(6)
% subplot(2,1,1)
% plot(f,onesided_FE_Mon_BF(1,:),'LineWidth',2); hold on
% title('Frequency Response of Body (No Filt)')
% xlabel('Frequency (Hz)','Fontname','Timesnewroman')
% xlim([0 25])
% plotformat
% subplot(2,1,2)
% plot(f,onesided_FE_Mon_BF(2,:),'LineWidth',2); hold on
% title('Frequency Response of Front Wheel (No Filt)')
% xlabel('Frequency (Hz)','Fontname','Timesnewroman')
% xlim([0 25])
% plotformat


% Storage matrices and cell arrays (used to store information for machine learning)
Monitor_Vehicle_Time{Day,ii}=T_Mon;
Monitor_Vehicle_Acceleration{Day,ii}=za_u_Mon;
Monitor_Vehicle_Frequency_Amp_Data.Original{Day,ii}=onesided_FE_Mon_Original(:, FFT_SELECT_START:FFT_SELECT_END);
%Monitor_Vehicle_Frequency_Amp_Data.Filtered{Day,ii}=onesided_FE_Mon_BF;
Monitor_Vehicle_Frequency_Data{Day,ii}=f(FFT_SELECT_START:FFT_SELECT_END);

Other_Vehicle_Time{Day,ii}=T_Sec;
Other_Vehicle_Acceleration{Day,ii}=za_u_Sec;
Other_Vehicle_Frequency_Amp_Data.Original{Day,ii}=onesided_FE_Sec_Original(:, FFT_SELECT_START:FFT_SELECT_END);
%Other_Vehicle_Frequency_Amp_Data.Filtered{Day,ii}=onesided_FE_Sec_BF;
Other_Vehicle_Frequency_Data{Day,ii}=f(FFT_SELECT_START:FFT_SELECT_END);

%% End of second section of "If Multiple_Vehicle" loop

%% Beginning of last section of "If Multiple_Vehicle" loop (This section is for if there is only one vehicle ever being considered)
else
    % Vehicle parameters
mv_Mon=normrnd(VehicleVariables(row(Day,ii,1),1),VehicleVariables(row(Day,ii,1),2));% sprung mass of vehicle kg (Randomly selects 1 of 8 vehicles)
mw_Mon=normrnd(VehicleVariables(row(Day,ii,1),3),VehicleVariables(row(Day,ii,1),4)); % wheel mass of vehicle kg
kv_Mon=VehicleVariables(row(Day,ii,1),5); %Stiffness of vehicle spring N/m
kw_Mon=VehicleVariables(row(Day,ii,1),6); %Stiffness of vehicle tire N/m
cv_Mon=VehicleVariables(row(Day,ii,1),7); % Damping of vehicle spring N*s/m
cw_Mon=VehicleVariables(row(Day,ii,1),8); % Damping of vehicle tire N*s/m

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
muu=diag([mv_Mon mw_Mon]); muw=[0;0];  mwu=[0,0];  mww=0;
M_Mon=[muu,muw;mwu,mww]; % Vehicle mass matrix
kuu=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
kuw=[0;-kw_Mon]; kwu=[0,-kw_Mon]; kww=kw_Mon;
K_Mon=[kuu,kuw;kwu,kww]; % Vehicle stiffness matrix
cuu=[cv_Mon, -cv_Mon; -cv_Mon, cv_Mon+cw_Mon];
cuw=[0;-cw_Mon]; cwu=[0,-cw_Mon]; cww=cw_Mon;
fue_t_dt=[0;0]; fwe_t_dt=-9.81*mv_Mon-9.81*mw_Mon;
lw=1;

ei_Mon=eig(kuu,muu); % eigenvalues
ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s]
Monitorfv{Day,ii}=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass

% Initial Condition Bridge
UB(:,1)=zeros(2*NumberElements,1); % Initial global displacements
VB(:,1)=zeros(2*NumberElements,1); % Initial global velocities
AB(:,1)=zeros(2*NumberElements,1); % Initial global accelerations
UB1=zeros(2*(NumberElements+1),Kt);
AB1=zeros(2*(NumberElements+1),Kt);
VB1=zeros(2*(NumberElements+1),Kt);

% Initial Condition Vehicle
zu=zeros(2,Kt-1); % Initial displacement upper Vehicle (vertical, rotational, wheel 1, wheel 2)
zw=zeros(1,Kt-1); % Initial displacement lower Vehicle (veritcal at both tires)
zv_u=zeros(2,Kt-1);  % Initial velocity of upper Vehicle
zv_w=zeros(1,Kt-1);  % Initial velocity of lower Vehicle
za_u=zeros(2,Kt-1);  % Initial acceleration of upper Vehicle
za_w=zeros(1,Kt-1);  % Initial acceleration of lower Vehicle

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
mc=lw\(mww-PSIwu/PSIuu*muw); %Mass contact matrix
cc=lw\(cww-PSIwu/PSIuu*cuw); %Damping contact matrix
kc=lw\(kww-PSIwu/PSIuu*kuw); %Stiffness contact matrix

%% Pre-allocate matrix sizes
rc=zeros(2*(NumberElements+1),1);

%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
cor_Mon(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1 % This loop is for calculating the accleration response for each vehicle crossing

xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
Nc=[1-3*xb.^2+2*xb.^3, xc.*(1-2*xb+xb.^2), 3*xb.^2-2*xb.^3, xc.*(xb.^2-xb)];
Ncd=[-6*(xc/l^2)+6*(xc.^2/l^3), 1-4*xb+3*xb.^2, 6*(xc/l^2)-6*(xc.^2/l^3), 3*xb.^2-2*xb];
Ncdd=[-6/l^2+12*xc/l^3, -4/l+6*xc/l^2, 6/l^2-12*xc/l^3, 6*xc/l^2-2/l];
Nct=transpose(Nc); % Column Vector

% Surface Rougness
if Surface==1
difference = abs(RoadMatrix(1,:)-xg);
[~, Column] = min(difference);
rx=RoadMatrix(2,Column(1,1));
drx=RoadMatrix(3,Column(1,1));
ddrx=RoadMatrix(4,Column(1,1));
rc_tdt_st=Nct*(V(Day,ii,1)^2*mc*ddrx+V(Day,ii,1)*cc*drx+kc*rx); %Contact force from road irregularities
RC=rc;
rch=KInsert2(rc_tdt_st,cor_Mon(j,:),2*(NumberElements+1));
RC=RC+rch;
end

% Load Vectors
qu_t=muu*(a1*zv_u(:,i)+a2*za_u(:,i))+cuu*(a6*zv_u(:,i)+a7*za_u(:,i))-kuu*zu(:,i);
qw_t=mwu*(a1*zv_u(:,i)+a2*za_u(:,i))+cwu*(a6*zv_u(:,i)+a7*za_u(:,i))-kwu*zu(:,i);
pc_tdt=lw\(PSIwu/PSIuu*fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu/PSIuu*qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct*mc*Nc; %Mass
cc_st=Nct*(2*V(Day,ii,1)*mc*Ncd+cc*Nc); %Damping
kc_st=Nct*(V(Day,ii,1)^2*mc*Ncdd+V(Day,ii,1)*cc*Ncd+kc*Nc); %Stiffness

% Equivalent nodal loads
pc_tdt_st=Nct*pc_tdt;
qc_t_st=Nct*qc_t;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh=KInsert(kc_st,cor_Mon(j,:),2*(NumberElements+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor_Mon(j,:),2*(NumberElements+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

CB1=CB; % resets global matrix each time
cch=KInsert(cc_st,cor_Mon(j,:),2*(NumberElements+1));
CB1=CB1+cch; % Updated beam stiffness matrix

% Apply boundary conditions
[M,K,C]=boundarycondition(MB1,KB1,NumberElements,CB1);

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

if Surface==1
RS=-PC-QC-RC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
else
RS=-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);
end

% Finding unknown displacements
du=LS\RS;

% Future Disp, Vel, Acc
UB(:,i+1)=UB(:,i)+du; % Future displacement
AB(:,i+1)=a0*du-a1*VB(:,i)-a2*AB(:,i); % Future Acceleration
VB(:,i+1)=VB(:,i)+a3*AB(:,i)+a4*AB(:,i+1); % Future Velocity

% Add back constrained dof
UB1(2:2*NumberElements,i)=UB(1:2*NumberElements-1,i+1);
UB1(2*(NumberElements+1),i)=UB(2*NumberElements,i+1);
AB1(2:2*NumberElements,i)=AB(1:2*NumberElements-1,i+1);
AB1(2*(NumberElements+1),i)=AB(2*NumberElements,i+1);
VB1(2:2*NumberElements,i)=VB(1:2*NumberElements-1,i+1);
VB1(2*(NumberElements+1),i)=VB(2*NumberElements,i+1);

% Element displacement, velocity, acceleration
db=UB1(cor_Mon(j,:),i); % Local Displacement
vb=VB1(cor_Mon(j,:),i); % Local Velocity
ab=AB1(cor_Mon(j,:),i); % Local Acceleration

% Contact points
dc=Nc*db; % Contact displacement
vc=V(Day,ii,1)*Ncd*db+Nc*vb; % Contact velocity
ac=Nc*ab+2*V(Day,ii,1)*Ncd*vb+(V(Day,ii,1)^2)*Ncdd*db; % Contact acceleration

if Surface==1
% Future displacement, velocity and acceleration in lower vehicle
zw(:,i+1)=dc+rx; % Future displacement
zv_w(:,i+1)=vc+V(Day,ii,1)*drx; % Future Velocity
za_w(:,i+1)=ac+V(Day,ii,1)^2*ddrx; % Future Acceleration
else
% Future displacement, velocity and acceleration in lower vehicle
zw(:,i+1)=dc; % Future displacement
za_w(:,i+1)=ac; % Future Acceleration
zv_w(:,i+1)=vc; % Future Velocity
end

% Contact Force
quc_tdt=muw*za_w(:,i+1)+cuw*zv_w(:,i+1)+kuw*zw(:,i+1);

% Change in upper vehicle displacement
dz=PSIuu\(fue_t_dt-quc_tdt+qu_t);

% Future displacement, vnocity and accneration in upper vehicle
zu(:,i+1)=zu(:,i)+dz; % Future displacement
za_u(:,i+1)=a0*dz-a1*zv_u(:,i)-a2*za_u(:,i); % Future Acceleration
zv_u(:,i+1)=zv_u(:,i)+a3*za_u(:,i)+a4*za_u(:,i+1); % Future Velocity

if xc>=l
j=j+1;
J=J+1;
xo=nodes(J-1,2); % Updated for new element
cor_Mon(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V(Day,ii,1);
end % end single vehicle loop

% Shifting acceleration data out of time domain and into frequency domain
Fs = 1/dT; % Sampling frequency
nFFT = Fs/FFT_STEP;
if rem(nFFT,2)>0
    nFFT = nFFT+1;
end
f = Fs*(0:(nFFT/2))/(nFFT); % Frequency domain

za_u = lowpass(za_u',100,Fs,'Steepness',.99)';

% Executing FFt for Upper Vehicle
fft_V=abs(fft(za_u,nFFT,2)/nFFT);
onesided_FE_Original = fft_V(:,1:nFFT/2+1); % Single-sided spectrum
onesided_FE_Original(:,2:end-1) = 2*onesided_FE_Original(:,2:end-1); % Scale Power by 2

% Apply Filters
% Bandpass
%fpass=[4 9];
%Filt_BF = bandpass(za_u',fpass,Fs,'Steepness',.99);
%fft_V = abs(fft(Filt_BF',nFFT,2)/nFFT);
%onesided_FE_BF = fft_V(:,1:nFFT/2+1); % Single-sided spectrum
%onesided_FE_BF(:,2:end-1) = 2*onesided_FE_BF(:,2:end-1); % Scale Power by 2

% % Storage matrices and cell arrays (used to store information for machine
% % learning)
Monitor_Vehicle_Time{Day,ii}=T;
Monitor_Vehicle_Acceleration{Day,ii}=za_u;
Monitor_Vehicle_Frequency_Amp_Data.Original{Day,ii}=onesided_FE_Original(:, FFT_SELECT_START:FFT_SELECT_END);
%Monitor_Vehicle_Frequency_Amp_Data.Filtered{Day,ii}=onesided_FE_BF;
Monitor_Vehicle_Frequency_Data{Day,ii}=f(FFT_SELECT_START:FFT_SELECT_END);

end

end % end day loop

save(outPath)
exitCode = 0;
end
