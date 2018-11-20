clc; clear all;
%Bridge Parameters
L=25; % Length m
n=20; % Number of elements bridge is divided into
l=L/n; % Length of each individual element
CN=(n+2); % Central node of bridge
A=2; % Cross-sectional area m^2
E=35e+9; % Modulus of elasticity of bridge N/m^2
I=.12; % Moment of Inertia m^4
mu=28125; % mass per unit length kg/m
beta=0.05; % Damping of bridge
kb=KBeam(E,I,l); % Stiffness matrix for bridge
% mb=MBeam(mu,l); % Consistent mass matrix for bridge
mb=(mu*l/2)*[1 0 0 0;0 1/12*(l)^2 0 0;0 0 1 0; 0 0 0 1/12*(l)^2]; % Lumped mass matrix
wn=sqrt((pi^4/L^4)*(E*I/mu)) % First natural frequency
w2=sqrt((2^4*pi^4/L^4)*(E*I/mu)); % Second natural frequency
[ele, nodes]=element(n,l,L);
% StartDamage=randsample(nodes(1:(n-1),2),2); % Where damage locations begin
% EndDamage=StartDamage+((l*.5*l)*rand(2,1)+.5*l); % Where damage locations ends
x=L/2; % Location of center of bridge

% Vehicle parameters
VehicleVariables=load('VehicleData.dat');
V=10; % Speed of vehicle m/s
mv=1200;% sprung mass of vehicle kg (Randomly selects 1 of 3 vehicles)
mw=0; % wheel mass of vehicle kg
P=-9.81*(mv+mw); % Point load of vehicle N
kv=500000; %Stiffness of vehicle spring N/m
kw=0;%VehicleVariables(row,6); %Stiffness of vehicle tire N/m
cs=0;%VehicleVariables(row,4); % Damping of vehicle spring N*s/m
cw=0;%VehicleVariables(row,5); % Damping of vehicle tire N*s/m
w=pi*V/L; % circular frequency
wv=sqrt(kv/mv); % Natural frequency of sprung mass
wb=wn*beta; % circular frequency of damping
fr=wn/wv; %Frequency ratio
sf1=wn-pi*V/L; %Shifted frequency 1
sf2=wn+pi*V/L; %Shifted frequency 2
omg=2*pi*V/L; % Circular (driving) frequency of vehicle (force) rad/s
f1=wn/(2*pi); % 1st natural frequency of bridge (hz)
fv=sqrt(kv/mv)/(2*pi); % Natural frequency of sprung mass
cr=2*f1*L; % Critical speed of bridge m/s
S=w/wn; % Speed parameter for bridge
vo=2*P*L^3/(pi^4*E*I); % Static displacement at midspan

% Vehicle Acceleration Components
A1=1-(1/(1-(2*fr*S)^2))-(S/(1-fr^2*(1-S)^2))+(S/(1-fr^2*(1+S)^2));
A2=(2*fr*S)^2/(1-(2*fr*S)^2);
A3=(S*fr^2*(1-S)^2)/(1-fr^2*(1-S)^2);
A4=(S*fr^2*(1+S)^2)/(1-fr^2*(1+S)^2);

% Time and position
dT=.001; % Time Step
T=0:dT:L/V; % Total time to cross bridge
Kt=length(T); % Number of Data points
xg=0; % Initial global position
j=1; % Initial row for elemental matrix
J=2; % Initial row for nodal matrix

%% Matricies and Initial Conditions
% Beam Matricies
KB=zeros(2*(n+1),2*(n+1));
MB=zeros(2*(n+1),2*(n+1));
for i=1:n
   cor(i,:)=ele(i,2:5);
   kk(:,:,i)=KInsert(kb,cor(i,:),2*(n+1));
   KB=KB+kk(:,:,i); % Beam stiffness matrix
   mm(:,:,i)=KInsert(mb,cor(i,:),2*(n+1));
   MB=MB+mm(:,:,i); % Beam mass matrix
end

% Remove boundary conditions only so we can calculate damping matirx 
M=zeros(2*n,2*n);
K=zeros(2*n,2*n);
[M,K]=boundarycondition(MB,KB,M,K,n);
ei=eig(K,M); % eigenvalues
ef=sort(real(sqrt(ei))); % sorted natural angular frequencies [rad/s] 
wn_FEA=ef; % sorted natural angular frequencies

% Beam damping matrix
al=2*beta*wn_FEA(1,1)*wn_FEA(2,1)/(wn_FEA(1,1)+wn_FEA(2,1)); % Alpha for Rayleigh Damping
be=2*beta/(wn_FEA(1,1)+wn_FEA(2,1)); % Beta for Rayleigh Damping
CB=al*MB+be*KB; % Damping Matrix for beam

% Load Matricies
pc=zeros(2*(n+1),1);
qc=zeros(2*(n+1),1);

% Vehicle matricies
muu=mv; muw=0;  mwu=0;  mww=mw;
Mv=[muu,muw;mwu,mww]; % Vehicle mass matrix
kuu=kv; kuw=-kv;    kwu=-kv;    kww=kv+kw;
Kv=[kuu,kuw;kwu,kww]; % Vehicle stiffness matrix
cuu=cs; cuw=-cs;    cwu=-cs;    cww=cs+cw;
Cv=[cuu,cuw;cwu,cww]; % Vehicle damping matrix
fue_t_dt=0; fwe_t_dt=-9.81*mv-9.81*mw;
Fe=[fue_t_dt;fwe_t_dt]; % External forces
lw=1;

% Initial Condition Bridge
UB(:,1)=zeros(2*n,1); % Initial global displacements
VB(:,1)=zeros(2*n,1); % Initial global velocities
AB(:,1)=zeros(2*n,1); % Initial global accelerations
UB1=zeros(2*(n+1),Kt+1);
AB1=zeros(2*(n+1),Kt+1);
VB1=zeros(2*(n+1),Kt+1);

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

%% Pre-allocate matrix Sizes
dz=zeros(1,Kt-1); % Change in upper vehicle displacement
quc_tdt=zeros(1,Kt);% Uppdated contact Force
ac=zeros(1,Kt-1);  % Contact acceleration
dc=zeros(1,Kt-1);  % Contact displacement
vc=zeros(1,Kt-1);  % Contact velocity
ab=zeros(4,Kt-1);  % Local element acceleration
db=zeros(4,Kt-1);  % Local element displacement
vb=zeros(4,Kt-1);  % Local element velocity
du=zeros(n*2,Kt-1); % Change in displacement of global matrices
Nc=zeros(Kt-1,4);   % Shape function
Ncd=zeros(Kt-1,4);  % Derivative of shape function
Ncdd=zeros(Kt-1,4); % Double derivative of shape function
wn_FEA=zeros(2*n,Kt-1);

%% Calc global position and shape functions
xo=nodes(J-1,2); % element start location
cor(j,:)=ele(j,2:5); % Current element coordinate

for i=1:Kt-1

xc=(xg-xo); % Local position
xb=xc/l; % Local coordiante
t=xc/V; % Local time
Nc(i,:)=[1-3*xb^2+2*xb^3, xc*(1-2*xb+xb^2), 3*xb^2-2*xb^3, xc*(xb^2-xb)]; % Shape function Row Vector
Ncd(i,:)=[-6*V^2*t/l^2+6*V^3*t^2/l^3, V-4*V^2*t/l+3*V^3*t^2/l^2, 6*V^2*t/l^2-6*V^3*t^2/l^3, 3*V^3*t^2/l^2-2*V^2*t/l];
Ncdd(i,:)=[-6*V^2/l^2+12*V^3*t/l^3, -4*V^2/l+6*V^3*t/l^2, 6*V^2/l^2-12*V^3*t/l^3, 6*V^3*t/l^2-2*V^2/l];
Nct=transpose(Nc); % Column Vector

% Load Vectors    
qu_t=muu*(a1*zv(1,i)+a2*za(1,i))+cuu*(a6*zv(1,i)+a7*za(1,i))-kuu*zu(1,i);
qw_t=mwu*(a1*zv(1,i)+a2*za(1,i))+cwu*(a6*zv(1,i)+a7*za(1,i))-kwu*zu(1,i);
pc_tdt=lw\(PSIwu*PSIuu\fue_t_dt-fwe_t_dt);
qc_t=lw\(PSIwu*PSIuu\qu_t-qw_t);

% Interaction contact matrices
mc_st=Nct(:,i)*mc*Nc(i,:); %Mass
cc_st=Nct(:,i)*(2*V*mc*Ncd(i,:)+cc*Nc(i,:)); %Damping
kc_st=Nct(:,i)*(V^2*mc*Ncdd(i,:)+V*cc*Ncd(i,:)+kc*Nc(i,:)); %Stiffness

% Equivalent nodal loads
pc_tdt_st=Nct(:,i)*pc_tdt;
qc_t_st=Nct(:,i)*qc_t;

% Global Stiffness
KB1=KB; % resets global matrix each time
kkh=KInsert(kc_st,cor(j,:),2*(n+1));
KB1=KB1+kkh; % Updated beam stiffness matrix

% Global Mass
MB1=MB; % resets global matrix each time
mmh=KInsert(mc_st,cor(j,:),2*(n+1));
MB1=MB1+mmh; % Updated beam stiffness matrix

% Global Damping
CB1=CB; % resets global matrix each time
cch=KInsert(cc_st,cor(j,:),2*(n+1));
CB1=CB1+cch; % Updated beam stiffness matrix

% Global Contact Loads
PC1=pc; 
pch=KInsert2(pc_tdt_st,cor(j,:),2*(n+1));
PC1=PC1+pch;

QC1=qc;
qch=KInsert2(qc_t_st,cor(j,:),2*(n+1));
QC1=QC1+qch;

% Apply Boundary Conditions for Beam Matrices
% Delete first and next to last row/column
C=zeros(2*n,2*n);
[M,K,C]=boundarycondition(MB1,KB1,M,K,n,CB1,C);

% Apply Boundary Conditions for Load Matrices
% Delete first and next to last row
PC=zeros(2*n,1);
QC=zeros(2*n,1);
[PC,QC]=bcloads(n,PC1,PC,QC1,QC);

% Apply Numark Method
RS=-PC-QC-K*UB(:,i)-(-a1*M+C-a1*a4*C)*VB(:,i)-(-a2*M+a3*C-a4*a2*C)*AB(:,i);
LS=(a0*M+a4*a0*C+K);

% Finding unknown displacements
du(:,i)=LS\RS;

% Future Disp, Vel, Acc
UB(:,i+1)=UB(:,i)+du(:,i); % Future displacement
AB(:,i+1)=a0*du(:,i)-a1*VB(:,i)-a2*AB(:,i); % Future Acceleration
VB(:,i+1)=VB(:,i)+a3*AB(:,i)+a4*AB(:,i+1); % Future Velocity


% Add back constrained dof
UB1(2:2*n,i)=UB(1:2*n-1,i);
UB1(2*(n+1),i)=UB(2*n,i);
AB1(2:2*n,i)=AB(1:2*n-1,i);
AB1(2*(n+1),i)=AB(2*n,i);
VB1(2:2*n,i)=VB(1:2*n-1,i);
VB1(2*(n+1),i)=VB(2*n,i);

% Element displacement, velocity, acceleration
db(:,i)=UB1(cor(j,:),i); % Local Displacement
vb(:,i)=VB1(cor(j,:),i); % Local Velocity
ab(:,i)=AB1(cor(j,:),i); % Local Acceleration

% Contact points
dc(:,i)=Nc(i,:)*db(:,i); % Contact displacement
vc(:,i)=V*Ncd(i,:)*db(:,i)+Nc(i,:)*vb(:,i); % Contact velocity
ac(:,i)=Nc(i,:)*ab(:,i)+2*V*Ncd(i,:)*vb(:,i)+(V^2)*Ncdd(i,:)*db(:,i); % Contact acceleration

% Future displacement, velocity and acceleration in lower vehicle
zu(2,i+1)=dc(:,i); % Future displacement
za(2,i+1)=vc(:,i); % Future Acceleration
zv(2,i+1)=ac(:,i); % Future Velocity

% Contact Force
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
cor(j,:)=ele(j,2:5); % New element coordinates
end
% Update global x position for next loop
xg=xg+dT*V;
end

dispV=zeros(1,Kt);
accV1=zeros(1,Kt);
accV2=zeros(1,Kt);
dispB1=zeros(1,Kt);
for i=1:Kt; % Time for vehicle to enter and exit bridge
    if beta > 0;
        J=1;
        wb=beta*wn;
        % Displacement of Bridge with Bridge Damping
       dispB1(i)=vo*sum(sin(J*pi*x/L)*1/(J^2*(J^2-S^2))*(sin(J*w*T(i))-(S/J)*exp(-wb*T(i))*sin(wn*T(i))));
    else
        % Displacement of Bridge with no bridge damping
        dispB1(i)=(vo/(1-S^2))*(sin(w*T(i))-S*sin(wn*T(i))); % Equation from Yang 2004 
    end
    % Acceleration of Bridge with no damping
    accB1(i)=(vo/(1-S^2))*sin(pi/2)*(-w^2*sin(w*T(i))+w*wn*sin(wn*T(i)));
    % Displacement of Vehicle with no vehicle damping
    dispV(i)=(vo/(2*(1-S^2)))*((1-cos(wv*T(i)))-((cos(omg*T(i))-cos(wv*T(i)))/(1-(2*fr*S)^2))-(S*(cos(sf1*T(i))-cos(wv*T(i)))/(1-fr^2*(1-S)^2))+(S*(cos(sf2*T(i))-cos(wv*T(i)))/(1-fr^2*(1+S)^2)));
    
    % Acceleration of Vehicle
    accV1(i)=((vo*wv^2)/(2*(1-S^2)))*(cos(wv*T(i))+((2*fr*S)^2*cos(omg*T(i))-cos(wv*T(i)))/(1-(2*fr*S)^2)+(S*(fr^2*(1-S)^2*cos(sf1*T(i))-cos(wv*T(i)))/(1-fr^2*(1-S)^2))-(S*(fr^2*(1+S)^2*cos(sf2*T(i))-cos(wv*T(i)))/(1-fr^2*(1+S)^2)));
%     accveh(i)=((vo*wv^2)/(2*(1-S^2)))*A1*cos(wv*T(i)); % Acceleration related to vehicle frequency
%     accspeed(i)=((vo*wv^2)/(2*(1-S^2)))*A2*cos(omg*T(i)); % Accleration related to vehicle speed
%     accbri(i)=((vo*wv^2)/(2*(1-S^2)))*(A3*cos(sf1*T(i))-A4*cos(sf2*T(i))); % Acceleration associated with bridge frequency
%     accV2(i)=accveh(i)+accspeed(i)+accbri(i);
end


%% Plots
% figure(1)
% set(gcf,'color','white')
% subplot(2,1,1);
% plot(0:dT:L/V,UB(CN,:),'b','linewidth',3);hold on
% plot(0:dT:L/V,dispB1,'-.r','linewidth',3); 
% title('Displacement of Bridge')
% xlabel('Time (s) ');
% ylabel('Displacement (m)');
% legend('FEM','Analytical','location','southeast','FontName','Timesnewroman')
% plotformat
% subplot(2,1,2);
% plot(0:dT:L/V,zu(1,:),'b','linewidth',3);hold on
% plot(0:dT:L/V,dispV,'-.r','linewidth',3); 
% title('Displacement of Vehicle')
% xlabel('Time (s) ');
% ylabel('Displacement (m)');
% legend('FEM','Analytical','location','southeast','FontName','Timesnewroman')
% plotformat

figure(2)
set(gcf,'color','white')
% subplot(2,1,1)
% plot(0:dT:L/V,AB(CN,:),'b', 'LineWidth',3);hold on
% plot(0:dT:L/V,accB1,'-.r', 'LineWidth',3); 
% set(gca,'fontsize',16);
% xlabel('Time (s) ');
% ylabel('Acceleration (m/s^2)');
% legend('FEM','Analytical','location','southeast','FontName','Timesnewroman')
% title('Acceleration of Bridge')
plotformat
plot(0:dT:L/V, za(1,:),'b','LineWidth',3);hold on
% plot(0:dT:L/V, accV1,'-.r','LineWidth',3); 
set(gca,'fontsize',16);
xlabel('Time (s) ');
ylabel('Acceleration (m/s^2)');
% legend('FEM','Analytical','location','southeast')
% title('Acceleration of Vehicle')
plotformat

%% Fast Fourier Transform
Fs = 1/dT; % Sampling frequency                          
t = (0:Kt-1)*dT; % Time vector
f = Fs*(0:(Kt))/(Kt*2); % Frequency domain

% Executing FFT for Vehicle
fftV_FE=abs(fft(za(1,:),2*Kt));
fftV_Ana=abs(fft(accV1,2*Kt));

Twosided_FE = fftV_FE/Kt; % two-sided spectrum
Twosided_Ana= fftV_Ana/Kt; 
onesided_FE = Twosided_FE(1:Kt+1); % Single-sided spectrum
onesided_Ana=Twosided_Ana(1:Kt+1);
onesided_FE(2:end-1) = onesided_FE(2:end-1);
onesided_Ana(2:end-1) = onesided_Ana(2:end-1);
% % 
figure(3)
set(gcf,'color','white')
% [pks_Ana,locs_Ana] = findpeaks(onesided_Ana);
plot(f,onesided_FE,'color','b','LineWidth',2);hold on
% [pks_FE,locs_FE] = findpeaks(onesided_FE);
% semilogx(f,onesided_Ana,'-.r','LineWidth',2);
axis([0 10 0 .016]);
set(gca,'fontsize',16);
% title('Power Spectrum of Vehicle Acceleration','Fontname','Timesnewroman')
% legend('FEM','Analytical','location','northeast','FontName','Timesnewroman')
xlabel('Frequency (Hz)','Fontname','Timesnewroman')
ylabel('Spectrum Amplitude','Fontname','Timesnewroman')
plotformat
% 
% % Executing FFt for Bridge
% fftB_FE=fft(AB(CN,:),2*Kt);
% fftB_Ana=fft(accB1,2*Kt);
% 
% Twosided_FE = abs(fftB_FE)/Kt; % two-sided spectrum
% Twosided_Ana=abs(fftB_Ana)/Kt; 
% onesided_FE = Twosided_FE(1:Kt+1); % Single-sided spectrum
% onesided_Ana=Twosided_Ana(1:Kt+1);
% onesided_FE(2:end-1) = onesided_FE(2:end-1);
% onesided_Ana(2:end-1) = onesided_Ana(2:end-1);
% 
% %%
% figure(4)
% semilogx(f(1:100),onesided_FE(1:100),'color','b','LineWidth',2);hold on
% semilogx(f(1:100),onesided_Ana(1:100),'-.r','LineWidth',2); 
% % axis([0 100 0 .01]);
% title('Power Spectrum of Bridge Acceleration')
% legend('FEM','Analytical','location','northeast','FontName','Timesnewroman')
% xlabel('Frequency (Hz)','Fontname','Timesnewroman')
% ylabel('Spectrum Amplitude','Fontname','Timesnewroman')
% plotformat
% % 
% % 
% % % Obtaining Bridge Frequency
% % [PeakValue_FE, maxIndex_FE] = max(onesided_FE);
% % [PeakValue_Ana, maxIndex_Ana] = max(onesided_Ana);
% % bridgefrequency_FE=f(maxIndex_FE)
% % bridgefrequency_Ana=f(maxIndex_Ana);
% % % line([bridgefrequency_FE,bridgefrequency_FE],[0,PeakValue_FE],'Linestyle','--','color','k');
% % % 
% % % %Adjust Bridge Frequency
% % % FE_wn=bridgefrequency_FE+(w/(2*pi));
% % % Ana_en=bridgefrequency_Ana+(w/(2*pi));
