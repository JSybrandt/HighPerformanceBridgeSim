function [Mon_Vehicle,Oth_Vehicle_1,Oth_Vehicle_2,Kt]=Vehicle_Parameters(number_vehicles,L,dT)
global xg j J zu_Monitor zv_Monitor za_Monitor zu_Oth_Vehicle_1 zv_Oth_Vehicle_1 za_Oth_Vehicle_1 zu_Oth_Vehicle_2 zv_Oth_Vehicle_2 za_Oth_Vehicle_2
VehicleVariables=load('VehicleData.dat');

if number_vehicles==1
V1=10; % Speed of vehicle m/s
mv1=1200;% sprung mass of vehicle kg, VehicleVariables column1
mw1=0; % wheel mass of vehicle kg, VehicleVariables column2
kv1=500000; %Stiffness of vehicle spring N/m, VehicleVariables column3
kw1=0;%VehicleVariables column6; %Stiffness of vehicle tire N/m
cs1=0;%VehicleVariables column4; % Damping of vehicle spring N*s/m
cw1=0;%VehicleVariables column5; % Damping of vehicle tire N*s/mm

muu1=mv1; muw1=0;    mwu1=0;    mww1=mw1;
kuu1=kv1; kuw1=-kv1; kwu1=-kv1; kww1=kv1+kw1;
cuu1=cs1; cuw1=-cs1; cwu1=-cs1; cww1=cs1+cw1;
fue_t_dt1=0;         fwe_t_dt1=-9.81*mv1-9.81*mw1;
lw1=1;

Mon_Vehicle=[V1,fue_t_dt1,fwe_t_dt1,lw1;
             muu1,muw1,mwu1,mww1;
             kuu1,kuw1,kwu1,kww1;
             cuu1,cuw1,cwu1,cww1];
Oth_Vehicle_1=[];
Oth_Vehicle_2=[];

T1=0:dT:L/V1; % Total time to cross bridge
Kt=length(T1);% Number of Data points (Mon_Vehicle)

xg=0; % Initial global position (Mon_Vehicle)
j=1; % Initial row for elemental matrix (Mon_Vehicle)
J=2; % Initial row for nodal matrix (Mon_Vehicle)

% Initial Condition Vehicle
zu_Monitor(:,1)=[0;0]; % Initial displacement vehicle
zv_Monitor(:,1)=[0;0]; % Initial velocity of vehicle
za_Monitor(:,1)=[0;0]; % Initial acceleration of vehicle

% Contact Matricies
PSIuu=a0*muu1+a5*cuu1+kuu1;
PSIwu=a0*mwu1+a5*cwu1+kwu1;
mc_Monitor=lw1\(mww1-PSIwu*PSIuu\muw1); %Mass contact matrix
cc_Monitor=lw1\(cww1-PSIwu*PSIuu\cuw1); %Damping contact matrix
kc_Monitor=lw1\(kww1-PSIwu*PSIuu\kuw1); %Stiffness contact matrix

elseif number_vehicles==2
V1=10; % Speed of vehicle m/s
V2=10;
mv1=1200;% sprung mass of vehicle kg, VehicleVariables column1
mv2=1000;
mw1=0; % wheel mass of vehicle kg, VehicleVariables column2
mw2=0;
kv1=500000; %Stiffness of vehicle spring N/m, VehicleVariables column3
kv2=550000;
kw1=0;%VehicleVariables column6; %Stiffness of vehicle tire N/m
kw2=0;
cs1=0;%VehicleVariables column4; % Damping of vehicle spring N*s/m
cs2=0;
cw1=0;%VehicleVariables column5; % Damping of vehicle tire N*s/mm
cw2=0;

muu1=mv1; muw1=0;    mwu1=0;    mww1=mw1;
kuu1=kv1; kuw1=-kv1; kwu1=-kv1; kww1=kv1+kw1;
cuu1=cs1; cuw1=-cs1; cwu1=-cs1; cww1=cs1+cw1;
fue_t_dt1=0;         fwe_t_dt1=-9.81*mv1-9.81*mw1;
lw1=1;

muu2=mv2; muw2=0;    mwu2=0;    mww2=mw2;
kuu2=kv2; kuw2=-kv2; kwu2=-kv2; kww2=kv2+kw2;
cuu2=cs2; cuw2=-cs2; cwu2=-cs2; cww2=cs2+cw2;
fue_t_dt2=0;         fwe_t_dt2=-9.81*mv2-9.81*mw2;
lw2=1;

Mon_Vehicle=[V1,fue_t_dt1,fwe_t_dt1,lw1;
             muu1,muw1,mwu1,mww1;
             kuu1,kuw1,kwu1,kww1;
             cuu1,cuw1,cwu1,cww1];
Oth_Vehicle_1=[V2,fue_t_dt2,fwe_t_dt2,lw2;
             muu2,muw2,mwu2,mww2;
             kuu2,kuw2,kwu2,kww2;
             cuu2,cuw2,cwu2,cww2];
Oth_Vehicle_2=[];

T1=0:dT:L/V1; % Total time to cross bridge
T2=0:dT:L/V2;
Kt=round(length(T1)/2)+length(T2);% Number of Data points (Mon_Vehicle,Oth_Vehicle_1)

xg=[0,L]; % Initial global position (Mon_Vehicle,Oth_Vehicle_1)
j=[1,1]; % Initial row for elemental matrix (Mon_Vehicle,Oth_Vehicle_1)
J=[2,2]; % Initial row for nodal matrix (Mon_Vehicle,Oth_Vehicle_1)

% Initial Condition Vehicle
zu_Monitor(:,1)=[0;0]; % Initial displacement vehicle
zv_Monitor(:,1)=[0;0]; % Initial velocity of vehicle
za_Monitor(:,1)=[0;0]; % Initial acceleration of vehicle
zu_Oth_Vehicle_1(:,1)=[0;0]; % Initial displacement vehicle
zv_Oth_Vehicle_1(:,1)=[0;0]; % Initial velocity of vehicle
za_Oth_Vehicle_1(:,1)=[0;0]; % Initial acceleration of vehicle

elseif number_vehicle2==3
V1=10; % Speed of vehicle m/s
V2=10;
V3=10;
mv1=1200;% sprung mass of vehicle kg, VehicleVariables column1
mv2=1000;
mv3=1400;
mw1=0; % wheel mass of vehicle kg, VehicleVariables column2
mw2=0;
mw3=0;
kv1=500000; %Stiffness of vehicle spring N/m, VehicleVariables column3
kv2=550000;
kv3=450000;
kw1=0;%VehicleVariables column6; %Stiffness of vehicle tire N/m
kw2=0;
kw3=0;
cs1=0;%VehicleVariables column4; % Damping of vehicle spring N*s/m
cs2=0;
cs3=0;
cw1=0;%VehicleVariables column5; % Damping of vehicle tire N*s/mm
cw2=0;
cw3=0;

muu1=mv1; muw1=0;    mwu1=0;    mww1=mw1;
kuu1=kv1; kuw1=-kv1; kwu1=-kv1; kww1=kv1+kw1;
cuu1=cs1; cuw1=-cs1; cwu1=-cs1; cww1=cs1+cw1;
fue_t_dt1=0;         fwe_t_dt1=-9.81*mv1-9.81*mw1;
lw1=1;

muu2=mv2; muw2=0;    mwu2=0;    mww2=mw2;
kuu2=kv2; kuw2=-kv2; kwu2=-kv2; kww2=kv2+kw2;
cuu2=cs2; cuw2=-cs2; cwu2=-cs2; cww2=cs2+cw2;
fue_t_dt2=0;         fwe_t_dt2=-9.81*mv2-9.81*mw2;
lw2=1;

muu3=mv3; muw3=0;    mwu3=0;    mww3=mw3;
kuu3=kv3; kuw3=-kv3; kwu3=-kv3; kww3=kv3+kw3;
cuu3=cs3; cuw3=-cs3; cwu3=-cs3; cww3=cs3+cw3;
fue_t_dt3=0;         fwe_t_dt3=-9.81*mv3-9.81*mw3;
lw3=1;


Mon_Vehicle=[V1,fue_t_dt1,fwe_t_dt1,lw1;
             muu1,muw1,mwu1,mww1;
             kuu1,kuw1,kwu1,kww1;
             cuu1,cuw1,cwu1,cww1];
Oth_Vehicle_1=[V2,fue_t_dt2,fwe_t_dt2,lw2;
             muu2,muw2,mwu2,mww2;
             kuu2,kuw2,kwu2,kww2;
             cuu2,cuw2,cwu2,cww2];
Oth_Vehicle_2=[V3,fue_t_dt3,fwe_t_dt3,lw3;
             muu3,muw3,mwu3,mww3;
             kuu3,kuw3,kwu3,kww3;
             cuu3,cuw3,cwu3,cww3];
         
T1=0:dT:L/V1; % Total time to cross bridge
T2=0:dT:L/V2;
T3=0:dT:L/V3;
Kt=round(length(T1)/3)+round(length(T2)/3)+length(T3);% Number of Data points (Mon_Vehicle,Oth_Vehicle_1,Oth_Vehicle_2)

xg=[0,L,L]; % Initial global position (Mon_Vehicle,Oth_Vehicle_1,Oth_Vehicle_2)
j=[1,1,1]; % Initial row for elemental matrix (Mon_Vehicle,Oth_Vehicle_1,Oth_Vehicle_2)
J=[2,2,2]; % Initial row for nodal matrix (Mon_Vehicle,Oth_Vehicle_1,Oth_Vehicle_2) 

% Initial Condition Vehicle
zu_Monitor(:,1)=[0;0]; % Initial displacement vehicle
zv_Monitor(:,1)=[0;0]; % Initial velocity of vehicle
za_Monitor(:,1)=[0;0]; % Initial acceleration of vehicle
zu_Oth_Vehicle_1(:,1)=[0;0]; % Initial displacement vehicle
zv_Oth_Vehicle_1(:,1)=[0;0]; % Initial velocity of vehicle
za_Oth_Vehicle_1(:,1)=[0;0]; % Initial acceleration of vehicle
zu_Oth_Vehicle_2(:,1)=[0;0]; % Initial displacement vehicle
zv_Oth_Vehicle_2(:,1)=[0;0]; % Initial velocity of vehicle
za_Oth_Vehicle_2(:,1)=[0;0]; % Initial acceleration of vehicle
end
end
    