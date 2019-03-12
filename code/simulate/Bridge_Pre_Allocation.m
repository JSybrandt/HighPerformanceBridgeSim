clc; clear all;
% Pre-Processing (Change values manually)
Multiple_Vehicles=1; % 1 if multiple vehicles are considered, 0 if just 1 vehicle is being considered
Surface=1; % 1 if surface is considered, 0 otherwise
Temp=1; % 1 if temp effects are considered, 0 otherwise
Damage=1; % 1 if damage effects are considered, 0 otherwise
Damage_Case=2; % Determines which damage case being analyzed
Bridge=1; % Indicates which bridge is being tested
lim=361; % Number of days monitoring subject bridge
NumberElements=10; % Number of elements bridge is divided into

% Load variable arrays
BridgeVariables=load('BridgeVariables.dat');
VehicleVariables=load('VehicleData.dat');
load('USW00003870.hourly.mat'); 

%% Pre Analysis Calculations (Sets up environmental and vehicle parameters)
Td=0:.016667:24; %Time of day record was taken (Assumes vehicles cross bridge every min of a day)
n=length(Td); % Number of times a monitoring vehicle crossed the bridge

% Bridge Parameters
L=BridgeVariables(Bridge,1); % Length m
W=BridgeVariables(Bridge,6);% Width m
l=L/NumberElements; % Length of each individual element
CN=(NumberElements+2); % Central node of bridge
E=BridgeVariables(Bridge,3); % Modulus of elasticity of bridge N/m^2
I=BridgeVariables(Bridge,4); % Moment of Inertia m^4
mu(1:(lim-1))=BridgeVariables(Bridge,2); % mass per unit length kg/m
bbeta=BridgeVariables(Bridge,5); % Damping of bridge
Rclass=BridgeVariables(Bridge,7);% Values For ISO Road Roughness Classification, from 3 to 9
[ele, ele2, nodes, nodes2]=element(NumberElements,l,L,Multiple_Vehicles);

% Surface Roughness
if Surface==1
[RoadMatrix]=CompleteSurfaceRoughness(Rclass,L);
end

% Damage Variables
if Damage==1 && Damage_Case==1
DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
ED1=[((.01+rand(1)*(.01-.05))*E),.001*rand(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1
elseif Damage==1 && Damage_Case==2
DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.3-lim*.25)); % The day damage is iniciated on bridge
DayDamage2=round(lim*.35)+round(rand(1)*(lim*.45-lim*.35)); % The day second damage is iniciated on bridge
DayDamage3=round(lim*.55)+round(rand(1)*(lim*.65-lim*.55)); % The day damage is iniciated on bridge
DayDamage4=round(lim*.7)+round(rand(1)*(lim*.75-lim*.7)); % The day damage is iniciated on bridge
DayDamage5=round(lim*.8)+round(rand(1)*(lim*.9-lim*.8)); % The day damage is iniciated on bridge
ED1=((.01+rand(1)*(.02-.01))*E); % Damaged Modulus 1
ED2=((.01+rand(1)*(.02-.01))*E); % Damaged Modulus 2   
ED3=((.01+rand(1)*(.02-.01))*E); % Damaged Modulus 3
ED4=((.01+rand(1)*(.02-.01))*E); % Damaged Modulus 4 
ED5=((.01+rand(1)*(.02-.01))*E); % Damaged Modulus 5  
end


if Temp==1
    Tact=zeros(lim,n); %Actual temperature degrees Celsius
end


% Randomized vehicle selection and speed
if Multiple_Vehicles==1
number_vehicles=randi([1 2],lim,n); % Randomly determines how many vehicles will cross bridge
row=zeros(lim,n,3);
V=zeros(lim,n,3);

MonitorVehicleMass=zeros(lim,n);
MonitorWheelMass=zeros(lim,n);
MonitorSuspensionStiffness=zeros(lim,n);
MonitorWheelStiffness=zeros(lim,n);
MonitorSuspensionDamping=zeros(lim,n);
MonitorWheelDamping=zeros(lim,n);
Monitorfv=cell(lim,n);

SecondVehicleMass=zeros(lim,n);
SecondWheelMass=zeros(lim,n);
SecondSuspensionStiffness=zeros(lim,n);
SecondWheelStiffness=zeros(lim,n);
SecondSuspensionDamping=zeros(lim,n);
SecondWheelDamping=zeros(lim,n);
Secondfv=cell(lim,n);

% Storage matrices and cell arrays (used to store information for machine
% learning)
Monitor_Vehicle_Time=cell(lim,n);
Monitor_Vehicle_Acceleration=cell(lim,n);
Monitor_Vehicle_Frequency_Amp_Data.Original=cell(lim,n);
Monitor_Vehicle_Frequency_Amp_Data.Filtered=cell(lim,n);
Monitor_Vehicle_Frequency_Data=cell(lim,n);

Other_Vehicle_Time=cell(lim,n);
Other_Vehicle_Acceleration=cell(lim,n);
Other_Vehicle_Frequency_Amp_Data.Original=cell(lim,n);
Other_Vehicle_Frequency_Amp_Data.Filtered=cell(lim,n);
Other_Vehicle_Frequency_Data=cell(lim,n);

Start_Time_Following_Vehicle=zeros(lim,n);
Order_of_Vehicles=cell(lim,n);

for gg=1:lim
    for hh=1:n
if number_vehicles(gg,hh)==1
    row(gg,hh,1)=randi([1 8]);
    row(gg,hh,2)=0;
%     row(gg,hh,3)=0;
    
    V(gg,hh,1)=randi([10 25]);
    V(gg,hh,2)=0;
%     V(gg,hh,3)=0;
elseif number_vehicles(gg,hh)==2
    row(gg,hh,1)=randi([1 8]);
    row(gg,hh,2)=randi([1 8]);
%     row(gg,hh,3)=0;
    
    V(gg,hh,1)=randi([10 25]);
    V(gg,hh,2)=randi([10 25]);
%     V(gg,hh,3)=0;
% elseif number_vehicles(gg,hh)==3
%     row(gg,hh,1)=randi([1 10]);
%     row(gg,hh,2)=randi([1 10]);
%     row(gg,hh,3)=randi([1 10]);
%     
%     V(gg,hh,1)=randi([10 25]);
%     V(gg,hh,2)=randi([10 25]);
%     V(gg,hh,3)=randi([10 V(gg,hh,2)]);
end
    end
end
else
    
MonitorVehicleMass=zeros(lim,n);
MonitorWheelMass=zeros(lim,n);
MonitorSuspensionStiffness=zeros(lim,n);
MonitorWheelStiffness=zeros(lim,n);
MonitorSuspensionDamping=zeros(lim,n);
MonitorWheelDamping=zeros(lim,n);
Monitorfv=cell(lim,n);

row=randi([1 8],lim,n); % Randomly selects which vehicle is crossing bridge
V=randi([10 25],lim,n); % Speed of vehicle m/s

% Storage matrices and cell arrays (used to store information for machine
% learning)
Monitor_Vehicle_Time=cell(lim,n);
Monitor_Vehicle_Acceleration=cell(lim,n);
Monitor_Vehicle_Frequency_Amp_Data.Original=cell(lim,n);
Monitor_Vehicle_Frequency_Amp_Data.Filtered=cell(lim,n);
Monitor_Vehicle_Frequency_Data=cell(lim,n);
Monitor_Vehicle_Road_Profile=cell(lim,n);
Monitor_Vehicle_Derivative_Road_Profile=cell(lim,n);
Monitor_Vehicle_Other_Derivative_Road_Profile=cell(lim,n);
end

% Vehicle Frequency Loop
% for i=1:11
% mv_Mon=VehicleVariables(i,1);
% mw_Mon=VehicleVariables(i,2);
% kv_Mon=VehicleVariables(i,3); %Stiffness of vehicle spring N/m
% kw_Mon=VehicleVariables(i,6); %Stiffness of vehicle tire N/m
% cs_Mon=VehicleVariables(i,4); % Damping of vehicle spring N*s/m
% cw_Mon=VehicleVariables(i,5); % Damping of vehicle tire N*s/m
% 
% 
% K_Mon=[kv_Mon, -kv_Mon; -kv_Mon, kv_Mon+kw_Mon];
% M_Mon=[mv_Mon, 0; 0, mw_Mon];
% ei_Mon=eig(K_Mon,M_Mon); % eigenvalues
% ef_Mon=sort(real(sqrt(ei_Mon))); % sorted natural angular frequencies [rad/s] 
% Vehicle(:,i)=ef_Mon/(2*pi); % 1st and second natural frequencies of sprung mass
% 
% end