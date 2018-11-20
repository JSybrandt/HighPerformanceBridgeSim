function exitCode = Bridge_Pre_Allocation(output_path, ...
                                          Multiple_Vehicles, ...
                                          Damage_Case, ...
                                          Environmental_Effects)

% Parse input
Environmental_Effects = str2num(Environmental_Effects);
Multiple_Vehicles = str2num(Multiple_Vehicles);
Damage_Case = str2num(Damage_Case);

assert((Environmental_Effects == 0) || ...
       (Environmental_Effects == 1));
% Set enviromental flags based on input
Temp = Environmental_Effects;
RainEffects = Environmental_Effects;
Surface = Environmental_Effects;
Wind=0;

% 1 if multiple vehicles are considered, 0 if just 1 vehicle is being considered
assert((Multiple_Vehicles == 0) || ...
       (Multiple_Vehicles == 1));
% Determines which damage case being analyzed
assert((Damage_Case == 1) || ...
       (Damage_Case == 2));

% Pre-Processing (Change values manually)
Damage=1; % 1 if damage effects are considered, 0 otherwise
Bridge=1; % Indicates which bridge is being tested
lim=720; % Number of days monitoring subject bridge
NumberElements=10; % Number of elements bridge is divided into

% Load variable arrays
BridgeVariables=load('data/BridgeVariables.dat');
VehicleVariables=load('data/VehicleData.dat');
load('data/USW00003870.hourly.double.mat');

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

if Damage == 1
  DamageClass=zeros(lim);
end
% Damage Variables
if Damage==1 && Damage_Case==1
DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
ED1=[((.05+rand(1)*(.1-.05))*E),.0025*rand(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1

elseif Damage==1 && Damage_Case==2
DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.3-lim*.25)); % The day damage is iniciated on bridge
DayDamage2=round(lim*.35)+round(rand(1)*(lim*.45-lim*.35)); % The day second damage is iniciated on bridge
DayDamage3=round(lim*.55)+round(rand(1)*(lim*.65-lim*.55)); % The day damage is iniciated on bridge
DayDamage4=round(lim*.7)+round(rand(1)*(lim*.75-lim*.7)); % The day damage is iniciated on bridge
DayDamage5=round(lim*.8)+round(rand(1)*(lim*.9-lim*.8)); % The day damage is iniciated on bridge
ED1=((.05+rand(1)*(.1-.05))*E); % Damaged Modulus 1
ED2=((.05+rand(1)*(.1-.05))*E); % Damaged Modulus 2
ED3=((.05+rand(1)*(.1-.05))*E); % Damaged Modulus 3
ED4=((.05+rand(1)*(.1-.05))*E); % Damaged Modulus 4
ED5=((.05+rand(1)*(.1-.05))*E); % Damaged Modulus 5
end

% Environment matrices (Rain, Temp, Wind)
if RainEffects==1
Rain=[0,0,randi([0 5],1,lim+1)];
Rain(Rain<4)=0;
Rain(Rain>0)=1;
end
if Temp==1
    Tact=zeros(lim,n); %Actual temperature degrees Celsius
    hour=0;
end
if Wind==1
    WindVelocity=zeros(lim,n);
    ForceWind=zeros(lim,n);
    FW=zeros(2*NumberElements,1);
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
Monitor_Vehicle_Frequency_Amp_Data=cell(lim,n);
Monitor_Vehicle_Frequency_Data=cell(lim,n);
Monitor_Vehicle_Road_Profile=cell(lim,n);
Monitor_Vehicle_Derivative_Road_Profile=cell(lim,n);
Monitor_Vehicle_Other_Derivative_Road_Profile=cell(lim,n);

Other_Vehicle_Time=cell(lim,n);
Other_Vehicle_Acceleration=cell(lim,n);
Other_Vehicle_Frequency_Amp_Data=cell(lim,n);
Other_Vehicle_Frequency_Data=cell(lim,n);
Other_Vehicle_Road_Profile=cell(lim,n);
Other_Vehicle_Derivative_Road_Profile=cell(lim,n);
Other_Vehicle_Other_Derivative_Road_Profile=cell(lim,n);

Start_Time_Following_Vehicle=zeros(lim,n);
Order_of_Vehicles=cell(lim,n);

for gg=1:lim
    for hh=1:n
if number_vehicles(gg,hh)==1
    row(gg,hh,1)=randi([1 10]);
    row(gg,hh,2)=0;
%     row(gg,hh,3)=0;

    V(gg,hh,1)=randi([10 25]);
    V(gg,hh,2)=0;
%     V(gg,hh,3)=0;
elseif number_vehicles(gg,hh)==2
    row(gg,hh,1)=randi([1 10]);
    row(gg,hh,2)=randi([1 10]);
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

row=randi([1 10],lim,n); % Randomly selects which vehicle is crossing bridge
V=randi([10 25],lim,n); % Speed of vehicle m/s

% Storage matrices and cell arrays (used to store information for machine
% learning)
Monitor_Vehicle_Time=cell(lim,n);
Monitor_Vehicle_Acceleration=cell(lim,n);
Monitor_Vehicle_Frequency_Amp_Data=cell(lim,n);
Monitor_Vehicle_Frequency_Data=cell(lim,n);
Monitor_Vehicle_Road_Profile=cell(lim,n);
Monitor_Vehicle_Derivative_Road_Profile=cell(lim,n);
Monitor_Vehicle_Other_Derivative_Road_Profile=cell(lim,n);
end

if RainEffects==1
    for i=1:lim

        if i>=3
if Rain(i)==1 && mu(i-1)>=1.01*BridgeVariables(Bridge,2)
       mu(i)=1.01*BridgeVariables(Bridge,2);
elseif Rain(i)==1
       mu(i)=mu(i-1)+.001*rand(1,1)*BridgeVariables(Bridge,2);
elseif Rain(i)==0 && Rain(i-1)==1
   mu(i)=mu(i-1)-.001*rand(1,1)*BridgeVariables(Bridge,2);
elseif Rain(i)==0 && mu(i-1)>BridgeVariables(Bridge,2)
           mu(i)=mu(i-1)-.001*rand(1,1)*BridgeVariables(Bridge,2);
elseif Rain(i)==0 && mu(i-1)<=.99*BridgeVariables(Bridge,2)
           mu(i)=.99*BridgeVariables(Bridge,2);
end
        end
    end
end

save(output_path)
exitCode = 0;
end
