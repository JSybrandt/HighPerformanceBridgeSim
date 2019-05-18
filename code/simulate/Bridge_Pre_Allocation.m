function exitCode = Bridge_Pre_Allocation(output_path, ...
                                          Multiple_Vehicles, ...
                                          Damage_Case, ...
                                          Environmental_Effects, ...
                                          num_days, ...
                                          weather_data_path, ...
                                          bridge_data_path, ...
                                          vehicle_data_path)

% Parse input
Environmental_Effects = str2num(Environmental_Effects);
Multiple_Vehicles = str2num(Multiple_Vehicles);
Damage_Case = str2num(Damage_Case);
num_days = str2num(num_days);

assert((Environmental_Effects == 0) || ...
       (Environmental_Effects == 1));
Temp = Environmental_Effects;
Wind=0;

% 1 if multiple vehicles are considered, 0 if just 1 vehicle is being considered
assert((Multiple_Vehicles == 0) || ...
       (Multiple_Vehicles == 1));
% Determines which damage case being analyzed
assert((Damage_Case == 1) || ...
       (Damage_Case == 2));

assert(num_days > 0);
lim=num_days;


% Pre-Processing (Change values manually)
Surface=1; % 1 if surface is considered, 0 otherwise
Damage=1; % 1 if damage effects are considered, 0 otherwise
Bridge=1; % Indicates which bridge is being tested
NumberElements=11; % Number of elements bridge is divided into

% Load variable arrays
BridgeVariables=load(bridge_data_path);
VehicleVariables=load(vehicle_data_path);
load(weather_data_path);

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
% Continuous
DamageLocation = 5;
%DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
%DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
DayDamage1=120;
damaged_days = lim - DayDamage1 - 120;
init_damage=0.05;
total_damage=0.25;
damage_per_day = (total_damage-init_damage)/damaged_days;
ED1=[init_damage*E, ones(1, damaged_days)*damage_per_day*E];
%ED1=[((.02+rand(1)*(.02-.01))*E), .001*rand(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1
elseif Damage==1 && Damage_Case==2
DamageLocation=5;
%DamageLocation=randsample(ele(:,1),1); % Where damage locations begin
DayDamage1=120;
%DayDamage1=round(lim*.25)+round(rand(1)*(lim*.3-lim*.25)); % The day damage is iniciated on bridge
DayDamage2=240;
%DayDamage2=round(lim*.35)+round(rand(1)*(lim*.45-lim*.35)); % The day second damage is iniciated on bridge
DayDamage3=360;
%DayDamage3=round(lim*.55)+round(rand(1)*(lim*.65-lim*.55)); % The day damage is iniciated on bridge
DayDamage4=480;
%DayDamage4=round(lim*.7)+round(rand(1)*(lim*.75-lim*.7)); % The day damage is iniciated on bridge
DayDamage5=600;
%DayDamage5=round(lim*.8)+round(rand(1)*(lim*.9-lim*.8)); % The day damage is iniciated on bridge
ED1=0.05*E;
ED2=0.05*E;
ED3=0.05*E;
ED4=0.05*E;
ED5=0.05*E;
%ED1=((.03+rand(1)*(.05-.03))*E); % Damaged Modulus 1
%ED2=((.03+rand(1)*(.05-.03))*E); % Damaged Modulus 2
%ED3=((.03+rand(1)*(.05-.03))*E); % Damaged Modulus 3
%ED4=((.03+rand(1)*(.05-.03))*E); % Damaged Modulus 4
%ED5=((.03+rand(1)*(.05-.03))*E); % Damaged Modulus 5
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

    V(gg,hh,1)=randi([2 25]);
    V(gg,hh,2)=0;
elseif number_vehicles(gg,hh)==2
    row(gg,hh,1)=randi([1 8]);
    row(gg,hh,2)=randi([1 8]);

    V(gg,hh,1)=randi([2 25]);
    V(gg,hh,2)=randi([2 25]);
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
V=randi([2 25],lim,n); % Speed of vehicle m/s

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

save(output_path)
exitCode = 0;
end % End function
