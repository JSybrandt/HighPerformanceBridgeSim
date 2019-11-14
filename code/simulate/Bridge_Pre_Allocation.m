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
Surface=0; % 1 if surface is considered, 0 otherwise
Damage=1; % 1 if damage effects are considered, 0 otherwise
Bridge=1; % Indicates which bridge is being tested
NumberElements=11; % Number of elements bridge is divided into

% Load variable arrays

if Damage==1 && Damage_Case==1
NumberElements=21; % Number of elements bridge is divided into
else
NumberElements=40;
end

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
mu=BridgeVariables(Bridge,2); % mass per unit length kg/m
bbeta=BridgeVariables(Bridge,5); % Damping of bridge
Rclass=BridgeVariables(Bridge,7);% Values For ISO Road Roughness Classification, from 3 to 9
[ele, ele2, nodes, nodes2]=element(NumberElements,l,L,Multiple_Vehicles);

% Surface Roughness
if Surface==1
[RoadMatrix]=CompleteSurfaceRoughness(Rclass,L);
end

% Damage Variables
if Damage==1 && Damage_Case==1
DamageLocation=round(NumberElements/2,0); % Where damage locations begin
DayDamage1=round(lim*.25,0); % The day damage is iniciated on bridge
DayDamage2=round(lim*.33,0); % The day second damage is iniciated on bridge
DayDamage3=round(lim*.5,0); % The day damage is iniciated on bridge
DayDamage4=round(lim*.66,0); % The day damage is iniciated on bridge
DayDamage5=round(lim*.75,0); % The day damage is iniciated on bridge
% Experiment 1
%ED1=(.169*E); % Damaged Modulus 1
%ED2=(.292*E); % Damaged Modulus 2
%ED3=(.386*E); % Damaged Modulus 3
%ED4=(.46*E); % Damaged Modulus 4
%ED5=(.496*E); % Damaged Modulus 5
% Experiment 2
ED1=(.097*E); % Damaged Modulus 1
ED2=(.17*E); % Damaged Modulus 2   
ED3=(.245*E); % Damaged Modulus 3
ED4=(.293*E); % Damaged Modulus 4 
ED5=(.354*E); % Damaged Modulus 5  


elseif Damage==1 && Damage_Case==2
DayDamage1=round(lim*.25,0); % The day damage is iniciated on bridge
DayDamage2=round(lim*.33,0); % The day second damage is iniciated on bridge
DayDamage3=round(lim*.5,0); % The day damage is iniciated on bridge
DayDamage4=round(lim*.66,0); % The day damage is iniciated on bridge
DayDamage5=round(lim*.75,0); % The day damage is iniciated on bridge
ED=(.2918*E); % Damaged Modulus 1

end

if Temp==1
    Tact=zeros(lim,n); %Actual temperature degrees Celsius
end


% Randomized vehicle selection and speed
if Multiple_Vehicles==1
number_vehicles=randi([1 2],lim,n); % Randomly determines how many vehicles will cross bridge
row=zeros(lim,n,2);
V=zeros(lim,n,2);

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
    V(gg,hh,1)=randi([2 25]);


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
