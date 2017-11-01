function exitCode = Bridge_Pre_Allocation(path)

% Pre-Processing (Change values manually)
Surface=1; % 1 if surface is considered, 0 otherwise
Wind=0; % 1 if wind effects are considered, 0 otherwise
Temp=1; % 1 if temp effects are considered, 0 otherwise
RainEffects=1; % 1 if rain effects are considered, 0 otherwise
Damage=1; % 1 if damage effects are considered, 0 otherwise
Bridge=1; % Indicates which bridge is being tested
lim=361; % Number of days monitoring subject bridge
NumberElements=10; % Number of elements bridge is divided into

% Load variable arrays
BridgeVariables=load('BridgeVariables.dat');
VehicleVariables=load('VehicleData.dat');
load('USW00003870.hourly.mat'); 

%% Pre Analysis Calculations (Sets up environmental and vehicle parameters)
Td=0:.016667:24; %Time of day record was taken (Assumes vehicles cross bridge every min of a day)
n=length(Td); % Number of cases looking at a single car crossing bridge

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
[ele, nodes]=element(NumberElements,l,L);

% Surface Roughness
if Surface==1
  RoadMatrix=CompleteSurfaceRoughness(Rclass,L);
end

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
row=randi([1 3],lim,n); % Randomly selects which vehicle is crossing bridge
V=randi([10 25],lim,n); % Speed of vehicle m/s

% Storage matrices and cell arrays (used to store information for machine
% learning)
Time=cell(lim,n);
VehicleMass=zeros(lim,n);
WheelMass=zeros(lim,n);
SuspensionStiffness=zeros(lim,n);
WheelStiffness=zeros(lim,n);
SuspensionDamping=zeros(lim,n);
WheelDamping=zeros(lim,n);
AccelerationVehicle=cell(lim,n);
AllFrequencyData=cell(lim,n);
fv=zeros(lim,n);


for i=1:lim
  
if RainEffects==1;
if i>=3
if Rain(i)==1 && mu(i-1)>=1.01*BridgeVariables(Bridge,2);
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

save(path)
exitCode = 0;
end
