function [RoadMatrix]=CompleteSurfaceRoughness(RoughnessClassification,Length)
Velocity=1;
dT=.001; % Time Step
Time=0:dT:(Length+dT)/Velocity; % Total time to cross bridge
DataPoints=length(Time); % Number of Data points

if RoughnessClassification==3
    gd=(.001+rand(1,1)*(32-.001))*10^-6;
elseif RoughnessClassification==4
    gd=(32+rand(1,1)*(128-32))*10^-6;
elseif RoughnessClassification==5
    gd=(128+rand(1,1)*(512-128))*10^-6;
else 
    gd=32e-6;  
end

SamplingInterval=(Length+dT)/DataPoints;
DeltaN=(1/SamplingInterval-1)/(DataPoints-1);
SpatialFrequencyBand=1:DeltaN:(1/(SamplingInterval)); % Spatial Frequency Band
PhaseAngle=2*pi*rand(1,DataPoints); % Random Phase Angle
n0 = 0.1; % Spatial Frequency (cycles/m)
xx = 0:SamplingInterval:(Length+dT)-SamplingInterval; % Abscissa Variable from 0 to L 
GD=gd*(SpatialFrequencyBand./n0).^-2;
d=sqrt(2*GD.*DeltaN);
rx = zeros(size(xx));
for i=1:DataPoints
rx(i)=sum(d.*cos(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
drx(i)=sum(-SpatialFrequencyBand.*d.*sin(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
ddrx(i)=sum((-(SpatialFrequencyBand).^2).*d.*cos(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
end
RoadMatrix=[round(xx,3); rx/(2*pi); drx/(2*pi); ddrx/(2*pi)];

% figure(1)
% set(gcf,'color','white')
% plot(Time,rx,'color','r','Linewidth',2); hold on
% xlabel('Distance Along Bridge (m)')
% xlim([0 Length])
% ylabel('Profile Height (m)')
% title ('Class A Surface Profile')
% legend ('Simulation 1', 'Simulation 2')
% plotformat
end





