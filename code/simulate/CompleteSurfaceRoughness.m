function [RoadMatrix]=CompleteSurfaceRoughness(RoughnessClassification,Length)
Velocity=1;
dT=.001; % Time Step
Time=0:dT:Length/Velocity; % Total time to cross bridge
DataPoints=length(Time); % Number of Data points

if RoughnessClassification==3
    gd= (.001+rand(1,1)*(sqrt(16)-.001))*10^-6;
elseif RoughnessClassification==4
    gd=128*10^-6;
elseif RoughnessClassification==5
    gd=512*10^-6;
else 
    gd=32e-6;  
end

SamplingInterval=dT;
DeltaN=1/Length;
SpatialFrequencyBand=0:DeltaN:(1/(SamplingInterval)); % Spatial Frequency Band
PhaseAngle=2*pi*rand(1,DataPoints); % Random Phase Angle
n0 = 0.1; % Spatial Frequency (cycles/m)
xx = 0:SamplingInterval:Length; % Abscissa Variable from 0 to L 
GD=gd*(SpatialFrequencyBand./n0).^-2;
d=sqrt(2*GD.*DeltaN);
d(1)=0;
rx = zeros(size(xx));
drx = zeros(size(xx));
ddrx = zeros(size(xx));
for i=1:DataPoints
rx(i)=sum(d.*cos(2*pi*SpatialFrequencyBand.*xx(i)+PhaseAngle));
drx(i)=sum(-2*pi*SpatialFrequencyBand.*d.*sin(2*pi*SpatialFrequencyBand.*xx(i)+PhaseAngle));
ddrx(i)=sum((-(2*pi*SpatialFrequencyBand).^2).*d.*cos(2*pi*SpatialFrequencyBand.*xx(i)+PhaseAngle));
end

RoadMatrix=[xx; rx; drx; ddrx];

% figure(1)
% set(gcf,'color','white')
% plot(xx,rx,'color','k','Linewidth',2); hold on
% xlabel('Distance Along Bridge (m)')
% xlim([0 Length])
% ylabel('Profile Height (m)')
% title ('Class A Surface Profile')
% legend ('Simulation 1', 'Simulation 2')
% plotformat
end
