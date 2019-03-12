function [RoadMatrix]=CompleteSurfaceRoughness(RoughnessClassification,Length)
Velocity=1;
dT=.001; % Time Step
Time=0:dT:Length/Velocity; % Total time to cross bridge
DataPoints=length(Time); % Number of Data points

if RoughnessClassification==3
    gd= (.001+rand(1,1)*(sqrt(32)-.001))*10^-6;
elseif RoughnessClassification==4
    gd=128*10^-6;
elseif RoughnessClassification==5
    gd=512*10^-6;
else 
    gd=32e-6;  
end

SamplingInterval=dT;
DeltaN=1/Length;
SpatialFrequencyBand=DeltaN:DeltaN:(1/(SamplingInterval)); % Spatial Frequency Band
PhaseAngle=2*pi*rand(1,DataPoints); % Random Phase Angle
n0 = 0.1; % Spatial Frequency (cycles/m)
xx = 0:SamplingInterval:Length-SamplingInterval; % Abscissa Variable from 0 to L 
GD=gd*(SpatialFrequencyBand./n0).^-2;
d=sqrt(2*GD.*DeltaN);
rx = zeros(size(xx));
drx = zeros(size(xx));
ddrx = zeros(size(xx));
for i=1:DataPoints-1
rx(i)=sum(d.*cos(2*pi*SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints-1))));
% drx(i)=sum(-SpatialFrequencyBand.*d.*sin(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints-1))));
% ddrx(i)=sum((-(SpatialFrequencyBand).^2).*d.*cos(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints-1))));
end
% 
drx=diff(rx);
drx(DataPoints-1)=drx(1);
ddrx=diff(rx,2);
ddrx(DataPoints-2:DataPoints-1)=ddrx(1:2);
RoadMatrix=[round(xx,3); rx; drx; ddrx];

% figure(1)
% set(gcf,'color','white')
% plot(Time,drx,'color','k','Linewidth',2); hold on
% xlabel('Distance Along Bridge (m)')
% xlim([0 Length])
% ylabel('Profile Height (m)')
% title ('Class A Surface Profile')
% legend ('Simulation 1', 'Simulation 2')
% plotformat
end





