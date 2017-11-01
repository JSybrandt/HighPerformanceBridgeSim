function [RoadMatrix]=CompleteSurfaceRoughness(RoughnessClassification,Length)
Velocity=1;
dT=.001; % Time Step
Time=0:dT:Length/Velocity; % Total time to cross bridge
DataPoints=length(Time); % Number of Data points

if RoughnessClassification==3
    gd=(.001+rand(1,1)*(sqrt(32)-.001))*10^-6;
else if RoughnessClassification==4
    gd=(sqrt(32)+rand(1,1)*(sqrt(128)-sqrt(32)))*10^-6;
else if RoughnessClassification==5
    gd=(sqrt(128)+rand(1,1)*(sqrt(512)-sqrt(128)))*10^-6;
else 
    gd=32e-6;  
    end
    end
end

SamplingInterval=Length/DataPoints;
DeltaN=(1/(SamplingInterval)-1)/(DataPoints-1);
SpatialFrequencyBand=1:DeltaN:(1/(SamplingInterval)); % Spatial Frequency Band
PhaseAngle=2*pi*rand(1,DataPoints); % Random Phase Angle
n0 = 0.1; % Spatial Frequency (cycles/m)
xx = 0:SamplingInterval:Length-SamplingInterval; % Abscissa Variable from 0 to L (+3*Length/DataPoints)
GD=gd*(SpatialFrequencyBand./n0).^-2;
d=sqrt(2*GD.*DeltaN);
rx = zeros(size(xx));
for i=1:DataPoints
rx(i)=sum(d.*cos(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
drx(i)=sum(-SpatialFrequencyBand.*d.*sin(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
ddrx(i)=sum((-(SpatialFrequencyBand).^2).*d.*cos(SpatialFrequencyBand.*xx(i)+PhaseAngle(1:(DataPoints))));
end
RoadMatrix=[round(xx,3); rx/(2*pi); drx/(2*pi); ddrx/(2*pi)];
end





