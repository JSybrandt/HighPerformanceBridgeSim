function [rx,drx,ddrx,d]=surfaceroughness(k,Length,Velocity,SampleInterval,DataPoints,SpatialFrequencyBand,RandomPhaseAngle)
if k==3
    gd=.001e-6;
else if k==4
    gd=8e-6;
else if k==5
    gd=16e-6;
else 
    gd=32e-6; 
    end
    end
end
SI=Velocity*SampleInterval;
n0 = 0.1; % Spatial Frequency (cycles/m)
xx = 0:Length/DataPoints:Length+Length/DataPoints; % Abscissa Variable from 0 to L
GD=gd*(SpatialFrequencyBand./n0).^-2;
d=sqrt(2*GD.*SI);
rx = zeros(size(xx));
for i=1:length(xx)
rx(i)=sum(d.*cos(SpatialFrequencyBand.*xx(i)+RandomPhaseAngle(1:(DataPoints+2))));
end
drx=diff(rx);
ddrx=diff(rx,2);
% dGD=diff(GD);
% ddGD=diff(GD,2);
% figure(1)
% plot(1:(NumberofDataPoints+2),rx,'color','b'); hold on
% plot(1:(NumberofDataPoints+1),drx,'color','r'); hold on
% plot(1:(NumberofDataPoints),ddrx,'color','g');
% % figure(2)
% [px,w]=pwelch(rx,[],[],[],B);
% [dpx]=pwelch(drx,[],[],[],B);
% [ddpx]=pwelch(ddrx,[],[],[],B);
% % loglog(w,pxx); hold on

% [q,C] = psd_1D(hx, B, 'x')  % B is Sampling Interval (m); for the case that I have explained above it will be 250/45000 = 5.55e-3
% lambda = (2*pi)./q; % wavelengths
% f = q/(2*pi); % spatial frequency
% PSD = 2*pi*C; % power spectrum
% loglog(f,PSD)
end





