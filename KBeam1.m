function [kb]=KBeam1(E,I,L);
k1=12*E*I/L^3;
k2=6*E*I/L^2;
k3=4*E*I/L;
k4=2*E*I/L;
kb=[k1,-k1;
    -k1,k1];
end
