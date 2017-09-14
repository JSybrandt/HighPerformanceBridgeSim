function [kb]=KBeam(E,I,l);
k1=12*E*I/l^3;
k2=6*E*I/l^2;
k3=4*E*I/l;
k4=2*E*I/l;
kb=[k1,k2,-k1,k2;
    k2,k3,-k2,k4;
    -k1,-k2,k1,-k2;
    k2,k4,-k2,k3];
end
