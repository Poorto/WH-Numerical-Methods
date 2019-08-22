function [F0]=F0(alpha,sigma,dvect,phi)
    n0=-sigma*cos(phi);
    C=1/(alpha-n0);
    Z0=376.7;
    
    F0=C*[2*i*sin(phi)/Z0;0;0];
end