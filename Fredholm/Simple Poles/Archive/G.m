function [G]=G(alpha,sigma,dvect)
    Z0=376.7;
    g=gam(alpha,sigma);
    a1=a1n(alpha,sigma,dvect);
    a2=a2n(alpha,sigma,dvect);
    
    Yc=g/(Z0*sigma);
    Y11=i*Yc*tan(g*dvect(1)/2);
    Y31=-i*Yc/sin(g*dvect(1));
    Y12=i*Yc*tan(g*dvect(2)/2);
    Y32=-i*Yc/sin(g*dvect(2));
    
    G=[(Yc+Y11+Y31) -Y31 0;-Y31 (Y11+Y12+Y31+Y32) -Y32;0 -Y32 (Yc+Y12+Y32)];
end