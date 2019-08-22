function [Co]=ABCDEF(alpha,sigma,dvect,phi,V1,V2,V3)
    g=gam(alpha,sigma); 
    a1=a1n(alpha,sigma,dvect); a2=a2n(alpha,sigma,dvect);
    c=exp(-i*sigma*dvect*sin(phi)); d=-i/(alpha+sigma*cos(phi));
    
    V1=V1;
    V2=V2;
    V3=V3;
    
    sing=10^(-10);
    if(abs(abs(alpha)-real(sigma))>sing)
    A=V1;
    B=(V1*a1-V2)/(a1-1/a1);
    C=(V1-a1*V2)/(1-a1^2);
    D=(a2*V2-V3)/(1-a1*a2);
    E=(V2/a2-V3)/(1-a1*a2);
    F=V3;
    
    Co=[A;B;C;D;E;F];
    else
    Co=[0;0;0;0;0;0];
    end
    
end