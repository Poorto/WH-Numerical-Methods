function [Vpm]=Vpm(alpha,A,h,sigma,dvect,phi,V1,V2,V3)
    N=length(V1(1,:));
    x=alpha;
    tpr=wp();
    
    VT=[V1;V2;V3];
    RT=inv(G(-sigma*cos(phi),sigma,dvect))*F0(x,sigma,dvect,phi);
    GT=G(x,sigma,dvect);
    
    MT=zeros(3,3);
    
    for n=1:N
        y=t((n-(N+1)/2)*h);
        MT=-h*Mp(x,y,sigma,dvect)*tpr/(2*pi()*i);
        T=VT(:,n);
        RT=RT+MT*T;
    end
    
    Vpm=RT;
end