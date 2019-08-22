function [Mp]=Mp(x,y,sigma,dvect)
    t=x; h=10^(-8);
    if (x==y)
        t=x+h;
    end
    G1=inv(G(t,sigma,dvect));
    G2=G(y,sigma,dvect);
    Mp=(eye(3)-G1*G2)/(t-y);
end