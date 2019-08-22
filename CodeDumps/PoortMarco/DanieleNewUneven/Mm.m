function [Mm]=Mm(x,y,sigma,dvect)
    t=x; h=10^(-8);
    if (x==y)
        t=x+h;
    end
    G1=G(t,sigma,dvect);
    G2=G(y,sigma,dvect);
    Mm=(G1-G2)/(t-y);
end