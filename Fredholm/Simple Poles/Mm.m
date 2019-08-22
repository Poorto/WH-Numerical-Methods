function [Mm]=Mm(Gfun,x,y)
%{
    Matrix Kernal for Function with Singularity in Lower Half Plane
    Form - G(a)*V+ = V- + F0
    Args:
        Gfun - G matrix (square matrix of size order MO)
        x - Solution for Fredholm IE Output of Solp
        x - Solution for Fredholm IE Output of Solp
    Returns:
        Mm - Matrix Kernal M- = (G(x)-G(y))/(x-y)
%}

    t=x; h=10^(-8);
    if (x==y)
        t=x+h;
    end
    G1=Gfun(t);
    G2=Gfun(y);
    Mm=(G1-G2)/(t-y);
end