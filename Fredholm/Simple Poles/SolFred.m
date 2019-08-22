function [SolFred]=SolFred(Gfun,Ffun,sing,A,h)
%{
    Solve Numerical Fredholm IE of Second Kind Using Quadrature Integration
    Form - G(a)*U+ = U- + F0
    Args:
        Gfun - G matrix (square matrix of size order MO)
        Ffun - F0 vector (length MO)
        sing - Location of 1st order Singularity in F0
        A - Length of Integration Line to Use
        h - Spacing Between Points for Quatrature Integration
    Returns:
        SolFred - Vector of Solutions at A/h Points
            Length = 2*MO*(A/h+1) 
%}

    if imag(sing)>0
        SolFred = Solp(Gfun,Ffun,sing,A,h);
    else
        SolFred = Solm(Gfun,Ffun,sing,A,h);    
    end
end
    