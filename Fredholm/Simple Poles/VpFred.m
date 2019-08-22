function [VpFred]=VpFred(alpha,Gfun,Ffun,A,h,sing,VT)
%{
    Cauchy Expansion for Solved Fredholm Integral Equation
    Form - G(a)*V+ = V- + F0
    Args:
        alpha - Complex value to Solve for Vp
        Gfun - G matrix (square matrix of size order MO)
        Ffun - F0 vector (length MO)
        sing - Location of 1st order Singularity in F0
        A - Length of Integration Line to Use
        h - Spacing Between Points for Quatrature Integration
        VT - Solution for Fredholm IE Output of Solp
    Returns:
        VpFred - Vector Solution at alpha
            Length = MO
%}

    Fcheck = Ffun(0);
    MO = length(Fcheck);
    N = length(alpha);
    VpFred = zeros(MO,N);
    
    if imag(sing)>0
        VpFun = @(x)Vpp(x,Gfun,Ffun,A,h,sing,VT);
    else
        VpFun = @(x)Vpm(x,Gfun,Ffun,A,h,sing,VT);
    end
    
    for n=1:N;
        VpFred(:,n)=VpFun(alpha(n));
    end
end

