function [Vpm]=Vpm(alpha,Gfun,Ffun,A,h,sing,VT)
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
        Vpm - Vector Solution at alpha
            Length = MO
%}
    t=@(y)exp(i*pi()/4)*y;
    N=length(VT(1,:));  %Identify Length of Solved Integral Equation
    x=alpha;            %Dummy Variable x=alpha for Shorthand
    tpr=exp(i*pi()/4);  %Warped Integration Line Derivative (Constant)
    
    RT=Ffun(x);         %Temporary F0(x);
    GT=Gfun(x);         %Temporary G(x)
    
    MT=zeros(3,3);
    
    for n=1:N
        y=t((n-(N+1)/2)*h);                 %Point Along Integral Path
        MT=-h*Mm(Gfun,x,y)*tpr/(2*pi()*i);  %Matrix Kernal M(x,y_n)
        T=VT(:,n);                          %Extract Vector from VT
        RT=RT+MT*T;                         %Increment Numerical Integral
    end
    
    Vpm=inv(GT)*RT; %Vector Solution for V+(alpha)
end