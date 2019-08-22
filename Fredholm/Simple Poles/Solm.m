function [Solm]=Solm(Gfun,Ffun,sing,A,h)
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
        Solm - Vector of Solutions at A/h Points
            Length = 2*MO*(A/h+1) 
%}
    t=@(y)exp(i*pi()/4)*y;
    N=round(A/h);       %Index Length of Path
    Ftemp = Ffun(0);    %Check Solvability of F0
    MO = length(Ftemp); %Detect Vector Length
    
    %Initialize Matrices for Quadrature Integral Approximation
    D=zeros(2*MO*N+3); 
    K=zeros(2*MO*N+3);
    S=zeros(2*MO*N+3,1);
    
    for m=-N:N;
        r=MO*(m+N); %Identify Row to Place Matrix Solutions
        
        x=t(m*h);   %Warp Intgration Path
        Gt=eye(MO); %Temporary Identity Matrix
              
        
        St=Ffun(x);
        for j=1:MO
            S(r+j,1)=St(j,1);   %Use Temp Solution to Populate S
        end
        
        for j=1:MO
            for k=1:MO
                D(r+j,r+k)=Gt(j,k); %Use Temp Solution to Populate D
            end
        end
    
        
        for n=-N:N
            c=MO*(n+N);             %Identify Column Location
            y=t(n*h);               %Warped Integration Path
            Tp=exp(i*pi()/4);       %Derivative of Warped Path
            Mt=Mm(Gfun,x,y);        %Solve Matrix Kernal M(x,y)
            Bt=h*Tp*Mt/(2*pi()*i);  %Numerical Step Spacing*M(x,y)
            
            for j=1:MO
                for k=1:MO
                    K(r+j,c+k)=Bt(j,k); %Use Temp Solution to Populate K
                end
            end
            
        end
    end
    Stemp=inv(D+K)*S; %Solve Numerical Quadrature Approximation

    Solm=zeros(MO,2*N+1);   %Initialize Vector Solution
    for n=1:2*N+1
        for m = 1:MO
            Solm(m,n)=S(MO*n-MO+m); %Reorder Vector to Rows of Vectors
        end
    end

end