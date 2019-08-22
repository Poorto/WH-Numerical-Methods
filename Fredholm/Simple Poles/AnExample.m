%{
Example solution to Wiener-Hopf equation using Fredholm Integral Equation
of the second kind
Form: G*V+ = V- - F0
Where F0 has a simple pole at x_0
%}

clear
clc

%Define Length and Spacing for Numerical Integration
h=.2;   %Spacing
A=10;   %Path Endpoints (from -A to A)

phi=3*pi()/12;
sigma=1-i*10^(-8);      %Note - Singularity must be in upper or lower plane
x0=-sigma*cos(phi);   %Location of Singularity of F0 (x_0)

%Define Function Handle of Matrix G
Gfun=@(x)[1 x*i x^2-1;2*x 1 i*x^2+1;1 x+4*i 0];
%Define Function Handle of Vector F0
Ffun=@(x)(1/(x-x0))*[2*i*sin(phi);0;0];

%Solve Fredholm IE Using Quadrature Integration
tic
S=SolFred(Gfun,Ffun,x0,A,h);
toc

resol=50;
alpha = linspace(-A,A,round(resol*A));
%Use Cauchy Expansion to Find Solution for Arbitrary Alpha
tic
Vsols = VpFred(alpha,Gfun,Ffun,A,h,x0,S);
toc

%Plot Real and Imaginary Part for Each Element of V+
for n=1:length(Vsols(:,1))
    figure(n)
    clf
    plot((alpha),real(Vsols(n,:)),'B',(alpha),imag(Vsols(n,:)),'R')
    legend('Real(V_n)','Imag(V_n)')
    title('V_+ for n=')
end