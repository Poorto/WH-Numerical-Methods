clear
clc

Testing=1;

tic
h=.1;
A=160*h;
N=round(A/h);

sigma=2*pi()-i*10^(-8);
d=4.1;
phi=9*pi()/12;

S=Solm(A,h,sigma,d,phi);

V1p=zeros(1,2*N+1);
V2p=zeros(1,2*N+1);
V3p=zeros(1,2*N+1);

b=imag((sigma+sigma*cos(phi))/2);
Alpha=linspace(-A,A,2*N+1);


for n=1:2*N+1
    V1p(1,n)=S(3*n-2);
    V2p(1,n)=S(3*n-1);
    V3p(1,n)=S(3*n-0);
end
toc

tic
N=round(2*sigma*d+0.5);
for n=1:N
    Gamm(n)=n*pi()/d;
end
Alpha=-sqrt(sigma^2-Gamm.^2);

FP=zeros(3,N);
for n=1:N
    FP(:,n)=Vpm(Alpha(n),A,h,sigma,d,phi,V1p,V2p,V3p);
    c1a(n)=i^(2*n-1)*n*pi()*FP(1,n)/(Alpha(n)*d^2);
    c1b(n)=i^(2*n-1)*n*pi()*FP(2,n)/(Alpha(n)*d^2);
    c2a(n)=i^(2*n-1)*n*pi()*FP(2,n)/(Alpha(n)*d^2);
    c2b(n)=i^(2*n-1)*n*pi()*FP(3,n)/(Alpha(n)*d^2);
    
end
toc

tic
M=200;
x=linspace(-2*d,0,M);
y1=linspace(-d,0,M);
y2=linspace(-2*d,-d,M);
[X1,Y1]=meshgrid(x,y1);
[X2,Y2]=meshgrid(x,y2);

U2=0*X1;
U3=0*X2;

for n=1:N
    U2a=(c1a(n)*sin(Gamm(n)*(y1+d)).')*exp(-i*Alpha(n)*x);
    U2b=(c1b(n)*sin(Gamm(n)*(y1)).')*exp(-i*Alpha(n)*x);
    U3a=(c2a(n)*sin(Gamm(n)*(y2+2*d)).')*exp(-i*Alpha(n)*x);
    U3b=(c2b(n)*sin(Gamm(n)*(y2+d)).')*exp(-i*Alpha(n)*x);
    
    U2=U2+U2a-U2b;
    U3=U3+U3a-U3b;
end

X=[X1;X2];Y=[Y1;Y2]; U=[U2;U3];

figure(1)
surf(X,Y,abs(U),'LineStyle','none')
colormap(hot);
view(2)
axis([-2*d 0 -2*d 0]);
toc

