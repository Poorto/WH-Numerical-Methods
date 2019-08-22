clear
clc

Testing=0;

tic
h=.1;
A=160*h;
N=round(A/h);

sigma=1-i*10^(-8);
dvect=[2.0,1.0];
phi=9*pi()/12;

S=Solm(A,h,sigma,dvect,phi);

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
N=2000;
M=5*real(sigma);
Alpha=linspace(-M,M,2*N+1)+b*i;

FP=zeros(3,2*N+1);
c=zeros(1,2*N+1);

for n=1:2*N+1
    c(1,n)=sqrt(sigma-Alpha(n));
    FP(:,n)=Vpm(Alpha(n),A,h,sigma,dvect,phi,V1p,V2p,V3p);
    COE(:,n)=ABCDEF(Alpha(n),sigma,dvect,phi,FP(1,n),FP(2,n),FP(3,n));
end

toc

if (Testing==0)
figure(1)
clf
plot(real(Alpha),abs(c(1,:).*(FP(1,:)+FP(3,:))),'K')
axis([-M M 0 120])
grid on

figure(2)
clf
plot(real(Alpha),abs(2*c.*FP(2,:)),'K')
axis([-M M 0 120])
grid on

end

if (Testing==1)
figure(1)
clf
plot(real(Alpha),abs(COE(1,:)),'K')
axis([-M M 0 10])
grid on

figure(2)
clf
plot(real(Alpha),abs(COE(2,:)),'K')
axis([-M M 0 120])
grid on

figure(3)
clf
plot(real(Alpha),abs(COE(3,:)),'K')
axis([-M M 0 120])
grid on

figure(4)
clf
plot(real(Alpha),abs(COE(4,:)),'K')
axis([-M M 0 120])
grid on

figure(5)
clf
plot(real(Alpha),abs(COE(5,:)),'K')
axis([-M M 0 120])
grid on

figure(6)
clf
plot(real(Alpha),abs(COE(6,:)),'K')
axis([-M M 0 120])
grid on
end

if (Testing==2)
tic
M=050; 
x=linspace(-5*sum(dvect),5*sum(dvect),M);
y=linspace(0.0*sum(dvect),10*sum(dvect),M);
[X,Y]=meshgrid(x,y);
U=-0*exp(i*sigma*(X*cos(phi)-Y*sin(phi)));
W=exp(i*sigma*(X*cos(phi)+Y*sin(phi)));

for m=1:M
    for n=1:M
        U(m,n)=U(m,n)+Us1(X(m,n),Y(m,n),Alpha,sigma,dvect,COE(1,:));
    end
end

figure(7)
clf
colormap(hot);
surf(X,Y,real(U),'LineStyle','none');
colorbar;
view(2);
figure(8)
clf
colormap(hot);
surf(X,Y,imag(U),'LineStyle','none');
colorbar;
view(2);
figure(9)
clf
colormap(hot);
surf(X,Y,abs(U),'LineStyle','none');
colorbar;
view(2);
toc
end

if (Testing==3)
tic
M=050; 
x=linspace(-5*sum(dvect),5*sum(dvect),M);
y=linspace(-12*sum(dvect),-2*sum(dvect),M);
[X,Y]=meshgrid(x,y);
U=0*exp(i*sigma*(X*cos(phi)+Y*sin(phi)));

for m=1:M
    for n=1:M
        U(m,n)=U(m,n)+Us3(X(m,n),Y(m,n),Alpha,sigma,dvect,COE(6,:));
    end
end

figure(7)
clf
colormap(hot);
surf(X,Y,real(U),'LineStyle','none');
colorbar;
view(2);
figure(8)
clf
colormap(hot);
surf(X,Y,imag(U),'LineStyle','none');
colorbar;
view(2);
figure(9)
clf
colormap(hot);
surf(X,Y,abs(U),'LineStyle','none');
colorbar;
view(2);
toc
end



