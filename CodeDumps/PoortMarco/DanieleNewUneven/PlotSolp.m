clear
clc

Testing=1;

tic
h=.500;
A=30;
N=round(A/h);

theta=4*pi()/12;
phi=3*pi()/12;
dvect=[2.25,0.75];

k0=2*pi()-i*10^(-8);
Z0=376.7; Y0=1/Z0;
sigma=k0*sin(theta);

S=Solp(A,h,sigma,dvect,phi);

V1p=zeros(1,2*N+1);
V2p=zeros(1,2*N+1);
V3p=zeros(1,2*N+1);

b=imag((sigma+sigma*cos(phi))/2);
Alpha1=linspace(-A,A,2*N+1);

for n=1:2*N+1
    V1p(1,n)=S(3*n-2);
    V2p(1,n)=S(3*n-1);
    V3p(1,n)=S(3*n-0);
end
toc

tic
N=2.0*round(2*sigma*max(dvect)+0.5);
for n=1:N
    Gamm1(n)=n*pi()/dvect(1);
    Gamm2(n)=n*pi()/dvect(2);
end
Alpha1=-sqrt(sigma^2-Gamm1.^2);
Alpha2=-sqrt(sigma^2-Gamm2.^2);

FP1=zeros(3,N);
FP2=zeros(3,N);
for n=1:N
    FP1(:,n)=Vpp(Alpha1(n),A,h,sigma,dvect,phi,V1p,V2p,V3p);
    FP2(:,n)=Vpp(Alpha2(n),A,h,sigma,dvect,phi,V1p,V2p,V3p);
    c1a(n)=i^(2*n-1)*n*pi()*FP1(1,n)/(Alpha1(n)*dvect(1)^2);
    c1b(n)=i^(2*n-1)*n*pi()*FP1(2,n)/(Alpha1(n)*dvect(1)^2);
    c2a(n)=i^(2*n-1)*n*pi()*FP2(2,n)/(Alpha2(n)*dvect(2)^2);
    c2b(n)=i^(2*n-1)*n*pi()*FP2(3,n)/(Alpha2(n)*dvect(2)^2);
    
end
toc

tic
M=500;
x=linspace(-sum(dvect),0,M);
z=linspace(0,sum(dvect),M);
y1=linspace(-dvect(1),0,M);
y2=linspace(-sum(dvect),-dvect(1),M);

[X1,Y1]=meshgrid(x,y1);
[X2,Y2]=meshgrid(x,y2);
[Xa,Za]=meshgrid(x,z);
[Xb,Zb]=meshgrid(x,z);

U2=0*X1;
U3=0*X2;

JX2=0*X1;
JX3=0*X2;

JY2=0*X1;
JY3=0*X2;

JZ2=0*Za;
JZ3=0*Zb;

for n=1:N
    U2a=(c1a(n)*sin(Gamm1(n)*(y1+dvect(1))).')*exp(-i*Alpha1(n)*x);
    U2b=(c1b(n)*sin(Gamm1(n)*(y1)).')*exp(-i*Alpha1(n)*x);
    U3a=(c2a(n)*sin(Gamm2(n)*(y2+sum(dvect))).')*exp(-i*Alpha2(n)*x);
    U3b=(c2b(n)*sin(Gamm2(n)*(y2+dvect(1))).')*exp(-i*Alpha2(n)*x);
    
    JX2a=-i*Alpha1(n)*(c1a(n)*sin(Gamm1(n)*(y1+dvect(1))).')*exp(-i*Alpha1(n)*x);
    JX2b=-i*Alpha1(n)*(c1b(n)*sin(Gamm1(n)*(y1)).')*exp(-i*Alpha1(n)*x);
    JX3a=-i*Alpha2(n)*(c2a(n)*sin(Gamm2(n)*(y2+sum(dvect))).')*exp(-i*Alpha2(n)*x);
    JX3b=-i*Alpha2(n)*(c2b(n)*sin(Gamm2(n)*(y2+dvect(1))).')*exp(-i*Alpha2(n)*x);
    
    JY2a=Gamm1(n)*(c1a(n)*cos(Gamm1(n)*(y1+dvect(1))).')*exp(-i*Alpha1(n)*x);
    JY2b=Gamm1(n)*(c1b(n)*cos(Gamm1(n)*(y1)).')*exp(-i*Alpha1(n)*x);
    JY3a=Gamm2(n)*(c2a(n)*cos(Gamm2(n)*(y2+sum(dvect))).')*exp(-i*Alpha2(n)*x);
    JY3b=Gamm2(n)*(c2b(n)*cos(Gamm2(n)*(y2+dvect(1))).')*exp(-i*Alpha2(n)*x);
    
    JZ2a=Gamm1(n)*(c1a(n)*cos(Gamm1(n)*(0)))*(cos(k0*z*cos(theta)).')*exp(-i*Alpha1(n)*x);
    JZ2b=Gamm1(n)*(c1b(n)*cos(Gamm1(n)*(-dvect(1))))*(cos(k0*z*cos(theta)).')*exp(-i*Alpha1(n)*x);
    JZ3a=Gamm2(n)*(c2a(n)*cos(Gamm2(n)*(dvect(2))))*(cos(k0*z*cos(theta)).')*exp(-i*Alpha2(n)*x);
    JZ3b=Gamm2(n)*(c2b(n)*cos(Gamm2(n)*(0)))*(cos(k0*z*cos(theta)).')*exp(-i*Alpha2(n)*x);    
    
    U2=U2+U2a-U2b;
    U3=U3+U3a-U3b;
    
    JX2=JX2+JX2a-JX2b;
    JX3=JX3+JX3a-JX3b;
    
    JY2=JY2+JY2a-JY2b;
    JY3=JY3+JY3a-JY3b;
    
    JZ2=JZ2-JZ2a+JZ2b;
    JZ3=JZ3-JZ3a+JZ3b;
end

JX2=2*i*Y0/sigma*JX2;
JX3=2*i*Y0/sigma*JX3;

JY2=2*i*Y0/sigma*JY2;
JY3=2*i*Y0/sigma*JY3;

JZ2=2*i*Y0/sigma*JZ2;
JZ3=-2*i*Y0/sigma*JZ3;

figure(1)
clf
surf(X1,Y1,abs(JX2),'LineStyle','none')
hold on
surf(X2,Y2,abs(JX3),'LineStyle','none')

title('|J_x| on Truncating Plane z=0');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(2)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');

figure(2)
clf
surf(X1,Y1,abs(JY2),'LineStyle','none')
hold on
surf(X2,Y2,abs(JY3),'LineStyle','none')

title('|J_y| on Truncating Plane z=0');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(2)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');

figure(3)
clf
surf(Xa,Za,abs(JZ2),'LineStyle','none')

title('|J_z| on Top of Plate z=-d');
xlabel('x/\lambda_0');
ylabel('z/\lambda_0');

colormap(hot);
view(2)
axis([-sum(dvect) 0 0 sum(dvect)]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');

figure(4)
clf
surf(Xb,Zb,abs(JZ3),'LineStyle','none')

title('|J_z| on Bottom of Pate z=-d');
xlabel('x/\lambda_0');
ylabel('z/\lambda_0');

colormap(hot);
view(2)
axis([-sum(dvect) 0 0 sum(dvect)]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');

% Continuity of Magnitude
figure(5)
clf
az=45;
el=30;

subplot(1,2,1)
surf(Xb,0*Xb-dvect(1),Zb,abs(JZ2),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,abs(JY2),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,abs(JY3),'LineStyle','none')

title('|J_{y/z}| on Truncating Plane z=0 and Plate y=d_1^+');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90+az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');

figure(5)
subplot(1,2,2)
surf(Xb,0*Xb-dvect(1),Zb,abs(JZ3),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,abs(JY3),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,abs(JY2),'LineStyle','none')

title('|J_{y/z}| on Truncating Plane z=0 and Plate y=d_1^-');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90-az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([0 .02]);
ylabel(h, 'A/m');


% Continuity of Real
figure(6)
clf
az=45;
el=30;

subplot(1,2,1)
surf(Xb,0*Xb-dvect(1),Zb,real(JZ2),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,real(JY2),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,real(JY3),'LineStyle','none')

title('Re(J_{y/z}) on Truncating Plane z=0 and Plate y=d_1^+');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90+az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([-.02 .02]); 
ylabel(h, 'A/m');

figure(6)
subplot(1,2,2)
surf(Xb,0*Xb-dvect(1),Zb,real(JZ3),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,real(JY3),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,real(JY2),'LineStyle','none')

title('Re(J_{y/z}) on Truncating Plane z=0 and Plate y=d_1^-');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90-az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([-.02 .02]); 
ylabel(h, 'A/m');


% Continuity of Imag
figure(7)
clf
az=45;
el=30;

subplot(1,2,1)
surf(Xb,0*Xb-dvect(1),Zb,imag(JZ2),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,imag(JY2),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,imag(JY3),'LineStyle','none')

title('Im(J_{y/z}) on Truncating Plane z=0 and Plate y=d_1^+');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90+az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([-.02 .02]);
ylabel(h, 'A/m');

figure(7)
subplot(1,2,2)
surf(Xb,0*Xb-dvect(1),Zb,imag(JZ3),'LineStyle','none')
hold on
surf(X2,Y2,0*Y2,imag(JY3),'LineStyle','none')
hold on
surf(X1,Y1,0*Y1,imag(JY2),'LineStyle','none')

title('Im(J_{y/z}) on Truncating Plane z=0 and Plate y=d_1^-');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');

colormap(hot);
view(90-az,el)
axis([-sum(dvect) 0 -sum(dvect) 0]);
h=colorbar;
caxis([-.02 .02]);
ylabel(h, 'A/m');


toc


