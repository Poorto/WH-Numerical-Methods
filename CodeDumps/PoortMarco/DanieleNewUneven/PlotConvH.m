clear
clc



theta=4*pi()/12;
phi=3*pi()/12;
dvect=[2.25,0.75];

k0=2*pi()-i*10^(-8);
Z0=376.7; Y0=1/Z0;
sigma=k0*sin(theta);

hrange=[9.0 5.0 2.0 1.0 0.5 0.1];
Nmax=2.0*round(2*sigma*max(dvect)+0.5);
v1h = zeros(length(hrange),Nmax);
v2h = zeros(length(hrange),Nmax);


for elem = 1:length(hrange)
    htry=hrange(elem);

    A=30;
    N=round(A/htry);

    S=Solp(A,htry,sigma,dvect,phi);

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
        FP1(:,n)=Vpp(Alpha1(n),A,htry,sigma,dvect,phi,V1p,V2p,V3p);
        FP2(:,n)=Vpp(Alpha2(n),A,htry,sigma,dvect,phi,V1p,V2p,V3p);
        v1h(elem,n)=FP1(1,n);
        v2h(elem,n)=FP1(2,n);
        c1a(n)=i^(2*n-1)*n*pi()*FP1(1,n)/(Alpha1(n)*dvect(1)^2);
        c1b(n)=i^(2*n-1)*n*pi()*FP1(2,n)/(Alpha1(n)*dvect(1)^2);
        c2a(n)=i^(2*n-1)*n*pi()*FP2(2,n)/(Alpha2(n)*dvect(2)^2);
        c2b(n)=i^(2*n-1)*n*pi()*FP2(3,n)/(Alpha2(n)*dvect(2)^2);
    end
    toc
end

M = length(hrange)
col=jet(M);
nindex = 1:Nmax;

figure(2)
clf
subplot(1,2,1)
for elem=1:M
    semilogy(nindex,abs(v1h(elem,:)),'x--','Color',col(elem,:),'LineWidth',3);
    hold on
    hold on
end
title('Convergence of first 4\sigma d_1 terms of V_{1+}(\alpha_{n,1}) vs. h')
xlabel('n','FontSize',14)
ylabel('|V_{1+}(\alpha_{n,1})|','FontSize',14)
for n=1:M
    temp=sprintf('%3.1f', hrange(n));
    entry(n,:)=strcat('h=',temp,' ');
end
h=legend(entry);
set(h,'FontSize',12);

figure(2)
subplot(1,2,2)
for elem=1:M
    semilogy(nindex,abs(v2h(elem,:)),'x--','Color',col(elem,:),'LineWidth',3);
    hold on
    hold on
end
title('Convergence of first 4\sigma d_1 terms of V_{2+}(\alpha_{n,1}) vs. h')
xlabel('n','FontSize',14)
ylabel('|V_{2+}(\alpha_{n,1})|','FontSize',14)
for n=1:M
    temp=sprintf('%3.1f', hrange(n));
    entry(n,:)=strcat('h=',temp,' ');
end
h=legend(entry);
set(h,'FontSize',12);


