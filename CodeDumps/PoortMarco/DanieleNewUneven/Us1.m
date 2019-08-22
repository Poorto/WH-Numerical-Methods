function [Us1]=Us1(x,y,Alpha,sigma,d,A)
    mm=min(real(Alpha)); MM=max(real(Alpha)); N=length(Alpha); delta=(MM-mm)/N;
    Alpha=Alpha.'; Gamm=zeros(N,1);
    for n=1:N
        Gamm(n,1)=gam(Alpha(n,1),sigma);
    end
    
    Summer=delta*exp(-i*Alpha*x-i*Gamm*y)/(2*pi());
    Us1=A*Summer;
    
end