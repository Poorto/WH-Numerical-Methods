function [Us1i]=Us1i(x,y,Alpha,sigma,d,A)
    mm=min(real(Alpha)); MM=max(real(Alpha)); N=length(Alpha); delta=(MM-mm)/N;
    Gamm=zeros(1,N);
    for n=1:N
        Gamm(1,n)=gam(Alpha(1,n),sigma);
    end
    
    F=A.*exp(-i*Alpha*x-i*Gamm*y)/(2*pi());
    fun=@(t) interp1(real(Alpha),F,t,'spline');
    Us1i=integral(fun,mm,MM);


end