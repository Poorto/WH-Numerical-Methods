function [Solm]=Solm(A,h,sigma,dvect,phi)
    N=round(A/h);
    
    D=zeros(6*N+3);
    K=zeros(6*N+3);
    S=zeros(6*N+3,1);
    
    for m=-N:N;
        r=3*(m+N);
        
        x=t(m*h);
        Gt=G(x,sigma,dvect);
        
        St=F0(x,sigma,dvect,phi);
        for j=1:3
            S(r+j,1)=St(j,1);
        end
        
        for j=1:3
            for k=1:3
                D(r+j,r+k)=Gt(j,k);
            end
        end
    
        
        for n=-N:N
            c=3*(n+N);
            y=t(n*h);
            Tp=tp();
            Mt=Mm(x,y,sigma,dvect);
            Bt=h*Tp*Mt/(2*pi()*i);
            
            for j=1:3
                for k=1:3
                    K(r+j,c+k)=Bt(j,k);
                end
            end
            
        end
    end
    Solm=inv(D+K)*S;
end