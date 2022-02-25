function [phi,phip,phipp]=Osplinepp(K,N0)
%Order K>3
%N0 Number of samples per cycle
delt=1/N0;
unit=(0:delt:1-delt)';
%Order K
knots=[flip(-[1:K]) [1:K]];
P=ones(N0,K+1);
Pp=ones(N0,K+1);
Ppp=ones(N0,K+1);
ir0=0;
for nint=1:K+1
    u=-(K+1)/2-1+nint+unit;
    coef= [1, -knots(ir0+1)]/(-knots(ir0+1));
    for k=1:K
        P(:,nint)=P(:,nint).*(u-knots(ir0+k))/(-knots(ir0+k));
        if k>1 coef = conv(coef, [1, -knots(ir0+k)]/(-knots(ir0+k)));
        end
    end
    %Derivative Polynomial coefficients
    C(:,nint)=coef';
    Cp(:,nint)=polyder(coef)';
    Cpp(:,nint)=polyder(Cp(:,nint))';
    %First derivative by Horner scheme
    Pp(:,nint)=Cp(1,nint)*Pp(:,nint);
    for k=1:K-1
        Pp(:,nint)=Pp(:,nint).*u + Cp(k+1,nint);
    end;
    %Second derivative by Horner scheme
    Ppp(:,nint)=Cpp(1,nint)*Ppp(:,nint);
    for k=1:K-2
        Ppp(:,nint)=Ppp(:,nint).*u + Cpp(k+1,nint);
    end;
    ir0=ir0+1;
end;
N=(K+1)*N0;
phi=reshape(P,[N,1]);
phip=reshape(Pp,[N,1]);
phipp=reshape(Ppp,[N,1]);







