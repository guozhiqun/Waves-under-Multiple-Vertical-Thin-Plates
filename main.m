% This is the programme for calculating the transmission coefficient of
% water waves under multiple 2D vertical thin plates.
%
% Publication or report using this programme should refer to: 
% Yu, Y.F.; Guo, Z.Q.; Ma, Q.W. Transmission of Water Waves under Multiple Vertical Thin Plates. Water 2018, 10(4), 517.

clc,clear,close 'all';

%-----------------------------Parameter Input-----------------------------%
% Input the depth of the water
h=0.321;

% Input the number of plates
N=2;

% Input the draft of plates (size 1*N)
d=0.06*ones(1,N);

% Input the distance from the plates to the origin (size 1*N)
b=0.585*[0:N-1];

% Input the wavelength (size 1*nw)
ld = 1./(1/6.85:0.01:1/0.27);

% Input the truncation number
maxm= 100;
%------------------------End Parameter Input------------------------------%

g=9.80665;
nw=length(ld);
om=sqrt(2*pi*g./ld);

R0=zeros(N-1,nw);
T0=zeros(N-1,nw);
A0=zeros(N-1,nw);
B0=zeros(N-1,nw);
degA0=zeros(N-1,nw);
degB0=zeros(N-1,nw);
degR0=zeros(N-1,nw);
degT0=zeros(N-1,nw);

bf0=1;
for jj=2:N
    for i5=1:nw
        jp0=om(i5);
        v=jp0^2/g;
        F1=@(x)x*tanh(x*h)-v;
        
        k=1:maxm+1;
        K0=fzero(F1,1);
        k(1)=-1i*abs(K0);
        
        incr=pi/h;
        for i2=2:maxm+1
            xb=incr*(i2-1);
            xa=xb-0.5*incr;
            fa=eff(xa,h,jp0);
            fb=eff(xb,h,jp0);
            for i1=1:200
                xc=(xa+xb)/2;
                if abs(xa-xc)<10E-6
                    break
                end
                fc=eff(xc,h,jp0);
                if fc*fb<0
                    xa=xc;
                else
                    xb=xc;
                end
                fb=eff(xb,h,jp0);
            end
            k(i2)=xc;
        end
        
        if jj==2
            am= -(g^2)*(bf0^2)/(jp0^2);
            hd1=h-d(1);
            hd2=h-d(2);
            Fmn11=zeros(maxm+1);
            Fmn12=zeros(maxm+1);
            Fmn21=zeros(maxm+1);
            Fmn22=zeros(maxm+1);
            E=1:maxm+1;C=zeros(maxm+1);
            D=zeros(maxm+1);
            G=zeros(maxm+1);
            H=zeros(maxm+1);
            A=1:maxm+1;
            B=1:maxm+1;
            P=zeros(2*(maxm+1));
            Q=1:2*(maxm+1);
            for n=1:maxm+1
                for m=1:maxm+1
                    cs=2*cos(k(n)*h)*cos(k(m)*h);
                    nm=k(n)+k(m);
                    mn=k(n)-k(m);
                    sin11=sin(nm*hd1);
                    sin21=sin(mn*hd1);
                    sin31=sin(nm*h);
                    sin41=sin(mn*h);
                    if m==n
                        Fmn11(n,m)=am*(hd1/cs)*(1+sin11/(nm*hd1));
                        Fmn21(n,m)=am*(d(1)/cs)*(1+(sin31-sin11)/(nm*d(1)));
                    else
                        Fmn11(n,m)=am*(1/cs)*(sin11/nm+sin21/mn);
                        Fmn21(n,m)=am*(1/cs)*((sin31-sin11)/nm+(sin41-sin21)/mn);
                    end
                end
            end
            
            for n=1:maxm+1
                for m=1:maxm+1
                    cs=2*cos(k(n)*h)*cos(k(m)*h);
                    nm=k(n)+k(m);
                    mn=k(n)-k(m);
                    sin12=sin(nm*hd2);
                    sin22=sin(mn*hd2);
                    sin32=sin(nm*h);
                    sin42=sin(mn*h);
                    if m==n
                        Fmn12(n,m)=am*(hd2/cs)*(1+sin12/(nm*hd2));
                        Fmn22(n,m)=am*(d(2)/cs)*(1+(sin32-sin12)/(nm*d(2)));
                    else
                        Fmn12(n,m)=am*(1/cs)*(sin12/nm+sin22/mn);
                        Fmn22(n,m)=am*(1/cs)*((sin32-sin12)/nm+(sin42-sin22)/mn);
                    end
                end
            end
            
            for m=1:maxm+1
                E(m)=k(1)*Fmn11(1,m);
                for n=1:maxm+1
                    C(n,m)=k(1)*Fmn11(n,m)+k(n)*Fmn21(n,m);
                    D(n,m)=-k(n)*exp(-k(n)*b(2))*Fmn21(n,m);
                    G(n,m)=k(n)*exp(-k(n)*b(2))*Fmn22(n,m);
                    H(n,m)=k(1)*Fmn12(n,m)-k(n)*Fmn22(n,m);
                end
            end
            
            for n=1:maxm+1
                Q(n)=E(n);
                for m=1:maxm+1
                    P(m,n)=C(n,m);
                    P(m+maxm+1,n)=G(n,m);
                    P(m,n+maxm+1)=D(n,m);
                    P(m+maxm+1,n+maxm+1)=H(n,m);
                end
            end
            
            for n=maxm+2:2*(maxm+1)
                Q(n)=0;
            end
            X=P\Q.';
            for n=1:maxm+1
                A(1,n)=X(n);
                B(1,n)=X(n+maxm+1);
            end
            r0=1-A(1,1)+B(1,1)*exp(-k(1)*b(2));
            t0=A(1,1)*exp(-k(1)*b(2))-B(1,1);
            
        else
            A=zeros(jj-1,maxm+1);
            B=zeros(jj-1,maxm+1);
            P=zeros(2*(jj-1)*(maxm+1),2*(jj-1)*(maxm+1));
            Q=zeros(1,2*(jj-1)*(maxm+1));
            Fmn11=zeros(maxm+1,maxm+1,jj);
            Fmn21=zeros(maxm+1,maxm+1,jj);
            am= -(g^2)*(bf0^2)/(jp0^2);
            hd=h-d;
            for i1=1:maxm+1
                for i2=1:maxm+1
                    for i3=1:jj
                        cs=2*cos(k(i1)*h)*cos(k(i2)*h);
                        nm=k(i1)+k(i2);
                        mn=k(i1)-k(i2);
                        sin11=sin(nm*hd(i3));
                        sin21=sin(mn*hd(i3));
                        sin31=sin(nm*h);
                        sin41=sin(mn*h);
                        if i2==i1
                            Fmn11(i1,i2,i3)=am*(hd(i3)/cs)*(1+sin11/(nm*hd(i3)));
                            Fmn21(i1,i2,i3)=am*(d(i3)/cs)*(1+(sin31-sin11)/(nm*d(i3)));
                        else
                            Fmn11(i1,i2,i3)=am*(1/cs)*(sin11/nm+sin21/mn);
                            Fmn21(i1,i2,i3)=am*(1/cs)*((sin31-sin11)/nm+(sin41-sin21)/mn);
                        end
                    end
                end
            end
            
            for j=2:(jj-1)
                P(j-1,(maxm+1)*(j-2)+1)=-exp(k(1)*(b(j-1)-b(j)));
                P(j-1,(maxm+1)*(j-1)+1)=1;
                P(j-1,(maxm+1)*(jj-1)+(maxm+1)*(j-2)+1)=1;
                P(j-1,(maxm+1)*(jj-1)+(maxm+1)*(j-1)+1)=-exp(-k(1)*(b(j+1)-b(j)));
                for j2=1:(maxm+1)
                    c1= Fmn21(j2,1,1);
                    c2= Fmn21(j2,1,j);
                    c3= Fmn11(j2,1,1);
                    c4= Fmn11(j2,1,j);
                    c8=(jj-1)*(maxm+1)+(maxm+1)*(j-2)+j2;
                    P(jj-2+maxm*(jj-2)+j2,1)=k(1)*c1+1i*k(1)*c3;
                    P(jj-2+maxm*(jj-2)+j2,(maxm+1)*(jj-1)+1)=-k(1)*c1*exp(-k(1)*b(2));
                    P(c8,(maxm+1)*(j-2)+1)=(k(1)*c2+1i*k(1)*c4)*exp(k(1)*(b(j-1)-b(j)));
                    P(c8,(maxm+1)*(j-1)+1)=-1i*k(1)*c4;
                    P(c8,(maxm+1)*(jj-1)+(maxm+1)*(j-2)+1)=-k(1)*c2+1i*k(1)*c4;
                    P(c8,(maxm+1)*(jj-1)+(maxm+1)*(j-1)+1)=-1i*k(1)*c4*exp(-k(1)*(b(j+1)-b(j)));
                    P((jj-1)*(maxm+1)+(maxm+1)*(jj-2)+j2,(maxm+1)*(jj-2)+1)=k(1)*Fmn21(j2,1,jj)*exp(k(1)*(b(jj-1)-b(jj)));
                    P((jj-1)*(maxm+1)+(maxm+1)*(jj-2)+j2,(maxm+1)*(jj-1)+(maxm+1)*(jj-2)+1)=-k(1)*Fmn21(j2,1,jj)+1i*k(1)* Fmn11(j2,1,jj);
                    Q(jj-2+maxm*(jj-2)+j2)=1i*k(1)*c3;
                    for j1=2:(maxm+1)
                        c5=Fmn21(j2,j1,j);
                        c6=Fmn11(j2,j1,j);
                        c7=Fmn21(j2,j1,1);
                        c9=jj-2+maxm*(j-2)+j1-1;
                        P(c9,(maxm+1)*(j-2)+j1)=exp(k(j1)*(b(j-1)-b(j)));
                        P(c9,(maxm+1)*(j-1)+j1)=-1;
                        P(c9,(maxm+1)*(jj-1)+(maxm+1)*(j-2)+j1)=-1;
                        P(c9,(maxm+1)*(jj-1)+(maxm+1)*(j-1)+j1)=exp(-k(j1)*(b(j+1)-b(j)));
                        P(jj-2+maxm*(jj-2)+j2,j1)=k(j1)*c7+1i*k(1)*Fmn11(j2,j1,1);
                        P(jj-2+maxm*(jj-2)+j2,(maxm+1)*(jj-1)+j1)=-c7*exp(-k(j1)*b(2));
                        P(c8,(maxm+1)*(j-2)+j1)=(k(j1)*c5+1i*k(1)*c6)*exp(k(j1)*(b(j-1)-b(j)));
                        P(c8,(maxm+1)*(j-1)+j1)=-1i*k(1)*c6;
                        P(c8,(maxm+1)*(jj-1)+(maxm+1)*(j-2)+j1)=1i*k(1)*c6-k(j1)*c5;
                        P(c8,(maxm+1)*(jj-1)+(maxm+1)*(j-1)+j1)=-1i*k(1)*c6*exp(-k(j1)*(b(j+1)-b(j)));
                        P((jj-1)*(maxm+1)+(maxm+1)*(jj-2)+j2,(maxm+1)*(jj-2)+j1)=k(j1)*Fmn21(j2,j1,jj)*exp(k(j1)*(b(jj-1)-b(jj)));
                        P((jj-1)*(maxm+1)+(maxm+1)*(jj-2)+j2,(maxm+1)*(jj-1)+(maxm+1)*(jj-2)+j1)=-1i*k(1)*Fmn11(j2,j1,jj)-k(j1)*Fmn21(j2,j1,jj);
                    end
                end
            end
            
            X=P\Q.';
            
            for j=1:(jj-1)
                A(j,:)=X((j-1)*(maxm+1)+1:(maxm+1)*j).';
                B(j,:)=X((maxm+1)*(jj-1)+(j-1)*(maxm+1)+1:(maxm+1)*(jj-1)+(maxm+1)*j).';
            end
            
            r0=1-A(1,1)+B(1,1)*exp(-k(1)*b(2));
            t0=A(jj-1,1)*exp(k(1)*(b(jj-1)-b(jj)))-B(jj-1,1);
        end
        
        %Results
        R0(jj-1,i5)=abs(r0);         % abs(R_0)
        T0(jj-1,i5)=abs(t0);         % abs(T_0)
        A0(jj-1,i5)=abs(A(1,1));     % abs(A_0)
        B0(jj-1,i5)=abs(B(1,1));     % abs(B_0)
        degR0(jj-1,i5)=angle(r0);    % angle(R_0)
        degT0(jj-1,i5)=angle(t0);    % angle(T_0£©
        degA0(jj-1,i5)=angle(A(1,1));    % angle(A_0£©
        degB0(jj-1,i5)=angle(B(1,1));    % angle(B_0)
    end
end

% Plot module of T0
bl=b(2)./ld;                            
T0=T0(N-1,:);                          
plot(bl,T0);
title(['N=',num2str(N)]);
xlabel('b/¦Ë');
ylabel('|T_0|');
box off;

% Plot phase of T0
figure;
degT0=degT0(N-1,:)*180/pi;
plot(bl,degT0);
title(['N=',num2str(N)]);
xlabel('b/¦Ë');
ylabel('¡ÏT_0');
box off;
