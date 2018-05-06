%mathematical model of a linear engine with a classical method of stacking a stator winding
%number of inductor grooves (z=6)
function lad_z12_6_sprint
    %output data of the asynchronous motor
Rb=0.1001*10^7;
rs=9.5;
Ls=0.037;
rr=4.6345*10^-5;
Lr=0.0372*10^-5;
dt=0.001;
tz=9.769*10^-3;
m=1.9;
v0=0;
wn=200;
f=50;
w=2*pi*f;
UA=wn/dt;
Um=310/2;
X=zeros(17,1);
F=0;
K=input('Duration of the cycle k=');
for k=1:(K+1)
    v(1,k)=v0;     %create vector string for speed graphing
    f(1,k)=sum(F); %create vector strings for graphing force
    
    Uab=Um*cos(w*(k-1)*dt+2*pi/3);
    Ubc=Um*cos(w*(k-1)*dt);
    
    %the formation of the matrix A.
   A=zeros(17);
B=2*Rb*(rr+Lr/dt)+1/dt;
B1=6*Rb*(rr+Lr/dt)+(-4*Rb)*Lr*v0/(2*tz)+1/dt;
B2=55*Rb*(rr+Lr/dt)+(-45*Rb)*Lr*v0/(2*tz)+1/dt;
B3=550*Rb*(rr+Lr/dt)+(-450*Rb)*Lr*v0/(2*tz)+1/dt;
B4=1000*Rb*(rr+Lr/dt)+1/dt;
B5=550*Rb*(rr+Lr/dt)+450*Rb*Lr*v0/(2*tz)+1/dt;
B6=55*Rb*(rr+Lr/dt)+(45*Rb)*Lr*v0/(2*tz)+1/dt;
B7=6*Rb*(rr+Lr/dt)+(4*Rb)*Lr*v0/(2*tz)+1/dt;

C=-Rb*(rr+Lr/dt)+(2*Rb*Lr+1)*v0/(2*tz);
C1=-Rb*(rr+Lr/dt)+(6*Rb*Lr+1)*v0/(2*tz);
C2=-5*Rb*(rr+Lr/dt)+(55*Rb*Lr+1)*v0/(2*tz);
C3=-50*Rb*(rr+Lr/dt)+(550*Rb*Lr+1)*v0/(2*tz);
C4=-500*Rb*(rr+Lr/dt)+(1000*Rb*Lr+1)*v0/(2*tz);
C5=-500*Rb*(rr+Lr/dt)+(550*Rb*Lr+1)*v0/(2*tz);
C6=-50*Rb*(rr+Lr/dt)+(55*Rb*Lr+1)*v0/(2*tz);
C7=-5*Rb*(rr+Lr/dt)+(6*Rb*Lr+1)*v0/(2*tz);

D=-Rb*Lr*v0/(2*tz);
D1=5*D;
D2=50*D;
D3=500*D;

E=-Rb*(rr+Lr/dt)-(2*Rb*Lr+1)*v0/(2*tz);
E1=-5*Rb*(rr+Lr/dt)-(6*Rb*Lr+1)*v0/(2*tz);
E2=-50*Rb*(rr+Lr/dt)-(55*Rb*Lr+1)*v0/(2*tz);
E3=-500*Rb*(rr+Lr/dt)-(550*Rb*Lr+1)*v0/(2*tz);
E4=-500*Rb*(rr+Lr/dt)-(1000*Rb*Lr+1)*v0/(2*tz);
E5=-50*Rb*(rr+Lr/dt)-(550*Rb*Lr+1)*v0/(2*tz);
E6=-5*Rb*(rr+Lr/dt)-(55*Rb*Lr+1)*v0/(2*tz);
E7=-Rb*(rr+Lr/dt)-(6*Rb*Lr+1)*v0/(2*tz);

T=-wn*Lr*v0/(2*tz);

Y=-wn*(rr+Lr/dt);

W1=-wn*Lr/dt;
P=-Rb*Lr/dt;
Q=(2*Rb*Lr+1)/dt;
KS=rs+Ls/dt;

Q1=(6*Rb*Lr+1)/dt;
Q2=(55*Rb*Lr+1)/dt;
Q3=(550*Rb*Lr+1)/dt;
Q4=(1000*Rb*Lr+1)/dt;
    for n=1:3
        A(n+3,n+14)=(-1)^(n+1)*T;
        A(n+4,n+14)=(-1)^(n+1)*Y;
        A(n+5,n+14)=(-1)^n*T;
        A(n+6,n+14)=(-1)^n*T;
        A(n+7,n+14)=(-1)^n*Y;
        A(n+8,n+14)=(-1)^(n+1)*T;
        A(17,n+14)=1;
    end;
    
    for n=1:2
        A(n+14,n+14)=(-1)^(n+1)*KS;
        A(n+14,17)=(-1)^n*KS;
        A(n+14,n+14)=UA;
        A(n+14,7)=(-1)^n*UA;
        A(n+14,n+7)=-UA;
        A(n+14,10)=(-1)^(n+1)*UA;
    end;
    
    for n=1:6
        A(n+4,n+4)=B;
        A(n+5,n+4)=E;
        A(n+3,n+4)=C;
    end;
    
    for n=1:7
        A(n+2,n+4)=D;
        A(n+5,n+3)=-D;
    end;
    
    A(1,1)=B4;
    A(1,2)=C5;
    A(1,3)=D2;
    A(2,1)=E4;
    A(2,2)=B5;
    B(2,3)=C6;
    A(2,4)=D1;
    A(3,1)=-D3;
    A(3,2)=E5;
    A(3,3)=B6;
    A(3,4)=C7;
    A(4,2)=-D2;
    A(4,3)=E6;
    A(4,4)=B7;
    A(5,3)=-D1;
    A(5,4)=E7;
    A(10,11)=C1;
    A(10,12)=D1;
    A(11,11)=B1;
    A(11,12)=C2;
    A(11,13)=D2;
    A(12,11)=E1;
    A(12,12)=B2;
    A(12,13)=C3;
    A(12,13)=D3;
    A(13,11)=-D1;
    A(13,12)=E2;
    A(13,13)=B3;
    A(13,14)=C4;
    A(14,12)=-D2;
    A(14,13)=E3;
    A(14,14)=B4;
    
    %matrix of free members
    
    S=[Q4*X(1)+P*(500*X(2)); %1
        Q3*X(2)+P*(500*X(1)+50*X(3)); %2
        Q2*X(3)+P*(50*X(2)+5*X(4)); %3
        Q1*X(4)+P*(5*X(3)+X(5)); %4
        W1*X(15)+Q*X(5)+P*(X(4)+X(6)); %5
        (-1)*W1*X(16)+Q*X(6)+P*(X(5)+X(7)); %6
        W1*X(17)+Q*X(7)+P*(X(6)+X(8)); %7
        (-1)*W1*X(15)+Q*X(8)+P*(X(7)+X(9)); %8
        W1*X(16)+Q*X(9)+P*(X(8)+X(10)); %9
        (-1)*W1*X(17)+Q*X(10)+P*(X(9)+X(11)); %10
        Q1*X(11)+P*(X(10)+5*X(12)); %11
        Q2*X(12)+P*(5*X(11)+50*X(13)); %12
        Q3*X(13)+P*(50*X(12)+500*X(14)); %13
        Q4*X(14)+P*500*X(13); %14
        UA*(X(5)-X(8)-X(7)+X(10))+Ls/dt*(X(15)-X(17))+Uab; %15
        UA*(X(7)-X(10)+X(6)-X(9))+Ls/dt*(X(17)-X(16))+Ubc; %16
        0];

    %solution by the Gauss-Jordan method
    
Z=rref([A S]); %transform the expanded matrix into triangular form
X=Z(1:17,18:18); %selecting the last column of the matrix

% current in the rotor

Ir=[1000*Rb*X(1)-Rb*(500*X(2)); %1
    550*Rb*X(2)-Rb*(500*X(1)+50*X(3)); %2
    55*Rb*X(3)-Rb*(50*X(2)+5*X(4)); %3
    6*Rb*X(4)-Rb*(5*X(3)+X(5)); %4
    -wn*X(15)+2*Rb*X(5)-Rb*(X(4)+X(6)); %5
    (-1)*(-wn)*X(16)+2*Rb*X(6)-Rb*(X(5)+X(7)); %6
    -wn*X(17)+2*Rb*X(7)-Rb*(X(6)+X(8)); %7
    (-1)*(-wn)*X(15)+2*Rb*X(8)-Rb*(X(7)+X(9)); %8
    -wn*X(16)+2*Rb*X(9)-Rb*(X(8)+X(10)); %9;
    (-1)*(-wn)*X(17)*2*Rb*X(10)-Rb*(X(9)+X(11)); %10
    6*Rb*X(11)-Rb*(X(10)+5*X(12)); %11
    55*Rb*X(12)-Rb*(5*X(11)+50*X(13)); %12
    550*Rb*X(13)-Rb*(50*X(12)+500*X(14)); %13
    1000*Rb*X(14)-Rb*(500*X(13))]; %14

%electromagnetic force

F(1)=X(2)*Ir(1)/(2*tz);
for n=1:12
    F(n+1)=(X(n+2)-X(n))*Ir(n+1)/(2*tz);
end;
    F(14)=-X(13)*Ir(14)/(2*tz);
    
    %speed
    
    v0=v0+(sum(f)/m)*dt;
end;

%drawing up charts

k=0:K;
subplot(2,1,1);
plot(k*dt,v);
title('Speed');
xlabel('t, c');
ylabel('v, m/s');
grid on;
subplot(2,1,2);
plot(k*dt,f);
title('electromagnetic force');
xlabel('t,c');
ylabel('F,H');
grid on;
end
        
