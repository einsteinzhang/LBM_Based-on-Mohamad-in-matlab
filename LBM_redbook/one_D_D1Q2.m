% Chapter 5 % LBM- 1-D, diffusion equation D1Q2 
% 2----------1

clear 
m=101; 
dx=1.0; 

f1=zeros(m);f2=zeros(m);%平衡分布函数
flux=zeros(m); 
x=zeros(m); 
x(1)=0.0; 
for i=1:m-1 
    x(i+1)=x(i)+dx;
end
alpha=0.25; 
omega=1/(alpha+0.5); 
twall=1.0; 
nstep=200; 

%初始化：
rho=zeros(m);%初始化温度
for i=1:m 
    f1(i)=0.5*rho(i); %初始化：分布函数feqi=wi*T
    f2(i)=0.5*rho(i); 
end
%Collision: 
for k1=1:nstep 
    for i=1:m 
        feq=0.5*rho(i); 
        f1(i)=(1-omega)*f1(i)+omega*feq; 
        f2(i)=(1-omega)*f2(i)+omega*feq; 
    end
% Streaming: 
    for i=1:m-1 
        f1(m-i+1)=f1(m-i); 
        f2(i)=f2(i+1); 
    end
    %Boundary condition: 
    %左边第一类BC，右边第二类BC。
    f1(1)=twall-f2(1); 
    f1(m)=f1(m-1); 
    f2(m)=f2(m-1); 
    for j=1:m
        rho(j)=f1(j)+f2(j);
    end
end
%Flux: 
for k=1:m 
    flux(k)=omega*(f1(k)-f2(k)); %q/k=omega/Cs^2*(f1-f2)
end
figure(1) 
plot(x,rho)
title('Temperature') 
xlabel('X') 
ylabel('T') 
figure(2) 
plot(x,flux,'o') 
title('Flux, time step=200')
xlabel('X') 
ylabel('Flux')