% Chapter 5, diffusion
% LBM- 1-D1Q3, diffusion equation
%2------0------1

clear
m=101;
%Ȩ��ϵ��
w0=4./6.;
w1=1./6.;
w2=w1;
%Cs^2
c2=1./3.;
dx=1.0;
rho=zeros(m);f0=zeros(m);f1=zeros(m);f2=zeros(m);
x=zeros(m);fluxq=zeros(m);
x(1)=0.0;
for i=1:m-1
x(i+1)=x(i)+dx;
end
alpha=0.25;
%�ɳ����ӣ�
omega=1/(3.*alpha+0.5);
twall=1.0;
nstep=200;
%��ʼ���ֲ�����
for i=1:m
    f0(i)=w0*rho(i);
    f1(i)=w1*rho(i);
    f2(i)=w1*rho(i);
end
%Collision:
for k1=1:nstep
    for i=1:m
        feq0=w0*rho(i);
        feq=w1*rho(i);
        f0(i)=(1-omega)*f0(i)+omega*feq0;
        f1(i)=(1-omega)*f1(i)+omega*feq;
        f2(i)=(1-omega)*f2(i)+omega*feq;
    end
% Streaming:
    for i=1:m-1
        f1(m-i+1)=f1(m-i);
        f2(i)=f2(i+1);
    end
    %Boundary condition:
    %��ߵ�һ��BC���ұߵڶ���BC��
    f1(1)=twall-f2(1)-f0(1);
    f1(m)=f1(m-1);
    f2(m)=f2(m-1);
    f0(m)=f0(m-1);
    for j=1:m
        rho(j)=f1(j)+f2(j)+f0(j);
    end
end
%Flux:
%��һ�ַ��������÷ֲ�������
for k=1:m
    flux(k)=omega*(f1(k)-f2(k))/c2;%q/k=omega/Cs^2*(f1-f2)
end
%�ڶ��ַ����������¶�ֵ��
for k=1:m-1
    fluxq(k)=rho(k)-rho(k+1);
end
fluxq(m)=fluxq(m-1);%��󽫾��ȶ˸�ֵ��
figure(1)
plot(x,rho)
title('Temperature, nstep=200')
xlabel('X')
ylabel('T')
figure(2)
plot(x,flux,'o',x,fluxq,'x')
title('Flux, time step=200')
xlabel('X')
ylabel('Flux')

