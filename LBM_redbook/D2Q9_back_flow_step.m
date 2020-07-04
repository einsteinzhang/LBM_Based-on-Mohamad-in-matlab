% Main.m
%LBM- 2-D2Q9, Flow in a channel with step, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear
nx=501;ny=81;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);
rho=ones(nx,ny);x=zeros(nx);y=zeros(ny);
Tm=zeros(nx);w(9)=zeros; Tvm=zeros(nx);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;
dx=1.0;dy=1.0;
xl=(nx-1)/(ny-1); yl=1.0;
dx=xl/(nx-1);
dy=yl/(ny-1);
x=(0:dx:xl);
y=(0:dy:yl);
uo=0.1;
alpha=0.01;
Re=uo*(ny-1)/alpha
omega=1./(3.*alpha+0.5);
count=0; tol=1.0e-4; error=10.;erso=0.0;
%setting velocity
for j=2:ny-1
u(1,j)=uo;
end
%Main Loop
while error>tol
% Collitions
[f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w);
% Streaming:
[f]=st(f);
% End of streaming
%Boundary condition:
[f]=boundary(nx,ny,f,uo,rho);
%Obsticale
[f]=obstc(nx,ny,f,uo,rho);
% Calculate rho, u, v
[rho,u,v]=ruv(nx,ny,f);
count=count+1;
ers=0.;
for i =1:nx
for j=1:ny
ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j);
end
end
error=abs(ers-erso);
erso=ers;
end
%Plotting data
result(nx,ny,x,y,u,v,uo,rho);



%Boudary conditions for Channel flow
function [f]=boundary(nx,ny,f,uo,rho)
    %right hand boundary
    for j=1:ny
    f(nx,j,3)=f(nx-1,j,3);
    f(nx,j,7)=f(nx-1,j,7);
    f(nx,j,6)=f(nx-1,j,6);
    end
    %bottom, and top boundary, bounce back
    for i=1:nx
    f(i,1,2)=f(i,1,4);
    f(i,1,5)=f(i,1,7);
    f(i,1,6)=f(i,1,8);
    f(i,ny,4)=f(i,ny,2);
    f(i,ny,7)=f(i,ny,5);
    f(i,ny,8)=f(i,ny,6);
    u(i,1)=0.0; v(i,1)=0.0;
    u(i,ny)=0.0; v(i,ny)=0.0;
    end
    %Left boundary, velocity is given= uo
    for j=2:ny-1
    f(1,j,1)=f(1,j,3)+2.*rho(1,j)*uo/3.;
    f(1,j,5)=f(1,j,7)-0.5*(f(1,j,2)-f(1,j,4))+rho(1,j)*uo/6.;
    f(1,j,8)=f(1,j,6)+0.5*(f(1,j,2)-f(1,j,4))+rho(1,j)*uo/6.;
    u(1,j)=uo; v(1,j)=0.0;
    end
% End of boundary conditions.
end
% Collition
function [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w)
    for j=1:ny
        for i=1:nx
            t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
            for k=1:9
                t2=u(i,j)*cx(k)+v(i,j)*cy(k);
                feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
                f(i,j,k)=(1.-omega)*f(i,j,k)+omega*feq(i,j,k);
            end
        end
    end
end
%Obsticale replace at the entrance, Back Fase Flow
function [f]=obstc(nx,ny,f,uo,rho)
    %length of obsticale= nx/4, heigh ny/2
    nxo=(nx-1)/4;
    nyl=(ny-1)/2;
    for i=1:nxo
        f(i,nyl,2)=f(i,nyl,4);
        f(i,nyl,5)=f(i,nyl,7);
        f(i,nyl,6)=f(i,nyl,8);
    end
    %bottom, and top boundary, bounce back
    for j=1:nyl
        f(nxo,j,1)=f(nxo,j,3);
    end
    for i=1:nxo
        for j=1:nyl
            u(i,j)=0.0;
            v(i,j)=0.0;
        end
    end
    % End
end
% Plots for channel flow
function result(nx,ny,x,y,u,v,uo,rho)
    nxo=(nx-1)/4;
    nyl=(ny-1)/2;
    for i=1:nxo
        for j=1:nyl
            u(i,j)=0.0;
            v(i,j)=0.0;
        end
    end
    for j=1:ny
        Tm1(j)=u(51,j)/uo;
        Tm2(j)=u(101,j)/uo;
        Tm3(j)=u(261,j)/uo;
        Tm4(j)=u(301,j)/uo;
        Tm5(j)=u(501,j)/uo;
    end
    figure
    plot(Tm1,y,Tm2,y,Tm3,y,Tm4,y,Tm5,y,'LineWidth',1.5)
    xlabel('U')
    ylabel('Y')
    %Stream function calculation
    for j=1:ny
        sx(:,j)=x(:);
    end
    for i=1:nx
        sy(i,:)=y(:);
    end
    str=zeros(nx,ny);
    for i=1:nx
        for j=2:ny
            str(i,j)=str(i,j-1)+0.5*(u(i,j)+u(i,j-1));
        end
    end
    figure
    contour(sx,sy,str)
end
function[rho,u,v]=ruv(nx,ny,f)
    rho=sum (f,3);
    for i=1:nx
    rho(i,ny)=f(i,ny,9)+f(i,ny,1)+f(i,ny,3)+2.*(f(i,ny,2)+f(i,ny,6)+f(i,ny,5));
    end
    %calculate velocity compnents
    u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) )./rho;
    v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) )./rho;
end
% Streaming:
function [f]=st(f)
    f(:,:,1)=circshift( squeeze(f(:,:,1)), [+1,+0] );
    f(:,:,2)=circshift( squeeze(f(:,:,2)), [+0,+1] );
    f(:,:,3)=circshift( squeeze(f(:,:,3)), [-1,+0] );
    f(:,:,4)=circshift( squeeze(f(:,:,4)), [+0,-1] );
    f(:,:,5)=circshift( squeeze(f(:,:,5)), [+1,+1] );
    f(:,:,6)=circshift( squeeze(f(:,:,6)), [-1,+1] );
    f(:,:,7)=circshift( squeeze(f(:,:,7)), [-1,-1] );
    f(:,:,8)=circshift( squeeze(f(:,:,8)), [+1,-1] );
end

% End of streaming
