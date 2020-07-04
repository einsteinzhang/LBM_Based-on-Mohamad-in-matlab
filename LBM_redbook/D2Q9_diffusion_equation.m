%LBM- 2-D2Q9, advection-diffusion equation, note that c2=1/3, w0=4/9, w1-w4, 1/9
% and w5-w8, 1/36
clear
m=101;n=101;
xl=1.0;yl=1.0;
w0=4./9.;
c2=1./3.;
dx=1.0;dy=1.0;
f0=zeros(m,n);f=zeros(m,n,8);feq=zeros(m,n,8);f0eq=zeros(m,n);
rho=zeros(m,n);x=zeros(m);y=zeros(n);fluxq=zeros(m);flux=zeros(m);
Tm=zeros(m);Z=zeros(n,m);w(8)=zeros;
x(1)=0.0; y(1)=0.0;
for i=1:m-1
    x(i+1)=x(i)+dx;
end
for j=1:n-1
    y(j+1)=y(j)+dy;
end
for k=1:4
    w(k)=1./9.;
end
for k=5:8
    w(k)=1./36.;
end
alpha=1.00;
u=0.1;v=0.4;
omega=1./(3.*alpha+0.5);
twall=1.0;
nstep=400;
%Collision:
for kk=1:nstep
    for j=1:n
        for i=1:m
            f0eq(i,j)=w0*rho(i,j);
            f0(i,j)=(1.-omega)*f0(i,j)+omega*f0eq(i,j);
            feq(i,j,1)=w(1)*rho(i,j)*(1.+3.*u);
            feq(i,j,2)=w(2)*rho(i,j)*(1.+3.*v);
            feq(i,j,3)=w(3)*rho(i,j)*(1.-3.*u);
            feq(i,j,4)=w(4)*rho(i,j)*(1.-3.*v);
            feq(i,j,5)=w(5)*rho(i,j)*(1.+3.*(u+v));
            feq(i,j,6)=w(6)*rho(i,j)*(1.+3.*(v-u));
            feq(i,j,7)=w(7)*rho(i,j)*(1.-3.*(u+v));
            feq(i,j,8)=w(8)*rho(i,j)*(1.+3.*(u-v));
            for k=1:8
                f(i,j,k)=(1.-omega)*f(i,j,k)+omega*feq(i,j,k);
            end
        end
    end
    % Streaming:
    f(:,:,1)=circshift( squeeze(f(:,:,1)), [+1,+0] );
    f(:,:,2)=circshift( squeeze(f(:,:,2)), [+0,+1] );
    f(:,:,3)=circshift( squeeze(f(:,:,3)), [-1,+0] );
    f(:,:,4)=circshift( squeeze(f(:,:,4)), [+0,-1] );
    f(:,:,5)=circshift( squeeze(f(:,:,5)), [+1,+1] );
    f(:,:,6)=circshift( squeeze(f(:,:,6)), [-1,+1] );
    f(:,:,7)=circshift( squeeze(f(:,:,7)), [-1,-1] );
    f(:,:,8)=circshift( squeeze(f(:,:,8)), [+1,-1] );
    % End of streaming
    %Boundary condition:
    %left boundary, twall=1.0
    for j=1:n
        f(1,j,1)=w(1)*twall+w(3)*twall-f(1,j,3);
        f(1,j,5)=w(5)*twall+w(7)*twall-f(1,j,7);
        f(1,j,8)=w(8)*twall+w(6)*twall-f(1,j,6);
    end
    %bottom boundary, T=0.0
    for i=1:m
        f(i,1,2)=-f(i,1,4);
        f(i,1,5)=-f(i,1,7);
        f(i,1,6)=-f(i,1,8);
    end
    %Top boundary, T=0.0
    for i=1:n
        f(i,m,7)=-f(i,m,5);
        f(i,m,4)=-f(i,m,2);
        f(i,m,8)=-f(i,m,6);
    end
    %right hand boundary
    for j=1:m
        f(n,j,3)=-f(n,j,1);
        f(n,j,7)=-f(n,j,5);
        f(n,j,6)=-f(n,j,8);
    end
% End of boundary conditions
    for j=1:n
        for i=1:m
            sumk=0.0;
            for k=1:8
                sumk=sumk+f(i,j,k);
            end
            rho(i,j)=f0(i,j)+sumk;
        end
    end
end
%rotating matrix for contour plotting
for j=1:n
    for i=1:m
        Z(j,i)=rho(i,j);
    end
end
%�м�y=0.5���¶�
for i=1:n
    Tm(i)=rho(i,(n-1)/2);
end
figure(1)
plot(x,Tm, 'LineWidth',2)
xlabel('X')
ylabel('T')
figure(2)
contour(Z,'showText','on','LineWidth',2)
% title(��Flux��)
% xlabel(��X��)
% ylabel(��Flux��)
