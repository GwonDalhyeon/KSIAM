% function Collinocorner(L)
clc
clear all
%% Constants and Variables
L=5;
Tlength=0.5;
Xlength = 1;
Ylength = 1;
dx = 1/100;
dy = 1/100;
dt = 1/5000;
Nx = round(Xlength/dx)+2 ;
Ny = round(Ylength/dy)+2 ;
Nt = round(Tlength/dt)+1;
alpha = zeros(L,1);
beta = zeros(L,1);
C = zeros(L);
A = zeros(L,1);
U = zeros(Nx,Ny,Nt); % solution matrix
angle = zeros(L,1);
M =zeros(2*L,2*L); % matrix corresponding to linear system for obtaining auxiliary function on the semi-boundary node
M1 = zeros(L) ; 
M2 = zeros(L) ;
b = zeros(2*L,1,4);
angle = [pi/(2*L+1): pi/(2*L+1) : L*pi/(2*L+1)];
alpha = (cos(angle)).^2;
beta =  2/(2*L+1) * (sin(angle)).^2;

for l=1:L
    for m=1:L
          C(l,m) = beta(m)/(alpha(l)+alpha(m)-alpha(l)*alpha(m));
    end
    A(l) = 1+ sum(beta./alpha)+ alpha(l)*C(l,m);
end

% model problem 
syms x y t
% u_tt = u_xx + u_yy 
% u(x,y,0) = @(x,y) exp(i*(pi*x+pi*y)) , u_t(x,y,0) = -sqrt(2)*pi*i*exp(i*(pi*x+pi*y)) : initial condition
u = @(x,y) exp(i*(-pi*x-pi*y));
v = @(x,y,t) exp(sqrt(-1)*(-pi*x-pi*y-sqrt(2)*pi*t));
u_t = @(x,y) -sqrt(2)*pi*i*exp(i*(-pi*x-pi*y));
%  u = @(x,y) exp(i*(pi*x+pi*y));
%  v = @(x,y,t) exp(sqrt(-1)*(pi*x+pi*y-sqrt(2)*pi*t));
%  u_t = @(x,y) -sqrt(2)*pi*i*exp(i*(pi*x+pi*y));

[X,Y] = meshgrid(-dx/2:dx:Xlength+dx/2,-dy/2:dy:Ylength+dy/2); 

% Psi_corner_old = zeros(L,4,4);
% Psi_corner_new = zeros(L,2,4);

%% get initial data for u0, phi0 

U(:,:,1) = u(X,Y);
U(:,:,2) = u(X,Y) + dt* u_t(X,Y); 


%% updating u at interior nodes
for n=2:Nt-1
%     U(1,:,n)=v(0,Y(:,1),(n-1)*dt);
%     U(:,1,n)=v(X(1,:),0,(n-1)*dt);
%     U(Nx,:,n)=v(1,Y(:,1),(n-1)*dt);
%     U(:,Ny,n)=v(X(1,:),1,(n-1)*dt);
    
for i=2:Nx-1
    for j=2:Ny-1
U(i,j,n+1) = 2*U(i,j,n) - U(i,j,n-1) +(dt^2)/(dx^2)*(U(i+1,j,n)-2*U(i,j,n) +U(i-1,j,n)) +...
    +(dt^2)/(dy^2) * ( U(i,j+1,n)-2*U(i,j,n) +U(i,j-1,n)); 
    end
end




end
% 
% for n=1:50:Nt
%     
%     surf(X,Y,real(U(:,:,n)))
%     pause(0.5)
% end

% for n=1:250:Nt
%     
%     surf(X,Y,real(v(X,Y,(n-1)*dt)))
%     pause(0.5)
% end

%
L2error = zeros(Nt,2);
for n=1:Nt
  %L2error(:,:,n) = abs(U(:,:,n)-v(X,Y,(n-1)*dt));
   L2error(n,1) = sqrt(norm(V(:,:,n)-v(X,Y,(n-1)*dt),'fro')*dx*dy);
   
end

for n =2:Nt-1
  L2error(n+1,2) = L2error(n+1,1)/L2error(n,1);
end
