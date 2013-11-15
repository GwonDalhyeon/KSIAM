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

Psi_corner_old = zeros(L,4,4);
Psi_corner_new = zeros(L,2,4);

%% get initial data for u0, phi0 

U(:,:,1) = u(X,Y);
U(:,:,2) = u(X,Y) + dt* u_t(X,Y); 

Phi_1 = zeros(L,Ny,2); % phi_1(l.ny,1) : 'n-1' time, phi_1(l.ny,2) : 'n' time   
Phi_2 = zeros(L,Nx,2);
Phi_3 = zeros(L,Ny,2); % phi_1(l.ny,1) : 'n-1' time, phi_1(l.ny,2) : 'n' time   
Phi_4 = zeros(L,Nx,2);
Psi_1 = zeros(L,Ny,2);
Psi_2 = zeros(L,Nx,2);
Psi_3 = zeros(L,Ny,2);
Psi_4 = zeros(L,Nx,2);

for l=1:L
    % Psi for initial data
    Psi_1(l,:,1) = alpha(l)*Phi_1(l,:,1) + 1/2 * (U(1,:,1)+U(2,:,1));
    Psi_2(l,:,1) = alpha(l)*Phi_2(l,:,1) + 1/2 * transpose(U(:,1,1)+U(:,2,1));
    Psi_3(l,:,1) = alpha(l)*Phi_3(l,:,1) + 1/2 * (U(Nx-1,:,1)+U(Nx,:,1));
    Psi_4(l,:,1) = alpha(l)*Phi_4(l,:,1) + 1/2 * transpose(U(:,Ny-1,1)+U(:,Ny,1));
    
    % Psi for time n=2 data
    Psi_1(l,:,2) = alpha(l)*Phi_1(l,:,2) + 1/2 * (U(1,:,2)+U(2,:,2));
    Psi_2(l,:,2) = alpha(l)*Phi_2(l,:,2) + 1/2 * transpose(U(:,1,2)+U(:,2,2));
    Psi_3(l,:,2) = alpha(l)*Phi_3(l,:,2) + 1/2 * (U(Nx-1,:,2)+U(Nx,:,2));
    Psi_4(l,:,2) = alpha(l)*Phi_4(l,:,2) + 1/2 * transpose(U(:,Ny-1,2)+U(:,Ny,2));
end


%% updating u at interior nodes
for n=2:Nt-1
    
for i=2:Nx-1
    for j=2:Ny-1
U(i,j,n+1) = 2*U(i,j,n) - U(i,j,n-1) +(dt^2)/(dx^2)*(U(i+1,j,n)-2*U(i,j,n) +U(i-1,j,n)) +...
    +(dt^2)/(dy^2) * ( U(i,j+1,n)-2*U(i,j,n) +U(i,j-1,n)); % at Nx-grid, the value of u is not changing. 
    end
end


%% update phi 


for i= 2: Nx-1

    Phi_2(:,i,1) = 2* Phi_2(:,i,2) - Phi_2(:,i,1) + (dt)^2/(dx)^2 * (Phi_2(:,i-1,2) - 2*Phi_2(:,i,2) + Phi_2(:,i+1,2)) ;
    Phi_4(:,i,1) = 2* Phi_4(:,i,2) - Phi_4(:,i,1) + (dt)^2/(dx)^2 * (Phi_4(:,i-1,2) - 2*Phi_4(:,i,2) + Phi_4(:,i+1,2)) ;

end

for j= Ny-1
    Phi_1(:,j,1) = 2* Phi_1(:,j,2) - Phi_1(:,j,1) + (dt)^2/(dx)^2 * (Phi_1(:,j-1,2) - 2*Phi_1(:,j,2) + Phi_1(:,j+1,2)) ;
    Phi_3(:,j,1) = 2* Phi_3(:,j,2) - Phi_3(:,j,1) + (dt)^2/(dx)^2 * (Phi_3(:,j-1,2) - 2*Phi_3(:,j,2) + Phi_3(:,j+1,2)) ;
end
%Phi(:,:,1) - n+1 time
%Phi(:,:,2) - n time
% Phi_1의 오리지널 데이터(n-1,n)을 n,n+1 타임에서의 데이터로 변환 

Phi_1temp = zeros(L,Ny);
Phi_2temp = zeros(L,Nx);
Phi_3temp = zeros(L,Ny);
Phi_4temp = zeros(L,Nx);

Phi_1temp = Phi_1(:,:,2) ;
Phi_1(:,:,2) = Phi_1(:,:,1);
Phi_1(:,:,1) = Phi_1temp;

Phi_2temp = Phi_2(:,:,2) ;
Phi_2(:,:,2) = Phi_2(:,:,1);
Phi_2(:,:,1) = Phi_2temp;

Phi_3temp = Phi_3(:,:,2) ;
Phi_3(:,:,2) = Phi_3(:,:,1);
Phi_3(:,:,1) = Phi_3temp;

Phi_4temp = Phi_4(:,:,2) ;
Phi_4(:,:,2) = Phi_4(:,:,1);
Phi_4(:,:,1) = Phi_4temp;

%% update u at boundary nodes

tempsum1 = zeros(1,Ny);
tempsum2 = zeros(1,Nx);
tempsum3 = zeros(1,Ny);
tempsum4 = zeros(1,Nx);

for l=1:L
    tempsum1 = tempsum1 + beta(l)* (Phi_1(l,:,2)-Phi_1(l,:,1));
    tempsum2 = tempsum2 + beta(l)* (Phi_2(l,:,2)-Phi_2(l,:,1));
    tempsum3 = tempsum3 + beta(l)* (Phi_3(l,:,2)-Phi_3(l,:,1));
    tempsum4 = tempsum4 + beta(l)* (Phi_4(l,:,2)-Phi_4(l,:,1));
end

U(:,1,n+1) = (2*dt*dy)/(dt+dy) *((1/(2*dy) - 1/(2*dt))*U(:,2,n+1) + (1/(2*dt) + 1/(2*dy))*U(:,2,n) + (1/(2*dt) - 1/(2*dy))*U(:,1,n)) +...
    +(2*dy)/(dt+dy) * transpose(tempsum1);

U(1,:,n+1) = (2*dt*dx)/(dt+dx) *((1/(2*dx) - 1/(2*dt))*U(2,:,n+1) + (1/(2*dt) + 1/(2*dx))*U(2,:,n) + (1/(2*dt) - 1/(2*dx))*U(1,:,n)) +...
    +(2*dx)/(dt+dx) *tempsum2;

U(:,Ny,n+1) = (2*dt*dy)/(dt+dy) *((1/(2*dy) - 1/(2*dt))*U(:,Ny-1,n+1) + (1/(2*dt) + 1/(2*dy))*U(:,Ny-1,n) + (1/(2*dt) - 1/(2*dy))*U(:,Ny,n)) +...
    +(2*dy)/(dt+dy) * transpose(tempsum3);

U(Nx,:,n+1) = (2*dt*dx)/(dt+dx) *((1/(2*dx) - 1/(2*dt))*U(Nx-1,:,n+1) + (1/(2*dt) + 1/(2*dx))*U(Nx-1,:,n) + (1/(2*dt) - 1/(2*dx))*U(Nx,:,n)) +...
    +(2*dx)/(dt+dx) *tempsum4;

%% update psi except at corner nodes

for l=1:L
        
    Psi_1(l,:,1) = alpha(l)*Phi_1(l,:,2) + 1/2 * (U(1,:,2)+U(2,:,2));
    Psi_2(l,:,1) = alpha(l)*Phi_2(l,:,2) + 1/2 * transpose(U(:,1,2)+U(:,2,2));
    Psi_3(l,:,1) = alpha(l)*Phi_3(l,:,2) + 1/2 * (U(Nx,:,2)+U(Nx-1,:,2));
    Psi_4(l,:,1) = alpha(l)*Phi_4(l,:,2) + 1/2 * transpose(U(:,Ny,2)+U(:,Ny-1,2));
end

%Psi(:,:,1) - n+1 time
%Psi(:,:,2) - n time
% Psi_1의 오리지널 데이터(n-1,n)을 n,n+1 타임에서의 데이터로 변환 

Psi_1temp = zeros(L,Ny);
Psi_2temp = zeros(L,Nx);
Psi_3temp = zeros(L,Ny);
Psi_4temp = zeros(L,Nx);

Psi_1temp = Psi_1(:,:,2) ;
Psi_1(:,:,2) = Psi_1(:,:,1);
Psi_1(:,:,1) = Psi_1temp;

Psi_2temp = Psi_2(:,:,2) ;
Psi_2(:,:,2) = Psi_2(:,:,1);
Psi_2(:,:,1) = Psi_2temp;

Psi_3temp = Psi_3(:,:,2) ;
Psi_3(:,:,2) = Psi_3(:,:,1);
Psi_3(:,:,1) = Psi_3temp;

Psi_4temp = Psi_4(:,:,2) ;
Psi_4(:,:,2) = Psi_4(:,:,1);
Psi_4(:,:,1) = Psi_4temp;

%% update psi at corner node by solving 2L equations.

% M1, M2, M : linear system of 2L unkowns(in order to find PSI1, PSI2 at corner)
for i=1:L
    M1(i,i) = 1/(2*dx) + A(i)/(2*dt);
    M3(i,i) = 1/(2*dy) + A(i)/(2*dt);
    for j=1:L
        M2(i,j) = -C(i,j)/(2*dt);
    end
end
        
M=[M1 M2 ; M2 M3];


% Corner 1
Psi_corner_old(:,1,1) = Psi_2(:,2,1); % A at time n-step
Psi_corner_old(:,2,1) = Psi_1(:,2,1); % B at time n-step
Psi_corner_old(:,3,1) = Psi_1(:,1,1); % D at time n-step
Psi_corner_old(:,4,1) = Psi_2(:,1,1); % E at time n-step
 
Psi_corner_new(:,1,1) = Psi_2(:,2,2); % A at time n+1-step
Psi_corner_new(:,2,1) = Psi_1(:,2,2); % B at time n+1-step
 
% Corner 2
Psi_corner_old(:,1,2) = Psi_3(:,2,1); % A at time n-step
Psi_corner_old(:,2,2) = Psi_2(:,Nx-1,1); % B at time n-step
Psi_corner_old(:,3,2) = Psi_2(:,Nx,1); % D at time n-step
Psi_corner_old(:,4,2) = Psi_3(:,1,1); % E at time n-step

Psi_corner_new(:,1,2) = Psi_3(:,2,2); % A at time n+1-step
Psi_corner_new(:,2,2) = Psi_2(:,Nx-1,2); % B at time n+1-step

% Corner 3
Psi_corner_old(:,1,3) = Psi_4(:,Nx-1,1); % A at time n-step
Psi_corner_old(:,2,3) = Psi_3(:,Ny-1,1); % B at time n-step
Psi_corner_old(:,3,3) = Psi_3(:,Ny,1); % D at time n-step
Psi_corner_old(:,4,3) = Psi_4(:,Nx,1); % E at time n-step

Psi_corner_new(:,1,3) = Psi_4(:,Nx-1,2); % A at time n+1-step
Psi_corner_new(:,2,3) = Psi_3(:,Ny-1,2); % B at time n+1-step

% Corner 4
Psi_corner_old(:,1,4) = Psi_1(:,Ny-1,1); % A at time n-step
Psi_corner_old(:,2,4) = Psi_4(:,2,1); % B at time n-step
Psi_corner_old(:,3,4) = Psi_4(:,1,1); % D at time n-step
Psi_corner_old(:,4,4) = Psi_1(:,Ny,1); % E at time n-step

Psi_corner_new(:,1,4) = Psi_1(:,Ny-1,2); % A at time n+1-step
Psi_corner_new(:,2,4) = Psi_4(:,2,2); % B at time n+1-step

if n>2
    Psi_E = zeros(L,4);
    Psi_D = zeros(L,4);
    Psi_E= Psi_ED(1:L,:);
    Psi_D= Psi_ED(L+1:2*L,:);
    
    Psi_2(:,1,1) = Psi_E(:,1);
    Psi_1(:,1,1) = Psi_D(:,1);
    
    Psi_3(:,1,1) = Psi_E(:,2);
    Psi_2(:,Nx,1) = Psi_D(:,2);
    
    Psi_4(:,Nx,1) = Psi_E(:,3);
    Psi_3(:,Ny,1) = Psi_D(:,3);
    
    Psi_1(:,Ny,1) = Psi_E(:,4);
    Psi_4(:,1,1) = Psi_D(:,4);
end



for l=1:L
    temp1 = zeros(1,1,4);
    temp2 = zeros(1,1,4);
    temp3 = zeros(1,1,4);
    
    b(l,1,:) = (A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_old(l,4,:) +...
        (-A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_new(l,1,:)+...
        (A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_old(l,1,:) ;
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,2,:);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,3,:);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,2,:);
    end
    b(l,1,:)=b(l,1,:) + temp1 + temp2 + temp3 ;
end


for l=1:L
    temp1 = zeros(1,1,4);
    temp2 = zeros(1,1,4);
    temp3 = zeros(1,1,4);
    
    b(l+L,1,:) = (A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_old(l,3,:) +...
        (-A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_new(l,2,:)+...
        (A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_old(l,2,:) ;
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,1,:);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,4,:);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,1,:);
    end
    b(l+L,1,:)=b(l+L,1,:) + temp1 + temp2 + temp3 ;
end

Psi_ED = zeros(2*L,4);

for i=1:4
Psi_ED(:,i) = M\b(:,1,i) ; % value of corner nodes 'E and D' at time ndt.
end



end

% for n=1:250:Nt
%     
%     surf(X,Y,real(U(:,:,n)))
%     pause(0.5)
% end
% 
% for n=1:250:Nt
%     
%     surf(X,Y,real(v(X,Y,(n-1)*dt)))
%     pause(0.5)
% end
% 
% 
% L2error = zeros(Nt,2);
% for n=1:Nt
%   % L2error(:,:,n) = abs(U(:,:,n)-v(X,Y,(n-1)*dt));
%    L2error(n,1) = sqrt(norm(U(:,:,n)-v(X,Y,(n-1)*dt),'fro')*dx*dy);
%    
% end
% 
% for n =2:Nt-1
%   L2error(n+1,2) = L2error(n+1,1)/L2error(n,1);
% end
