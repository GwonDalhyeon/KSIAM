% function Collinocorner(L)
clc
clear all
%% Constants and Variables
L=5;
Tlength=1;
Xlength = 1;
Ylength = 1;
dx = 1/50;
dy = 1/50;
dt = 1/5000;
Nx = round(Xlength/dx)+2 ;
Ny = round(Ylength/dy)+2 ;
Nt = round(Tlength/dt)+1;
alpha = zeros(L,1);
beta = zeros(L,1);
C = zeros(L);
A = zeros(L,1);
U = zeros(Ny,Nx,Nt); % solution matrix
angle = zeros(L,1);
% matrix corresponding to linear system for obtaining auxiliary function on the semi-boundary node
MM1 =zeros(2*L,2*L); %corner1
MM2 =zeros(2*L,2*L); %corner2
MM3 =zeros(2*L,2*L); %corner3
MM4 =zeros(2*L,2*L); %corner4
Mxp = zeros(L) ; 
Mxm = zeros(L) ;
Myp = zeros(L) ;
Mym = zeros(L) ;
M = zeros(L);
b = zeros(2*L,1,4);
angle = [pi/(2*L+1): pi/(2*L+1) : L*pi/(2*L+1)];
alpha = (cos(angle)).^2;
beta =  2/(2*L+1) * (sin(angle)).^2;

for l=1:L
    for m=1:L
          C(l,m) = beta(m)/(alpha(l)+alpha(m)-alpha(l)*alpha(m));
    end
    A(l) = 1+ sum(beta./alpha)- alpha(l)*sum(C(l,:)./transpose(alpha(:)));
end

% model problem 
syms x y t
% u_tt = u_xx + u_yy 
% u(x,y,0) = @(x,y) exp(i*(pi*x+pi*y)) , u_t(x,y,0) = -sqrt(2)*pi*i*exp(i*(pi*x+pi*y)) : initial condition

u = @(x,y) exp(1i*(pi*x+pi*y));
v = @(x,y,t) exp(1i*(pi*x+pi*y-sqrt(2)*pi*t));
u_t = @(x,y) -sqrt(2)*pi*1i*exp(1i*(pi*x+pi*y));

% u = @(x,y) exp(1i*(pi*x+0*y));
% u_t = @(x,y) -pi*1i * exp(1i*(pi*x+0*y));
% v = @(x,y,t) exp(1i*(pi*x+0*y-pi*t));

% u = @(x,y) exp(1i*(-pi*y+0*x));
% u_t = @(x,y) -pi*1i * exp(1i*(-pi*y+0*x));
% v = @(x,y,t) exp(1i*(-pi*y+0*x+pi*t));

%  u = @(x,y) exp(i*(pi*x+pi*y));
%  v = @(x,y,t) exp(sqrt(-1)*(pi*x+pi*y-sqrt(2)*pi*t));
%  u_t = @(x,y) -sqrt(2)*pi*i*exp(i*(pi*x+pi*y));

[X,Y] = meshgrid(-dx/2:dx:Xlength+dx/2,-dy/2:dy:Ylength+dy/2); 
% [Y,X] = meshgrid(-dy/2:dy:Ylength+dy/2,-dx/2:dx:Xlength+dx/2); 
%% get initial data for U0, PHI0 

U(:,:,1) = u(X,Y) - dt/2* u_t(X,Y);
U(:,:,2) = u(X,Y) + dt/2* u_t(X,Y); 

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
    Psi_1(l,:,1) = alpha(l)*Phi_1(l,:,1) + 1/2 * transpose(U(:,1,1)+U(:,1,1));
    Psi_2(l,:,1) = alpha(l)*Phi_2(l,:,1) + 1/2 * (U(1,:,1)+U(2,:,1));
    Psi_3(l,:,1) = alpha(l)*Phi_3(l,:,1) + 1/2 * transpose(U(:,Nx-1,1)+U(:,Nx,1));
    Psi_4(l,:,1) = alpha(l)*Phi_4(l,:,1) + 1/2 * (U(Ny-1,:,1)+U(Ny,:,1));
    
    % Psi for time n=2 data
    Psi_1(l,:,2) = alpha(l)*Phi_1(l,:,2) + 1/2 * transpose(U(:,1,2)+U(:,2,2));
    Psi_2(l,:,2) = alpha(l)*Phi_2(l,:,2) + 1/2 * (U(1,:,2)+U(2,:,2));
    Psi_3(l,:,2) = alpha(l)*Phi_3(l,:,2) + 1/2 * transpose(U(:,Nx-1,2)+U(:,Nx,2));
    Psi_4(l,:,2) = alpha(l)*Phi_4(l,:,2) + 1/2 * (U(Ny-1,:,2)+U(Ny,:,2));
end


%% Main iteration

for n=2:Nt-1
    qq=1;  

%% updating U at interior nodes
for i=2:Nx-1
    for j=2:Ny-1
U(j,i,n+1) = 2*U(j,i,n) - U(j,i,n-1) +(dt^2)/(dx^2)*(U(j,i+1,n)-2*U(j,i,n) +U(j,i-1,n)) +...
    +(dt^2)/(dy^2) * ( U(j+1,i,n)-2*U(j,i,n) +U(j-1,i,n)); % at Nx-grid, the value of u is not changing. 
    end
end


%% update PHI 


for i= 2: Nx-1

    Phi_2(:,i,1) = 2* Phi_2(:,i,2) - Phi_2(:,i,1) + (dt)^2/(dx)^2 * (Psi_2(:,i-1,2) - 2*Psi_2(:,i,2) + Psi_2(:,i+1,2)) ;
    Phi_4(:,i,1) = 2* Phi_4(:,i,2) - Phi_4(:,i,1) + (dt)^2/(dx)^2 * (Psi_4(:,i-1,2) - 2*Psi_4(:,i,2) + Psi_4(:,i+1,2)) ;

end

for j= 2:Ny-1
    Phi_1(:,j,1) = 2* Phi_1(:,j,2) - Phi_1(:,j,1) + (dt)^2/(dy)^2 * (Psi_1(:,j-1,2) - 2*Psi_1(:,j,2) + Psi_1(:,j+1,2)) ;
    Phi_3(:,j,1) = 2* Phi_3(:,j,2) - Phi_3(:,j,1) + (dt)^2/(dy)^2 * (Psi_3(:,j-1,2) - 2*Psi_3(:,j,2) + Psi_3(:,j+1,2)) ;
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

%% update U at boundary nodes

tempsum1 = zeros(1,Ny-2);
tempsum2 = zeros(1,Nx-2);
tempsum3 = zeros(1,Ny-2);
tempsum4 = zeros(1,Nx-2);

for l=1:L
    tempsum1 = tempsum1 + beta(l)* 1/dt* (Phi_1(l,2:Ny-1,2)-Phi_1(l,2:Ny-1,1));
    tempsum2 = tempsum2 + beta(l)* 1/dt* (Phi_2(l,2:Nx-1,2)-Phi_2(l,2:Nx-1,1));
    tempsum3 = tempsum3 + beta(l)* 1/dt* (Phi_3(l,2:Ny-1,2)-Phi_3(l,2:Ny-1,1));
    tempsum4 = tempsum4 + beta(l)* 1/dt* (Phi_4(l,2:Nx-1,2)-Phi_4(l,2:Nx-1,1));
end

U(2:Ny-1,1,n+1) = -2/(1/dt-1/dx) *((+1/(2*dx) + 1/(2*dt))*U(2:Ny-1,2,n+1) + (-1/(2*dt) + 1/(2*dx))*U(2:Ny-1,2,n) + (-1/(2*dt) - 1/(2*dx))*U(2:Ny-1,1,n)) +...
    +2/(1/dt-1/dx) *transpose(tempsum1);

U(1,2:Nx-1,n+1) = -2/(1/dt-1/dy) *((1/(2*dy) + 1/(2*dt))*U(2,2:Nx-1,n+1) + (-1/(2*dt) + 1/(2*dy))*U(2,2:Nx-1,n) + (-1/(2*dt) - 1/(2*dy))*U(1,2:Nx-1,n)) +...
    +2/(1/dt-1/dy) *tempsum2;

U(2:Ny-1,Nx,n+1) = -2/(1/dt+1/dx) *((-1/(2*dx) + 1/(2*dt))*U(2:Ny-1,Nx-1,n+1) + (-1/(2*dt) - 1/(2*dx))*U(2:Ny-1,Nx-1,n) + (-1/(2*dt) + 1/(2*dx))*U(2:Ny-1,Nx,n)) +...
    +2/(1/dt+1/dx) *transpose(tempsum3);

U(Ny,2:Nx-1,n+1) = -2/(1/dt+1/dy) *((-1/(2*dy) + 1/(2*dt))*U(Ny-1,2:Nx-1,n+1) + (-1/(2*dt) - 1/(2*dy))*U(Ny-1,2:Nx-1,n) + (-1/(2*dt) + 1/(2*dy))*U(Ny,2:Nx-1,n)) +...
    +2/(1/dt+1/dy) *tempsum4;

%% update PSI except at corner nodes

for l=1:L
        
    Psi_1(l,2:Ny-1,1) = alpha(l)*Phi_1(l,2:Ny-1,2) + 1/2 * transpose(U(2:Ny-1,1,2)+U(2:Ny-1,2,2));
    Psi_2(l,2:Nx-1,1) = alpha(l)*Phi_2(l,2:Nx-1,2) + 1/2 * (U(1,2:Nx-1,2)+U(2,2:Nx-1,2));
    Psi_3(l,2:Ny-1,1) = alpha(l)*Phi_3(l,2:Ny-1,2) + 1/2 * transpose(U(2:Ny-1,Nx,2)+U(2:Ny-1,Nx-1,2));
    Psi_4(l,2:Nx-1,1) = alpha(l)*Phi_4(l,2:Nx-1,2) + 1/2 * (U(Ny,2:Nx-1,2)+U(Ny-1,2:Nx-1,2));
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

%% update PSI at corner node by solving 2L equations.

% M1, M2, M : linear system of 2L unkowns(in order to find PSI1, PSI2 at corner)
for i=1:L
    Mxp(i,i) = 1/(2*dx) + A(i)/(2*dt);
    Mxm(i,i) = -1/(2*dx) + A(i)/(2*dt);
    Myp(i,i) = 1/(2*dy) + A(i)/(2*dt);
    Mym(i,i) = -1/(2*dy) + A(i)/(2*dt);
    for j=1:L
        M(i,j) = -C(i,j)/(2*dt);
    end
end
        
MM1 = [Mxm M ; M Mym];
MM2 = [Mxp M ; M Mym];
MM3 = [Mxp M ; M Myp];
MM4 = [Mxm M ; M Myp];

Psi_corner_old = zeros(L,4,4);
Psi_corner_new = zeros(L,4,4);


% Corner 1
Psi_corner_old(:,4,1) = Psi_2(:,2,1); % A at time n-step
Psi_corner_old(:,3,1) = Psi_1(:,2,1); % B at time n-step
Psi_corner_old(:,2,1) = Psi_1(:,1,1); % D at time n-step
Psi_corner_old(:,1,1) = Psi_2(:,1,1); % E at time n-step

Psi_corner_new(:,4,1) = Psi_2(:,2,2); % A at time n+1-step
Psi_corner_new(:,3,1) = Psi_1(:,2,2); % B at time n+1-step
Psi_corner_new(:,2,1) = Psi_1(:,1,2); % D at time n+1-step
Psi_corner_new(:,1,1) = Psi_2(:,1,2); % E at time n+1-step
 
% Corner 2
Psi_corner_old(:,1,2) = Psi_3(:,2,1); % A at time n-step
Psi_corner_old(:,3,2) = Psi_2(:,Nx-1,1); % B at time n-step
Psi_corner_old(:,2,2) = Psi_2(:,Nx,1); % D at time n-step
Psi_corner_old(:,4,2) = Psi_3(:,1,1); % E at time n-step

Psi_corner_new(:,1,2) = Psi_3(:,2,2); % A at time n+1-step
Psi_corner_new(:,3,2) = Psi_2(:,Nx-1,2); % B at time n+1-step
Psi_corner_new(:,2,2) = Psi_2(:,Nx,2); % D at time n+1-step
Psi_corner_new(:,4,2) = Psi_3(:,1,2); % E at time n+1-step

% Corner 3
Psi_corner_old(:,1,3) = Psi_4(:,Nx-1,1); % A at time n-step
Psi_corner_old(:,2,3) = Psi_3(:,Ny-1,1); % B at time n-step
Psi_corner_old(:,3,3) = Psi_3(:,Ny,1); % D at time n-step
Psi_corner_old(:,4,3) = Psi_4(:,Nx,1); % E at time n-step

Psi_corner_new(:,1,3) = Psi_4(:,Nx-1,2); % A at time n+1-step
Psi_corner_new(:,2,3) = Psi_3(:,Ny-1,2); % B at time n+1-step
Psi_corner_new(:,3,3) = Psi_3(:,Ny,2); % D at time n+1-step
Psi_corner_new(:,4,3) = Psi_4(:,Nx,2); % E at time n+1-step

% Corner 4
Psi_corner_old(:,4,4) = Psi_1(:,Ny-1,1); % A at time n-step
Psi_corner_old(:,2,4) = Psi_4(:,2,1); % B at time n-step
Psi_corner_old(:,3,4) = Psi_4(:,1,1); % D at time n-step
Psi_corner_old(:,1,4) = Psi_1(:,Ny,1); % E at time n-step

Psi_corner_new(:,4,4) = Psi_1(:,Ny-1,2); % A at time n-step
Psi_corner_new(:,2,4) = Psi_4(:,2,2); % B at time n-step
Psi_corner_new(:,3,4) = Psi_4(:,1,2); % D at time n-step
Psi_corner_new(:,1,4) = Psi_1(:,Ny,2); % E at time n-step




% plugging 1~Lth elements of b.
for l=1:L
    
        
    for i=1:4
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    if i==1
   
    b(l,1,1) = (A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_old(l,1,1) +...
        (-A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_new(l,4,1)+...
        (A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_old(l,4,1) ;
    
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,3,1);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,2,1);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,3,1);
    end
    b(l,1,1)=b(l,1,1) + temp1 + temp2 + temp3 ;
    
    end
    
    if i==2
        
    b(l,1,2) = (A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_old(l,1,2) +...
        (-A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_new(l,1,2)+...
        (A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_old(l,4,2) ;
    
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,3,2);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,2,2);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,3,2);
    end
    b(l,1,2)=b(l,1,2) + temp1 + temp2 + temp3 ;
        
    end
    
    if i==3
        
        b(l,1,3) = (A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_old(l,1,3) +...
        (-A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_new(l,1,3)+...
        (A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_old(l,4,3) ;
    
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,2,3);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,3,3);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,2,3);
    end
    b(l,1,3)=b(l,1,3) + temp1 + temp2 + temp3 ;
    end
    
    if i==4
        
        b(l,1,4) = (A(l)/(2*dt) + 1/(2*dx)) * Psi_corner_old(l,1,4) +...
        (-A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_new(l,4,4)+...
        (A(l)/(2*dt) - 1/(2*dx)) * Psi_corner_old(l,4,4) ;
    
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,2,4);
        temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,3,4);
        temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,2,4);
    end
    b(l,1,4)=b(l,1,4) + temp1 + temp2 + temp3 ;
    end
    
    end
end

% plugging L+1~2Lth elements of b. 
for l=1:L
    
    for i=1:4
        temp1 = 0; 
        temp2 = 0; 
        temp3 = 0; 
        % i denotes the location of corner.
        
    if i==1
        
    b(l+L,1,1) = (A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_old(l,2,1) +...
        (-A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_new(l,3,1)+...
        (A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_old(l,3,1) ;
        for m=1:L
            temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,4,1);
            temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,1,1);
            temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,4,1);
        end
    b(l+L,1,1)=b(l+L,1,1) + temp1 + temp2 + temp3 ;
    end
    
    if i==2
        
    b(l+L,1,2) = (A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_old(l,2,2) +...
        (-A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_new(l,3,2)+...
        (A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_old(l,3,2) ;
        for m=1:L
            temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,1,2);
            temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,4,2);
            temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,1,2);
        end
    b(l+L,1,2)=b(l+L,1,2) + temp1 + temp2 + temp3 ;
    end
    
     if i==3
        
    b(l+L,1,3) = (A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_old(l,3,3) +...
        (-A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_new(l,2,3)+...
        (A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_old(l,2,3) ;
        for m=1:L
            temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,1,3);
            temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,4,3);
            temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,1,3);
        end
    b(l+L,1,3)=b(l+L,1,3) + temp1 + temp2 + temp3 ;
     end
    
      if i==4
        
    b(l+L,1,4) = (A(l)/(2*dt) - 1/(2*dy)) * Psi_corner_old(l,3,4) +...
        (-A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_new(l,2,4)+...
        (A(l)/(2*dt) + 1/(2*dy)) * Psi_corner_old(l,2,4) ;
        for m=1:L
            temp1 = temp1 + C(l,m)/(2*dt) * Psi_corner_new(m,4,4);
            temp2 = temp2 - C(l,m)/(2*dt) * Psi_corner_old(m,1,4);
            temp3 = temp3 - C(l,m)/(2*dt) * Psi_corner_old(m,4,4);
        end
    b(l+L,1,4)=b(l+L,1,4) + temp1 + temp2 + temp3 ;
      end
    end
end

Psi_ED1 = zeros(2*L,1);
Psi_ED2 = zeros(2*L,1);
Psi_ED3 = zeros(2*L,1);
Psi_ED4 = zeros(2*L,1);

Psi_ED1 = MM1\b(:,1,1) ; % value of corner nodes 'E and D' at time ndt.
Psi_ED2 = MM2\b(:,1,2) ;
Psi_ED3 = MM3\b(:,1,3) ;
Psi_ED4 = MM4\b(:,1,4) ;
    
    Psi_2(:,1,2) = Psi_ED1(1:L);
    Psi_1(:,1,2) = Psi_ED1(L+1:2*L);
    
    Psi_2(:,Nx,2) = Psi_ED2(1:L);
    Psi_3(:,1,2) = Psi_ED2(L+1:2*L);
   
    Psi_4(:,Nx,2) = Psi_ED3(1:L);
    Psi_3(:,Ny,2) = Psi_ED3(L+1:2*L);
    
    
    Psi_4(:,1,2) = Psi_ED4(1:L);
    Psi_1(:,Ny,2) = Psi_ED4(L+1:2*L);

    Un = U(:,:,n);
    Vn = abs(U(:,:,n)-U(:,:,n-1));
    Wn = max(max(Vn));
    Qn = v(X,Y,(n-3/2)*dt);
    Rn = abs(U(:,:,n)-Qn);
end
% figure
% axis([-0.2 1.2 -0.2 1.2 -0.1 2]) 
% for n=1:20:Nt
%     
%      surf(X(10:Ny-1,10:Nx-1),Y(10:Ny-1,10:Nx-1),real(U(10:Ny-1,10:Nx-1,n)))
%     
%     
%     pause(0.5)
% end
for n=1:20:Nt
    
     surf(X(2:Ny-1,2:Nx-1),Y(2:Ny-1,2:Nx-1),real(U(2:Ny-1,2:Nx-1,n)))
    
    
    pause(0.5)
end

% for n=1:20:Nt
%     
%     surf(X,Y,real(v(X,Y,(n-3/2)*dt)))
%     pause(0.5)
% end

% 
% L2error = zeros(Nt,2);
% for n=1:Nt
%   % L2error(:,:,n) = abs(U(:,:,n)-v(X,Y,(n-1)*dt));
%    L2error(n,1) = sqrt(norm(U(:,:,n)-v(X,Y,(n-3/2)*dt),'fro')*dx*dy);
%    
% end
% 
% for n =2:Nt-1
%   L2error(n+1,2) = L2error(n+1,1)/L2error(n,1);
% end