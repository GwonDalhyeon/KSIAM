% function Collinocorner(L)
clc
clear all
%% Constants and Variables
% sdksajfdkl;safjdklsa;fjsad;lfjsdl;afjkl;asda;f


asdj;flkj;asdklf;
L=5;
T=0.1;
dx = 1/20;
dy = 1/20;
dt = 1/1000;
Nx = 21;
Ny = 21;
Nt = round(T/dt)
alpha = zeros(L,1);
beta = zeros(L,1);
C = zeros(L);
A = zeros(L,1);
u = zeros(Nx+1,Ny+1,Nt+1);5
angle = zeros(L,1);
M =zeros(2*L,2*L); % matrix corresponding to linear system for obtaining auxiliary function on the semi-boundary node
P1 = zeros(L) ; 
P2 = zeros(L) ;
psi_old = zeros(L,4);
psi_new = zeros(L,4);
b = zeros(2*L,1);
angle = [pi/(2*L+1): pi/(2*L+1) : L*pi/(2*L+1)];
alpha = (cos(angle)).^2;
beta =  2/(2*L+1) * (sin(angle)).^2;

%% Fourier transforming with initial and boundary conditions, find auxiliary function in a domain.


syms x y t k1 k2 w
f = @(x,y,t) exp(-x^2-y^2-t^2);
g = fourier(fourier(f,x,k1),t,w);
psi_cell = cell(4,L);
for l=1:5
    
%     eval(['phi1_temp', num2str(l), ' = @(k1,w)',char(g)]);
%     psi_cell{1,l} = eval(['phi2_temp', num2str(l), '=@(k2,w)', char(g)]);
    eval(['phi2_temp', num2str(l), '=@(k2,w)', char(g)]);
    psi_cell{1,l} = eval(['phi2_temp', num2str(l)]);
    
end
    

%% initial data pluginÇØ¾ßÇÔ

for l=1:L
    for m=1:L
          C(l,m) = beta(m)/(alpha(l)+alpha(m)-alpha(l)*alpha(m));
    end
    A(l) = 1+ sum(beta./alpha)+ alpha(l)*C(l,m);
end


%% solving 2L linear system of auxiliary variable psi on the semi-boundary node.

for i=1:L
    P1(i,i) = 1/(2*dx) + A(i)/(2*dt);
    for j=1:L
        P2(i,j)=-C(i,j)/(2*dt);
    end
end
        
M=[P1 P2 ; P2 P1];

temp1 = 0;
temp2 = 0;
temp3 = 0;

for l=1:L
    b(l) = (A(l)/(2*dt) - 1/(2*dx)) * psi_old(l,4) +...
        (-A(l)/(2*dt) + 1/(2*dx)) * psi_new(l,1)+...
        (A(l)/(2*dt) + 1/(2*dx)) * psi_old(l,1) ;
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * psi_new(m,2);
        temp2 = temp2 - C(l,m)/(2*dt) * psi_old(m,3);
        temp3 = temp3 - C(l,m)/(2*dt) * psi_old(m,2);
    end
    b(l)=b(l) + temp1 + temp2 + temp3 ;
end

for l=1:L
    b(l+L) = (A(l)/(2*dt) - 1/(2*dy)) * psi_old(l,3) +...
        (-A(l)/(2*dt) + 1/(2*dy)) * psi_new(l,2)+...
        (A(l)/(2*dt) + 1/(2*dy)) * psi_old(l,2) ;
    for m=1:L
        temp1 = temp1 + C(l,m)/(2*dt) * psi_new(m,1);
        temp2 = temp2 - C(l,m)/(2*dt) * psi_old(m,4);
        temp3 = temp3 - C(l,m)/(2*dt) * psi_old(m,1);
    end
    b(l)=b(l) + temp1 + temp2 + temp3 ;
end

x = M\b ; % value of corner nodes at time ndt.

%% Get the value of auxiliary variable phi on the regular grid

phi_1 = zeros(L,Ny,3); % phi_1(l.ny,1) : 'n-1' time, phi_1(l.ny,2) : 'n' time, phi_1(l.ny,3) : 'n+1' time,    
phi_2 = zeros(L,Nx,3);
psi_1 = zeros(L,Ny);
psi_2 = zeros(L,Nx);

for l=1:L
    for i=2:Nx
        psi_2(l,i) = alpha(l)*phi_2(l,i,2) + 1/2 * (u(i,1,n)+u(i,2,n)); % n-loop
        
    end
    
    for j=2:Ny
        psi_1(l,j) = alpha(l)*phi_1(l,j,2) + 1/2 * (u(1,j,n)+u(2,j,n));
    end
end

for l=1:L
    for i=2:Nx-1
        phi_2(l,i,3) = (dt)^2/(dx)^2 * (psi_2(l,i-1) - 2*psi_2(l,i) + psi_2(l,i+1)) + 2* phi_2(l,i,2) - phi_2(l,i,1);        
    end
    
    for j=2:Ny-1
        phi_1(l,j,3) = (dt)^2/(dy)^2 * (psi_1(l,j-1) - 2*psi_1(l,j) + psi_1(l,j+1)) + 2* phi_1(l,j,2) - phi_1(l,j,1);
    end
end

% for l=1:L
%    i=Nx;
%         phi_2(l,i,3) = (dt)^2/(dx)^2 * (psi_2(l,i-1) - 2*psi_2(l,i) + psi_2(l,i+1)) + 2* phi_2(l,i,2) - phi_2(l,i,1);        
%         
%     j=Ny;
%         phi_1(l,j,3) = (dt)^2/(dy)^2 * (psi_1(l,j-1) - 2*psi_1(l,j) + psi_1(l,j+1)) + 2* phi_1(l,j,2) - phi_1(l,j,1);
% end


%% Get the solution at the semi-boundary node.
tempsum = zeros(Nx,1);
for l=1:L
    tempsum= tempsum + beta(l)* (phi_2(l,:,n)-phi_2(l,:,n-1));
end

u(:,1,n) = (dy)/(dt+dy) * (-u(:,2,n)+u(:,2,n-1)+u(:,1,n-1)) + (dt)/(dt+dy) *(-u(:,1,n-1) + u(:,1,n)+u(:,1,n-1)) +...
    (2*dy)/(dt+dy) * tempsum;


%% Get the solution on entire domain.

for i=2:Nx
    for j=2:Nx
        u(i,j,n+1) = (dt)^2/(dx)^2 * (u(i+1,j,n)-2*u(i,j,n) + u(i-1,j,n)) + ...
            (dt)^2/(dy)^2 * (u(i,j+1,n)-2*u(i,j,n) + u(i,j-1,n)) + ...
            2*u(i,j,n) - u(i,j,n-1);
    end
end




