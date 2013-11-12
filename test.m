clc
clear all
close all

L=5;
T=0.1;
dx = 1/20;
dy = 1/20;
dt = 1/1000;
Nx = 21;
Ny = 21;
Nt = round(T/dt);
alpha = zeros(L,1);
beta = zeros(L,1);
C = zeros(L);
A = zeros(L,1);
u = zeros(Nx+1,Ny+1,Nt+1);
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


f = @(x,y,t) exp(-x.^2-y.^2-t.^2);

[X,Y] = meshgrid(-1:dx:1,-1:dy:1); 
Xsize=size(X,2);
Ysize=size(Y,2);
tGrid = [-0.5:0.1:0.5];
tGridSize=size(tGrid,2);
B = zeros(tGridSize,Xsize);
%% Find PHI2
for i=1:tGridSize
    B(i,:) = f(X(1,:),0,tGrid(i));
end

FB = fft2(B);

[K1,W] = meshgrid(-1:dx:1,-0.5:0.1:0.5);
[K2,W] = meshgrid(-1:dx:1,-0.5:0.1:0.5);
% W = [-0.5:0.1:0.5];
FB_temp = zeros(tGridSize, Xsize, L);

for l=1:L
    for i=1:tGridSize
        FB_temp(i,:,l) = K1(1,:).^2/(W(i,:).^2 - alpha(l).*K1(1,:).^2).*FB(i,:);
    end
end

IFB2 = ifft2(FB_temp); % inverse fourier transform
%% Find PHI1

C= f(0,Y,T);
FC = fft2(C);

for i=1:tGridSize
    C(i,:) = f(0,Y(:,1)',tGrid(i));
end

FC = fft2(C);

FC_temp = zeros(tGridSize, Ysize, L);

for l=1:L
    for i=1:tGridSize
        FC_temp(i,:,l) = K2(1,:).^2/(W(i,:).^2 - alpha(l).*K2(1,:).^2).*FC(i,:);
    end
end

IFC2 = ifft(FC_temp); % inverse fourier transform