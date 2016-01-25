clear all
close all
clc

spec = importdata('record');
N = spec(1); Nx = spec(2); Ny = spec(3); Nz = spec(4); Nt = spec(5);
Lx = spec(6); Ly = spec(7); Lz = spec(8);

fileID = fopen('xp.bin');
xp = fread(fileID,N*3*Nt,'double');
xp = reshape(xp,[N,3,Nt]);

fileID = fopen('vp.bin');
vp = fread(fileID,N*3*Nt,'double');
vp = reshape(vp,[N,3,Nt]);

fileID = fopen('Ex.bin');
Ex = fread(fileID,Nx*Ny*Nz*Nt,'double');
Ex = reshape(Ex,[Nx,Ny,Nz,Nt]);

%%
close all

x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny); z = linspace(0,Lz,Nz); [X,Y,Z]=meshgrid(x,y,z);

s=5;
t = 32;
for i=1:Nt
    figure(1)
    plot(xp(:,1,i),vp(:,1,i),'.k');
%     scatter3(squeeze(xp(:,1,i)),squeeze(xp(:,2,i)),squeeze(xp(:,3,i)),s);
    axis([0 Lx -0.6 0.6]);
    title('Spatial distribution');
    xlabel('$x$','Interpreter','Latex');
    ylabel('$v$','Interpreter','Latex');
    set(gca,'fontsize',25);
%     
%     figure(2)
%     mesh(squeeze(Y(:,:,t)),squeeze(X(:,:,t)),squeeze(Ex(:,:,t,i)));
% %     view([1 0 0]);
%     axis([0 Ly 0 Lx -1.1 1.1]);
%     xlabel('$x$','Interpreter','Latex'); ylabel('$y$','Interpreter','Latex'); zlabel('$E(x,y)$','Interpreter','Latex');
%     title('$E_x(x,y,z=\frac{L}{2})$','Interpreter','Latex');
%     set(gca,'fontsize',25);

%     figure(3)
%     plot(reshape(Y,[Nx*Ny*Nz,1]),reshape(Ex(:,:,:,1),[Nx*Ny*Nz,1]),'*k');
%     axis([0 Lx -1.1 1.1]);
%     title('$E_x(x)$ for one sheet','Interpreter','Latex');
%     xlabel('$x$','Interpreter','Latex');
%     ylabel('$E_x$','Interpreter','Latex');
%     set(gca,'fontsize',25);
    pause(.001);
end

%%
close all

t=16;
dx = Lx/Nx; dt = 0.2;
NFFT = 2^nextpow2(Nx);
Ekhistory = [];
for i=1:Nt
    Ek = fft(Ex(:,t,t,i),NFFT)/Nx;
    Ekhistory = [Ekhistory 2*abs(Ek(1:NFFT/2+1))];
end

time = (1:Nt)*dt;
semilogy(time, Ekhistory(2,:),'-k',time,0.00001*exp(1/2/sqrt(2)*time),'-r');