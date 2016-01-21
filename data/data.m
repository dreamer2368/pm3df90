clear all
close all
clc

spec = importdata('record.out');
N = spec(1); Nx = spec(2); Ny = spec(3); Nz = spec(4); Nt = spec(5);
L = importdata('L.out');

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

x = linspace(0,L,Nx); y = linspace(0,L,Ny); z = linspace(0,L,Nz); [X,Y,Z]=meshgrid(x,y,z);

s=5;
t = 32;
for i=1:Nt
%     figure(1)
%     plot(xp(:,1,i),vp(:,1,i),'.k');
%     axis([0 2.1 -0.6 0.6]);
    figure(1)
    scatter3(squeeze(xp(:,1,i)),squeeze(xp(:,2,i)),squeeze(xp(:,3,i)),s);
    axis([0 L 0 L 0 L]);
    title('Spatial distribution');
    xlabel('$x$','Interpreter','Latex');
    ylabel('$y$','Interpreter','Latex');
    zlabel('$z$','Interpreter','Latex');
    set(gca,'fontsize',25);
%     view([0 1 0]);
    
    figure(2)
    mesh(squeeze(Y(:,:,t)),squeeze(X(:,:,t)),squeeze(Ex(:,:,t,i)));
%     view([1 0 0]);
    axis([0 L 0 L -1.1 1.1]);
    xlabel('$x$','Interpreter','Latex'); ylabel('$y$','Interpreter','Latex'); zlabel('$E(x,y)$','Interpreter','Latex');
    title('$E_x(x,y,z=\frac{L}{2})$','Interpreter','Latex');
    set(gca,'fontsize',25);

    figure(3)
    plot(reshape(Y,[Nx*Ny*Nz,1]),reshape(Ex(:,:,:,1),[Nx*Ny*Nz,1]),'*k');
    axis([0 L -1.1 1.1]);
    title('$E_x(x)$ for one sheet','Interpreter','Latex');
    xlabel('$x$','Interpreter','Latex');
    ylabel('$E_x$','Interpreter','Latex');
    set(gca,'fontsize',25);
    pause();
end
