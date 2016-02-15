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

fileID = fopen('E.bin');
E = fread(fileID,N*3*Nt,'double');
E = reshape(E,[N,3,Nt]);

%%
close all

x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny); z = linspace(0,Lz,Nz); [X,Y,Z]=meshgrid(x,y,z);

s=5;
t = 32;
for i=1:Nt
    figure(1)
%     plot(xp(:,1,i),vp(:,1,i),'.k');
    scatter3(squeeze(xp(:,1,i)),squeeze(xp(:,2,i)),squeeze(xp(:,3,i)),50);
    axis([0 Lx 0 Ly 0 Lz]);
%     axis([0 Lx -0.6 0.6]);
%     title('Spatial distribution');
%     xlabel('$x$','Interpreter','Latex');
%     ylabel('$v$','Interpreter','Latex');
%     set(gca,'fontsize',25);
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
    pause(.01);
end

%%
close all

r1 = squeeze( abs(xp(1,:,:) - xp(2,:,:)) );
r1 = sqrt( sum( r1.^2 , 1 ) );
time = 0.2*(1:Nt);

figure(2)
plot(time,r1);

r = squeeze( min( abs(xp(1,:,:) - xp(2,:,:)), Lx - abs(xp(1,:,:) - xp(2,:,:)) ) );
r = sqrt( sum( r.^2 , 1 ) );
figure(3)
plot(time, r);

%%
close all

time = 0.2*(1:Nt);
KE = 0.5*squeeze( sum(sum( vp.^2, 1),2) );
PE = 0.5*squeeze( sum(sum( E.^2, 1),2) );

figure(4)
plot(time,KE,'-r',time,PE,'-b',time,KE+PE,'-k');

dEp = squeeze( sum(E,1) )/max(max(max(E)));
figure(5)
plot(time,dEp(3,:));