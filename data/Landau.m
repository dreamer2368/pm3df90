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
E = fread(fileID,Nx*Ny*Nz*3*Nt,'double');
E = reshape(E,[Nx,Ny,Nz,3,Nt]);

fileID = fopen('PE.bin');
PE = fread(fileID,Nt,'double');

fileID = fopen('KE.bin');
KE = fread(fileID,Nt,'double');

fileID = fopen('phi.bin');
phi = fread(fileID,Nx*Ny*Nz,'double');
phi = reshape(phi,[Nx,Ny,Nz]);

%%
close all

x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny); z = linspace(0,Lz,Nz); [X,Y,Z]=meshgrid(x,y,z);

s=5;
t = 32;

%video clip
writerObj = VideoWriter('close_encounter.avi');
writerObj.FrameRate = 40;
open(writerObj);

for i=1:Nt
    figure(1)
%     plot(xp(:,1,i),vp(:,1,i),'.k');
    scatter3(squeeze(xp(:,1,i)),squeeze(xp(:,2,i)),squeeze(xp(:,3,i)),50);
    axis([0 Lx 0 Ly 0 Lz]);
    title('2 Particles with $v_0=1.0$','Interpreter','Latex');
    set(gca,'fontsize',25);
%     axis([0 Lx -0.6 0.6]);
%     title('Spatial distribution');
%     xlabel('$x$','Interpreter','Latex');
%     ylabel('$v$','Interpreter','Latex');

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

    %videoclip
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
%     pause(.01);
end

% videoclip close
close(writerObj);

%%
close all

r1 = squeeze( abs(xp(1,:,:) - xp(2,:,:)) );
r1 = sqrt( sum( r1.^2 , 1 ) );
time = 0.01*(1:Nt)/2/pi;

figure(2)
plot(time,r1);

r = squeeze( min( abs(xp(1,:,:) - xp(2,:,:)), Lx - abs(xp(1,:,:) - xp(2,:,:)) ) );
r = sqrt( sum( r.^2 , 1 ) );
figure(3)
plot(time, r);
title(strcat('$\vert x_1 - x_2 \vert$ with $v_0 = 1.0$ : $\overline{r}=$',num2str(mean(r(500:Nt)))),'Interpreter','Latex');
xlabel('$T/T_p$','Interpreter','Latex');
ylabel('minimum distance');
set(gca,'fontsize',25);

%%
close all

time = 0.1*(1:Nt);
wp = 1; eps0 = 1;
tp = 2*pi/wp;
PE = eps0*PE;

figure(4)
plot(time/tp,KE,'-r',time/tp,PE,'-b',time/tp,KE+PE,'-k');
title('Two particles with $v_0=1.0$','Interpreter','Latex');
xlabel('$T/T_p$','Interpreter','Latex');
ylabel('Energy');
legend('Kinetic','Potential','Total');
set(gca,'fontsize',25);

% dEp = squeeze( sum(E,1) )/max(max(max(E)));
% figure(5)
% plot(time,dEp(3,:));

%%
close all

fileID = fopen('ek.bin');
% ek = fread(fileID,6*20,'double');
% ek = reshape(ek,[20,6]);
ek = fread(fileID,2*20,'double');
ek = reshape(ek,[20,2]);
time = [1/2/pi, 1, 2, 4, 10, 20]; time = time*2*pi;

fileID = fopen('dA.bin');
dA = fread(fileID,20,'double');

fileID = fopen('dJdA.bin');
ek0 = fread(fileID,6,'double');

figure(5)
loglog(dA,.1^3*dA,'-r');
hold on
for i=1:2
    loglog(dA,ek(:,i),'.-');
end
title('Landau-damping with $N=2^{6+6+6}$','Interpreter','Latex');
xlabel('$\Delta$','Interpreter','Latex');
ylabel('Normalized error');
h=legend('$\mathcal{O}(\Delta^{-1})$','$T_p/4\pi$','$T_p/2\pi$','$T_p$','$2T_p$','$4T_p$','$10T_p$','$20T_p$');
set(h,'Interpreter','Latex');
set(gca,'fontsize',25);

figure(6)
semilogy(time,ek0,'o-');
xlabel('$T$','Interpreter','Latex');
ylabel('Sensitivity');
title('Sensitivity divergence with $v_0=0.2$','Interpreter','Latex');
set(gca,'fontsize',25);

%%
close all
time = 0.1*(1:Nt);

x = Lx/Nx*(1:Nx); y = Ly/Ny*(1:Ny); z = Lz/Nz*(1:Nz);
[Y, Z] = meshgrid(y,z);

%This shows that efield is somewhat insensitive to y and z direction,
%and k=0 spectrum is not necessarily zero : average over y and z is an
%effective value.
figure(7)
mesh(Y,Z,squeeze(E(8,:,:,1,1)));

dy = Ly/Ny; dz = Lz/Nz;
Exb = dy*dz/Ly/Lz*squeeze( sum( sum(E(:,:,:,1,:),2), 3) );
NFFT = 2^nextpow2(Nx); histSpectrum = [];
for i=1:Nt
    Ekb = fft(Exb(:,i),NFFT)/Nx;
    histSpectrum = [histSpectrum 2*abs(Ekb(1:NFFT/2+1))];
end
figure(8)
semilogy(time,histSpectrum(2,:),'-k',time,0.1*exp(-0.1539*time),'-r');
title('Landau damping in 3D CIC');
xlabel('time');
ylabel('$\mathcal{F}\{\overline{E_x}\}$','Interpreter','Latex');
legend('3D CIC','Analytic');
set(gca,'fontsize',25);

%%
close all
fileID = fopen('Fpx.bin');
Fpx = fread(fileID,1000,'double');
fileID = fopen('xd.bin');
xd = fread(fileID,1000,'double');

figure(8)
plot(xd-Lx/2,Fpx,'.k',xd-Lx/2,1./(xd-Lx/2)./(xd-Lx/2),'-r');
axis([-Lx/2 Lx/2 min(Fpx) max(Fpx)]);