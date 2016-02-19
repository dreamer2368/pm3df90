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
writerObj = VideoWriter('Repulsive2.avi');
writerObj.FrameRate = 25;
open(writerObj);

for i=1:Nt
    figure(1)
%     plot(xp(:,1,i),vp(:,1,i),'.k');
    scatter3(squeeze(xp(:,1,i)),squeeze(xp(:,2,i)),squeeze(xp(:,3,i)),50);
    axis([0 Lx 0 Ly 0 Lz]);
    title('2 Particles with $v_0=0.9$','Interpreter','Latex');
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
time = 0.05*(1:Nt)/2/pi;

figure(2)
plot(time,r1);

r = squeeze( min( abs(xp(1,:,:) - xp(2,:,:)), Lx - abs(xp(1,:,:) - xp(2,:,:)) ) );
r = sqrt( sum( r.^2 , 1 ) );
figure(3)
plot(time, r);
title(strcat('$\vert x_1 - x_2 \vert$ with $v_0 = 0.9$ : $\overline{r}=$',num2str(mean(r(500:Nt)))),'Interpreter','Latex');
xlabel('$T/T_p$','Interpreter','Latex');
ylabel('minimum distance');
set(gca,'fontsize',25);

%%
close all

time = 0.05*(1:Nt);
wp = 1; eps0 = 1;
tp = 2*pi/wp;
PE = eps0*PE;

figure(4)
plot(time/tp,KE,'-r',time/tp,PE,'-b',time/tp,KE+PE,'-k');
title('Two particles with $v_0=0.9$','Interpreter','Latex');
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
ek = fread(fileID,6*20,'double');
ek = reshape(ek,[20,6]);
ek0 = [8.4764477945905060E-003, 2.2032771176303872E-003, 4.3627461166659406E-003, 3.1766359150327572, 38.088234823333075, 44491365.287056476];
time = [1/2/pi, 1, 2, 4, 10, 20]; time = time*2*pi;

fileID = fopen('dA.bin');
dA = fread(fileID,20,'double');

figure(5)
for i=1:6
    loglog(dA,ek(:,i)/ek0(i),'.-');
    hold on
end
title('Degradation in 2 particles with $v_0=0.2$','Interpreter','Latex');
xlabel('$\Delta$','Interpreter','Latex');
ylabel('Normalized error');
h=legend('$T_p/2\pi$','$T_p$','$2T_p$','$4T_p$','$10T_p$','$20T_p$');
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
time = 0.2*(1:Nt);

figure(7)
for i=1:2
    semilogy(time, squeeze(abs(E(i,2,:))) );
    hold on
    semilogy(time, squeeze(abs(E(i,3,:))) );
end
xlabel('time');
ylabel('$\vert E_{y,z}\vert$','Interpreter','Latex');
title('E-field acting on each particle in y,z');
set(gca,'fontsize',25);