clear all
close all
clc

spec0 = importdata('record0');
N0 = spec0(1); Nx0 = spec0(2); Ny0 = spec0(3); Nz0 = spec0(4); Nt0 = spec0(5);
Lx0 = spec0(6); Ly0 = spec0(7); Lz0 = spec0(8);

fileID = fopen('xp0.bin');
xp0 = fread(fileID,N0*3*Nt0,'double');
xp0 = reshape(xp0,[N0,3,Nt0]);

fileID = fopen('vp0.bin');
vp0 = fread(fileID,N0*3*Nt0,'double');
vp0 = reshape(vp0,[N0,3,Nt0]);

fileID = fopen('E0.bin');
E0 = fread(fileID,Nx0*Ny0*Nz0*3*Nt0,'double');
E0 = reshape(E0,[Nx0,Ny0,Nz0,3,Nt0]);

fileID = fopen('rho0.bin');
rho0 = fread(fileID,Nx0*Ny0*Nz0*Nt0,'double');
rho0 = reshape(rho0,[Nx0,Ny0,Nz0,Nt0]);

fileID = fopen('PE0.bin');
PE0 = fread(fileID,Nt0,'double');

fileID = fopen('KE0.bin');
KE0 = fread(fileID,Nt0,'double');

spec1 = importdata('record1');
N1 = spec1(1); Nx1 = spec1(2); Ny1 = spec1(3); Nz1 = spec1(4); Nt1 = spec1(5);
Lx1 = spec1(6); Ly1 = spec1(7); Lz1 = spec1(8);

fileID = fopen('xp1.bin');
xp1 = fread(fileID,N1*3*Nt1,'double');
xp1 = reshape(xp1,[N1,3,Nt1]);

fileID = fopen('vp1.bin');
vp1 = fread(fileID,N1*3*Nt1,'double');
vp1 = reshape(vp1,[N1,3,Nt1]);

fileID = fopen('E1.bin');
E1 = fread(fileID,Nx1*Ny1*Nz1*3*Nt1,'double');
E1 = reshape(E1,[Nx1,Ny1,Nz1,3,Nt1]);

fileID = fopen('rho1.bin');
rho1 = fread(fileID,Nx0*Ny0*Nz0*Nt0,'double');
rho1 = reshape(rho1,[Nx0,Ny0,Nz0,Nt0]);

fileID = fopen('PE1.bin');
PE1 = fread(fileID,Nt1,'double');

fileID = fopen('KE1.bin');
KE1 = fread(fileID,Nt1,'double');

spec2 = importdata('record2');
N2 = spec2(1); Nx2 = spec2(2); Ny2 = spec2(3); Nz2 = spec2(4); Nt2 = spec2(5);
Lx2 = spec2(6); Ly2 = spec2(7); Lz2 = spec2(8);

fileID = fopen('xp2.bin');
xp2 = fread(fileID,N2*3*Nt2,'double');
xp2 = reshape(xp2,[N2,3,Nt2]);

fileID = fopen('vp2.bin');
vp2 = fread(fileID,N2*3*Nt2,'double');
vp2 = reshape(vp2,[N2,3,Nt2]);

fileID = fopen('E2.bin');
E2 = fread(fileID,Nx2*Ny2*Nz2*3*Nt2,'double');
E2 = reshape(E2,[Nx2,Ny2,Nz2,3,Nt2]);

fileID = fopen('rho2.bin');
rho2 = fread(fileID,Nx0*Ny0*Nz0*Nt0,'double');
rho2 = reshape(rho2,[Nx0,Ny0,Nz0,Nt0]);

fileID = fopen('PE2.bin');
PE2 = fread(fileID,Nt2,'double');

fileID = fopen('KE2.bin');
KE2 = fread(fileID,Nt2,'double');

%% Configuration-space comparison
close all

%video clip
writerObj = VideoWriter('testcase.avi');
writerObj.FrameRate = 1;
open(writerObj);

for i=1:Nt1
    figure(1)
    plot(squeeze(xp0(:,1,i)),squeeze(vp0(:,1,i)),'.k');
    axis([0 Lx0 -.3 .3]);
    xlabel('x'); ylabel('v');
    title('Two-stream instability for $0.25T_p$','Interpreter','latex');
    set(gca,'fontsize',25);
    
%     figure(1)
%     scatter3(squeeze(xp1(:,1,i)),squeeze(xp1(:,2,i)),squeeze(xp1(:,3,i)),1);
%     axis([0 Lx1 0 Ly1 0 Lz1]);
%     title(strcat('Two-stream with $B=',num2str(exp(-5)),'$'),'Interpreter','Latex');
%     view([0 1 0]);
%     set(gca,'fontsize',25);
%     
%     figure(2)
%     scatter3(squeeze(xp2(:,1,i)),squeeze(xp2(:,2,i)),squeeze(xp2(:,3,i)),1);
%     axis([0 Lx2 0 Ly2 0 Lz2]);
%     title(strcat('Two-stream with $B=',num2str(exp(-6)),'$'),'Interpreter','Latex');
%     view([0 1 0]);
%     set(gca,'fontsize',25);
%     
%     figure(3)
%     scatter3(squeeze(vp1(:,1,i)),squeeze(vp1(:,2,i)),squeeze(vp1(:,3,i)),5);
%     axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
%     title(strcat('Two-stream with $B=',num2str(exp(-5)),'$'),'Interpreter','Latex');
%     view([0 1 0]);
%     set(gca,'fontsize',25);
%     
%     figure(4)
%     scatter3(squeeze(vp2(:,1,i)),squeeze(vp2(:,2,i)),squeeze(vp2(:,3,i)),5);
%     axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
%     title(strcat('Two-stream with $B=',num2str(exp(-6)),'$'),'Interpreter','Latex');
%     view([0 1 0]);
%     set(gca,'fontsize',25);

    %videoclip
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
%     pause();
end

% videoclip close
close(writerObj);

%% Efield Comparison in x direction
close all

dx1 = Lx1/Nx1; dy1 = Ly1/Ny1; dz1 = Lz1/Nz1; dx2 = Lx2/Nx2; dy2 = Ly2/Ny2; dz2 = Lz2/Nz2;
x1 = dx1*(1:Nx1); x2 = dx2*(1:Nx2);

Ebx0 = dy1*dz1/Ly1/Lz1*squeeze( sum( sum(E1(:,:,:,1,:),2), 3) );
Ebx1 = dy2*dz2/Ly2/Lz2*squeeze( sum( sum(E2(:,:,:,1,:),2), 3) );

for i=1:Nt1
    figure(1)
    plot(x1,Ebx0(:,i),'-k',x2,Ebx1(:,i),'-r');
    axis([0 Lx1 -1e-5 1e-5]);
    set(gca,'fontsize',25);
    
    pause();
end

%% Velocity perturbation comparison
close all

vb0 = abs(vp1); vb1 = abs(vp2);
vb0(:,1,:) = vb0(:,1,:)-0.2; vb1(:,1,:) = vb1(:,1,:)-0.2;

for i=1:Nt1 
    figure(3)
    plot(xp1(:,1,i),vb0(:,1,i),'ok',xp2(:,1,i),vb1(:,1,i),'.r');
%     axis([0 Lx0 -1e-5 1e-5]);
    title('velocity perturbation','Interpreter','Latex');
    set(gca,'fontsize',25);
    
    pause();
end

%% 1D1V phase-space comparison
close all

for i=1:Nt1
    figure(1)
    plot(xp1(:,1,i),vp1(:,1,i),'ok',xp2(:,1,i),vp2(:,1,i),'.r');
    axis([0 Lx1 -.3 .3]);
    set(gca,'fontsize',25);
    
    pause(1.0);
end

%% Deviation comparison
close all
dx1 = Lx1/Nx1; dy1 = Ly1/Ny1; dz1 = Lz1/Nz1; dx2 = Lx2/Nx2; dy2 = Ly2/Ny2; dz2 = Lz2/Nz2;
x1 = dx1*(1:Nx1); x2 = dx2*(1:Nx2);

dxp1 = xp1-xp0; dvp1 = vp1-vp0; dE1 = E1-E0; drho1 = rho1-rho0;
dxp2 = xp2-xp0; dvp2 = vp2-vp0; dE2 = E2-E0; drho2 = rho2-rho0;

dEbx1 = dy2*dz2/Ly1/Lz1*squeeze( sum( sum(dE1(:,:,:,1,:),2), 3) );
dEbx2 = dy2*dz2/Ly2/Lz2*squeeze( sum( sum(dE2(:,:,:,1,:),2), 3) );

drhobx1 = dy2*dz2/Ly1/Lz1*squeeze( sum( sum(drho1,2), 3) );
drhobx2 = dy2*dz2/Ly1/Lz1*squeeze( sum( sum(drho2,2), 3) );

for i=1:Nt0
    figure(1)
    plot(dxp1(:,1,i),dvp1(:,1,i),'.k',dxp2(:,1,i),dvp2(:,1,i),'.r');
    title('Phase-space deviation');
    xlabel('$\delta x_p$','Interpreter','latex');
    ylabel('$\delta v_p$','Interpreter','latex');
    h=legend('$\delta B=e^{-5}$','$\delta B=e^{-6}$');
    set(h,'Interpreter','latex');
    set(gca,'fontsize',25);
    axis([-3e-3 3e-3 -1.5e-3 1.5e-3]);
    
%     figure(3)
% %     hold off
% %     plot(xp1(N1/2,1,i),dxp1(N1/2,1,i),'ok',xp1(N1/2+1,1,i),dxp1(N1/2+1,1,i),'*k');
% %     plot(xp2(N2/2,1,i),dxp2(N2/2,1,i),'or',xp2(N2/2+1,1,i),dxp2(N2/2+1,1,i),'*k');
% %     hold on
%     plot(xp1(:,1,i),dxp1(:,1,i),'.k',xp2(:,1,i),dxp2(:,1,i),'.r');
%     axis([0 Lx0 -1e-3 1e-3]);
% %     axis([1-2e-1 1+2e-1 -1e-3 1e-3]);
%     
%     figure(4)
%     plot(xp1(:,1,i),vp1(:,1,i),'.k',xp2(:,1,i),vp2(:,1,i),'.r');
%     axis([0 Lx0 -.3 .3]);

    figure(2)
    plot(x1,dEbx1(:,i),'-k',x2,dEbx2(:,i),'-r');
    xlabel('x'); ylabel('$\delta \overline{E}_x$','Interpreter','latex');
    title('Deviation in E-field');
    h=legend('$\delta B=e^{-5}$','$\delta B=e^{-6}$');
    set(h,'Interpreter','latex');
    set(gca,'fontsize',25);
    
    figure(5)
    plot(x1,drhobx1(:,i),'-k',x2,drhobx2(:,i),'-r');
    xlabel('x'); ylabel('$\delta \overline{\rho}$','Interpreter','latex');
    title('Deviation in $\rho(x)$','Interpreter','latex');
    h=legend('$\delta B=e^{-5}$','$\delta B=e^{-6}$');
    set(h,'Interpreter','latex');
    set(gca,'fontsize',25);
    
    pause();
end

%% Particles that have a position deviation as large as Lx
close all
clc

dxp1 = xp1-xp0; dvp1 = vp1-vp0; dE1 = E1-E0;
dxp2 = xp2-xp0; dvp2 = vp2-vp0; dE2 = E2-E0;

%This number difference shows that particles of the position deviation=Lx
%are not the bounday particles.
nnz( (abs(dxp1)>Lx0/2) )/2^10
nnz( (abs(dxp2)>Lx0/2) )/2^10

for i=1:Nt0
    outliers1 = (abs(dxp1(:,1,i))>Lx0/2);
%     nnz(outliers1)
%     if( nnz(outliers1)==0 )
%         outliers1=ones(N0,1);
%     end
    sxp1 = squeeze(xp1(:,1,i)); svp1 = squeeze(vp1(:,1,i));
    figure(1)
    plot(sxp1(outliers1),svp1(outliers1),'.k');
    axis([0 Lx1 -.3 .3]);
    
    pause();
end

%% J comparison in time : MKE
close all
dx1 = Lx1/Nx1; dy1 = Ly1/Ny1; dz1 = Lz1/Nz1; dx2 = Lx2/Nx2; dy2 = Ly2/Ny2; dz2 = Lz2/Nz2;
x1 = dx1*(1:Nx1); x2 = dx2*(1:Nx2);
time = 0.2*(1:Nt0);

J0 = 1/N0/(Nt0-2)*sum( vp0(:,:,3:Nt0).^2 ); J0 = squeeze(J0);
J1 = 1/N1/(Nt1-2)*sum( vp1(:,:,3:Nt1).^2 ); J1 = squeeze(J1);
J2 = 1/N2/(Nt2-2)*sum( vp2(:,:,3:Nt2).^2 ); J2 = squeeze(J2);
dJ1 = J1-J0; dJ2 = J2-J0;

figure(1)
for i=1:3
    plot(time(3:Nt0), squeeze(dJ1(i,:)),'-k',time(3:Nt0),squeeze(dJ2(i,:)),'-r');
    hold on
end

J0x = 1/N0/(Nt0-2)*sum( vp0(:,1,3:Nt0).^2 ); J0x = squeeze(J0x);
J1x = 1/N1/(Nt1-2)*sum( vp1(:,1,3:Nt1).^2 ); J1x = squeeze(J1x);
J2x = 1/N2/(Nt2-2)*sum( vp2(:,1,3:Nt2).^2 ); J2x = squeeze(J2x);
dJ1x = J1x-J0x; dJ2x = J2x-J0x;

figure(2)
plot(time(3:Nt0),dJ1x,'-k',time(3:Nt0),dJ2x,'-r');
%% J comparison in time : MPE
close all
dx1 = Lx1/Nx1; dy1 = Ly1/Ny1; dz1 = Lz1/Nz1; dx2 = Lx2/Nx2; dy2 = Ly2/Ny2; dz2 = Lz2/Nz2;
x1 = dx1*(1:Nx1); x2 = dx2*(1:Nx2);
time = 0.2*(1:Nt0);

J0 = 1/Nx0/Ny0/Nz0/(Nt0-2)*sum(sum(sum(sum( E0(:,:,:,:,3:Nt0).^2 )))); J0 = squeeze(J0);
J1 = 1/Nx1/Ny1/Nz1/(Nt1-2)*sum(sum(sum(sum( E1(:,:,:,:,3:Nt1).^2 )))); J1 = squeeze(J1);
J2 = 1/Nx2/Ny2/Nz2/(Nt2-2)*sum(sum(sum(sum( E2(:,:,:,:,3:Nt2).^2 )))); J2 = squeeze(J2);
dJ1 = J1-J0; dJ2 = J2-J0;

figure(2)
plot(time(3:Nt0), J1-J0,'-k',time(3:Nt0),J2-J0,'-r');

%Here, the scale of MPE is 10^(-9~-12), which is significantly smaller than
%that of MKE : 10^(-2 ~ -3). Although both J changes with the scale of dJ :
%10^-15, because MKE is more large scale quantity, only MKE is susceptible
%to round-off error.

%% Assignment error?
close all

for i=1:Nt0
    figure(1)
    histogram(squeeze(xp1(:,1,i)));

    figure(2)
    histogram(squeeze(xp2(:,1,i)));
    
    pause()
end