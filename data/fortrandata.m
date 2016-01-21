clear all
close all
clc

xp = importdata('xp.out');  xp = xp';
vp = importdata('vp.out');  vp = vp';
xp0 = importdata('xp0.out');  xp0 = xp0';
vp0 = importdata('vp0.out');  vp0 = vp0';
% E = importdata('E.out');    E = E';
% fDA = importdata('fDA.out');
% ek = importdata('ek.out');

%%
close all

ytick = [-1e-3 0 1e-3];
yticklabel = {'-10^-^3','0','10^-^3'};

Fig = figure(1);
set(Fig,'Position',[100,100,650,550]);
plot(xp(:,285),vp(:,285),'.k');
axis([0 2*pi/3.0619 -1 1]);
xlabel('$X$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
h=title('Particle distribution\\','Interpreter','Latex');
% set(h,'Position',[pi/3.0619,1.3e-3]);
set(gca,'fontsize',45);
ax = gca;
% hold(ax,'on');
% keyboard
% set(ax,'ytick',ytick,'yticklabel',yticklabel);
% set(ax,'position',[0.2 0.2 0.6 0.6]);

%video clip
writerObj = VideoWriter('N2chaos.avi');
writerObj.FrameRate = 25;
open(writerObj);
figure(1)
for i = 1:size(xp,2)
    plot(xp(:,i),vp(:,i),'or',xp0(:,i),vp(:,i),'*k','markers',10);
    axis([0 2*pi/3.0619 -1 1]);
    xlabel('$X$','Interpreter','Latex');
    ylabel('$V$','Interpreter','Latex');
    h=title('Particle distribution','Interpreter','Latex');
%     set(h,'Position',[pi/3.0619,1.3e-3]);
    set(gca,'fontsize',45);
%     ytick = [-1e-3 0 1e-3];
%     yticklabel = {'-10^-^3','0','10^-^3'};
%     set(gca,'ytick',ytick,'yticklabel',yticklabel);
%     set(gca,'position',[0.2 0.2 0.7 0.6]);
    h=legend('$B+\delta B$','$B$','Location','east');
    set(h,'Interpreter','Latex','fontsize',30);    

    %videoclip
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
%     pause(.00001);
end
% videoclip close
close(writerObj);

% figure(1)
% plot(xp(:,1),vp(:,1),'.k');
% axis([0 2*pi/3.0619 -1e-3 1e-3]);
% xlabel('X');
% ylabel('V');
% title('Particle distribution');
% set(gca,'fontsize',25);

% Ng = 64;
% x = linspace(0,2*pi,Ng+1)'; x = x(1:Ng) + 2*pi/Ng*0.5;
% figure(2)
% for i = 1:size(xp,2)
%     plot(x,E(:,i));
%     axis([0 2*pi -0.1 0.1]);
%     pause(.1);
% end
% 
% figure(3)
% loglog(fDA,ek,'.-k',fDA,fDA,'-r');
% title('discrete adjoint error','fontsize',20)
% xlabel('\Delta A','fontsize',20)
% ylabel('error','fontsize',20)
% legend('double','O(\Delta A)')