clear all
close all
clc

fDA = importdata('fDA.out');
% ek = importdata('ek_twostream.out');
ek = importdata('ek.out');
ek = reshape(ek,[length(ek)/length(fDA), length(fDA)]);
Time = importdata('Tf.out');
dJdA = importdata('dJdA.out');
J0 = importdata('J0.out');

% figure(1)
% loglog(fDA,.0001*fDA,'-r');
% hold on
% loglog(fDA,ek0,'.-k');
% for i=1:length(Time);
%     loglog(fDA,ek(i,:),'.-');
% end
% title('Dicrete Adjoint degradation - Double,CIC','fontsize',20);
% xlabel('\Delta B','fontsize',20);
% ylabel('error','fontsize',20);
% h=legend('O(\Delta B)','langmuir10T_p','T_p/2\pi','0.5T_p','T_p','2T_p','4T_p','6T_p','8T_o','10T_p');
% set(h,'fontsize',12);

Fig = figure(1);
set(Fig,'Position',[100,100,650,550]);
loglog(fDA,.01*fDA,'-r','LineWidth',5);
hold on
% loglog(fDA,ek0,'.-k');
loglog(fDA,ek(1,:),'+-','LineWidth',5);
loglog(fDA,ek(length(Time),:),'o-','LineWidth',5);
for i=2:length(Time)-1;
    loglog(fDA,ek(i,:),'--');
end
% axis([1e-10 1e5 1e-15 1e5]);
axis([1e-10 1e1 1e-15 1e5]);
title('CIC, $N=2$, $v_0=0.6$','Interpreter','Latex');
xlabel('$\Delta$','Interpreter','Latex');
ylabel('error','Interpreter','Latex');
h=legend('$\mathcal{O}(\Delta)$','$T_p/2\pi$','$20T_p$');
% h=legend('$\mathcal{O}(\Delta B)$','$T_p/2\pi$','$10T_p$','$0.5T_p$','$T_p$','$2T_p$','$4T_p$','$6T_p$','$8T_p$');
set(h,'fontsize',25,'Interpreter','Latex');
set(gca,'fontsize',35);

ytick = [4e-10  9e-10];
yticklabel = {'-10^-^3','0','10^-^3'};

Fig = figure(2);
set(Fig,'Position',[100,100,650,550]);
semilogy(Time,abs(dJdA),'o-','LineWidth',5);
xlabel('$T_f-T_I$','Interpreter','Latex');
ylabel('$\frac{\partial J}{\partial B}$','Interpreter','Latex');
title('$Sensitivity$','Interpreter','Latex');
set(gca,'fontsize',35);
% set(gca,'ytick',ytick,'yticklabel',yticklabel);
% set(gca,'position',[0.2 0.2 0.7 0.6]);