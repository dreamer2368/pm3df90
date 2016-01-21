clear all
close all
clc

B = linspace(0.9,1.1,401);
Jdata = importdata('Jdata.out');

%%
close all
clc

Fig = figure(1);
set(Fig,'Position',[100,100,750,500]);
i = size(Jdata,2)-1;
h=plot(B,(Jdata(:,i)-Jdata(201,i))/Jdata(201,i),'--b','LineWidth',1 );
hold on
i = size(Jdata,2);
h=plot(B,(Jdata(:,i)-Jdata(201,i))/Jdata(201,i),'-r','LineWidth',3 );
for i = 1:size(Jdata,2)-2
    h = plot(B,(Jdata(:,i)-Jdata(201,i))/Jdata(201,i),'--k','LineWidth',1 );
end
% axis([0.9 1.1 -1.5e-3 1e-3]);
axis([0.9 1.1 -0.05 0.05]);
h=legend('$\le 10T_p$','$20T_p$');
set(h,'Interpreter','Latex');
xlabel('$B$','Interpreter','Latex');
ylabel('$\frac{\delta\mathcal{J}}{\mathcal{J}_0}$','Interpreter','Latex','Units','Normalized', 'Position', [-0.2, 0.5, 0]);
title('$\mathcal{J}$ fluctuation','Interpreter','Latex');
set(gca,'fontsize',40);

ytick = [-0.05 0 0.05];
yticklabel = {'$-0.05$','$\mathcal{J}_0$','$0.05$'};
ax = gca;
hold(ax,'on');
set(ax,'TickLabelInterpreter','Latex');
set(ax,'ytick',ytick,'yticklabel',yticklabel);
set(ax,'position',[0.3 0.2 0.6 0.65]);