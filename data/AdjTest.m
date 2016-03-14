close all
clear all
clc

%%
close all
clc

fileID = fopen('dJdA.bin');
dJdA = fread(fileID,1,'double');

fileID = fopen('dA.bin');
dA = fread(fileID,20,'double');

fileID = fopen('dJdAFD.bin');
dJdAFD = fread(fileID,20,'double');

error = abs( dJdA - dJdAFD )/abs(dJdA);

figure(1)
loglog( dA, error, '.-k', dA, .1^1*(dA), '-r' );
hold on
loglog( dA(5), error(5),'or', 'LineWidth',2,'markers', 24);
loglog( dA(6), error(6),'or', 'LineWidth',2,'markers', 24);
xlabel('$\Delta B$','Interpreter','Latex');
ylabel('Normalized error');
% axis([1e-9 1e2 1e-12 1e-4]);
title('Two-stream, $N=2^8\times2^5\times2^5$, $T\simeq0.25T_p$','Interpreter','Latex');
set(gca,'fontsize',25);