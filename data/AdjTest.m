close all
clear all
clc

%%
fileID = fopen('dJdA.bin');
dJdA = fread(fileID,1,'double');

fileID = fopen('dA.bin');
dA = fread(fileID,20,'double');

fileID = fopen('dJdAFD.bin');
dJdAFD = fread(fileID,20,'double');

figure(1)
loglog( dA(3:20), abs( dJdA - dJdAFD(3:20) ), '.-k', dA(3:20), .1^5*dA(3:20), '-r' );
xlabel('$\Delta x_{2,y}$','Interpreter','Latex');
ylabel('|Adjoint - FD|');
% axis([1e-9 1e2 1e-12 1e-4]);
title('single time-step adjoint');
set(gca,'fontsize',25);