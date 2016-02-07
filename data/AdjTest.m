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
loglog( dA, abs( dJdA - dJdAFD ), '.-k', dA, .1^1*dA, '-r' );
xlabel('$\Delta B$','Interpreter','Latex');
ylabel('|Adjoint - FD|');
% axis([1e-9 1e2 1e-12 1e-4]);
title('2 Particles, 8 time-steps Adjoint');
set(gca,'fontsize',25);