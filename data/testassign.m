clear all
close all
clc

N=2; Nx=64; Ny=64; Nz=64;
L = importdata('L.out');

fileID = fopen('testphi.bin');
testphi = fread(fileID,Nx*Ny*Nz,'double');
testphi = reshape(testphi,[Nx,Ny,Nz]);

fileID = fopen('testE.bin');
testE = fread(fileID,Nx*Ny*Nz*3,'double');
testE = reshape(testE,[Nx, Ny, Nz,3]);

%%
close all

x = linspace(0,L,Nx); y = linspace(0,L,Ny); z = linspace(0,L,Nz);
[X,Y,Z] = meshgrid(x,y,z);

s = 32;
mesh(X(:,:,s),Y(:,:,s),testphi(:,:,s));