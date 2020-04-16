
clear; clc; close all;

addpath(genpath('/Volumes/schillkg/MATLAB/anderson'))

scheme = dlmread('DWMRI_all/SingleVoxelSignals_norm/scheme.scheme')

dirs = dlmread('DWMRI_all/SingleVoxelSignals_norm/BVECS.bvec')

signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/VentralHorn.txt')

x = scheme(:,5);
b0 = scheme(:,4);
y = signal(b0==0);
TE = x(b0==0);

figure; plot(TE,y,'o');

ms=10;
figure; hold on;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   bval(i) = median(scheme(indices,4));
   indices_minus_b0 = indices; indices_minus_b0(1)=[];
   sig = signal(indices_minus_b0);
   bvals = median(scheme(indices,4)).*ones(length(sig),1);
   SIG(i) = mean(sig);
   plot(bvals,sig,'.','MarkerSize',ms);
end

 figure; plot(bval,SIG,'o');
 
 
   