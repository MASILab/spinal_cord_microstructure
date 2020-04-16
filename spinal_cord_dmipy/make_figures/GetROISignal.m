addpath(genpath('/Volumes/schillkg/MATLAB/NIFTI/'))

a = load_untouch_nii('dwmri.nii.gz');
dwi = a.img;

FC = load_untouch_nii('ROIs/FasciculusCuneatus.nii.gz')
FC = FC.img;

lCST = load_untouch_nii('ROIs/LateralCST.nii.gz')
lCST = lCST.img;

SL = load_untouch_nii('ROIs/SpinalLemniscus.nii.gz')
SL = SL.img;

vCST = load_untouch_nii('ROIs/VentralCST.nii.gz')
vCST = vCST.img;

VH = load_untouch_nii('ROIs/VentralHorn.nii.gz')
VH = VH.img;

figure; hold on;

name = 'FasciulusCuneatus'
var = FC;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite([name '.txt'],signal)
plot(signal);

name = 'LateralCST'
var = lCST;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite([name '.txt'],signal)
plot(signal);

name = 'SpinalLemniscus'
var = SL;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite([name '.txt'],signal)
plot(signal);

name = 'VentralCST'
var = vCST;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite([name '.txt'],signal)
plot(signal);

name = 'VentralHorn'
var = VH;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite([name '.txt'],signal)
plot(signal);

%% NORM
addpath(genpath('/Volumes/schillkg/MATLAB/NIFTI/'))
addpath(genpath('/Users/kurtschilling/MATLAB/NIFTI_20130306'))

a = load_untouch_nii_gz('dwmri_norm.nii.gz');
dwi = a.img;

FC = load_untouch_nii_gz('ROIs/FasciculusCuneatus.nii.gz')
FC = FC.img;

lCST = load_untouch_nii_gz('ROIs/LateralCST.nii.gz')
lCST = lCST.img;

SL = load_untouch_nii_gz('ROIs/SpinalLemniscus.nii.gz')
SL = SL.img;

vCST = load_untouch_nii_gz('ROIs/VentralCST.nii.gz')
vCST = vCST.img;

VH = load_untouch_nii_gz('ROIs/VentralHorn.nii.gz')
VH = VH.img;

figure; hold on;

name = 'FasciulusCuneatus'
var = FC;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite(['SingleVoxelSignals_norm' filesep name '.txt'],signal)
plot(signal);

name = 'LateralCST'
var = lCST;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite(['SingleVoxelSignals_norm' filesep name '.txt'],signal)
plot(signal);

name = 'SpinalLemniscus'
var = SL;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite(['SingleVoxelSignals_norm' filesep name '.txt'],signal)
plot(signal);

name = 'VentralCST'
var = vCST;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite(['SingleVoxelSignals_norm' filesep name '.txt'],signal)
plot(signal);

name = 'VentralHorn'
var = VH;
signal = [];
for i=1:936
    tmp = dwi(:,:,1,i);
    signal(i) = mean(tmp(var==1));
end
dlmwrite(['SingleVoxelSignals_norm' filesep name '.txt'],signal)
plot(signal);






























