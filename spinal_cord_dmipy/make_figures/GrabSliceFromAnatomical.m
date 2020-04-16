addpath(genpath('/Users/kurtschilling/MATLAB'))


nii = load_untouch_nii('mFFE_01.nii')
img = nii.img;

size(img)

figure; imagesc(squeeze(img(:,:,6)));

c = squeeze(img(:,:,6));

nii.img = c;
nii.hdr.dime.dim(4) = 1;
save_untouch_nii(nii,'C3_01.nii')


nii = load_untouch_nii('mFFE_02.nii')
img = nii.img;

size(img)

figure; imagesc(squeeze(img(:,:,6)));

c = squeeze(img(:,:,6));

nii.img = c;
nii.hdr.dime.dim(4) = 1;
save_untouch_nii(nii,'C3_02.nii')

% use 02




