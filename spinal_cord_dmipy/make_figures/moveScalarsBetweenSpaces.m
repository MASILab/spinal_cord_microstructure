% eddy b0 is standard space
% want to move C3_02 to standard space

clc; clear; close all;
addpath(genpath('/Users/kurtschilling/MATLAB'))

% for each, load eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs, dwmri_all.bval
nii = load_untouch_nii_gz('eddy_b0.nii.gz')
ref=nii.img; 

nii2 = load_untouch_nii('C3_02.nii')
moving = nii2.img;

[selectedMovingPoints,selectedFixedPoints] = cpselect(mat2gray(moving),mat2gray(ref),'Wait',true);
t_concord = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'nonreflectivesimilarity');
disp(['scale: ' num2str(t_concord.T(1,1))])
disp(['rot: ' num2str(asind(t_concord.T(2,1)))])
       
ortho_ref = imref2d(size(ref)); %relate intrinsic and world coordinates
   
registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
    
% mask registered and ref with cylinder
% save

mask = load_untouch_nii_gz('CylinderMask.nii.gz')
mask = mask.img;

A = ref.*single(mask);

nii.img = A;
save_untouch_nii_gz(nii,'eddy_b0_masked.nii.gz')

B = registered.*int16(mask);
nii.img = B;
save_untouch_nii_gz(nii,'C3_02_2b0_masked.nii.gz')


% also crop and mask the b0 and FA from dti_10 I guess
nii2 = load_untouch_nii_gz('diff_10/dti__FA.nii.gz')
C = nii2.img;
C = squeeze(C(:,:,4));
C = C.*single(mask);
nii.img = C;
save_untouch_nii_gz(nii,'dti_10_FA_slice4_masked.nii.gz')

nii2 = load_untouch_nii_gz('diff_10/eddy_unwarped_images.nii.gz')
D = nii2.img;
D = squeeze(D(:,:,4,37));
D = D.*single(mask);
nii.img = D;
save_untouch_nii_gz(nii,'dti_10_eddy_slice4_masked.nii.gz')



%% 
% load dti_10_FA and a single examp;le of all ROIs
% then load moving T1 again
% register and save

clc; clear; close all;
addpath(genpath('/Users/kurtschilling/MATLAB'))

% for each, load eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs, dwmri_all.bval
nii = load_untouch_nii_gz('eddy_b0.nii.gz')
ref=nii.img; 

ROI = load_untouch_nii_gz('ROIs/FasciculusCuneatus.nii.gz')
ROI = ROI.img;
ref(ROI==1)=0;
ROI = load_untouch_nii_gz('ROIs/LateralCST.nii.gz')
ROI = ROI.img;
ref(ROI==1)=0;
ROI = load_untouch_nii_gz('ROIs/SpinalLemniscus.nii.gz')
ROI = ROI.img;
ref(ROI==1)=0;
ROI = load_untouch_nii_gz('ROIs/VentralCST.nii.gz')
ROI = ROI.img;
ref(ROI==1)=0;
ROI = load_untouch_nii_gz('ROIs/VentralHorn.nii.gz')
ROI = ROI.img;
ref(ROI==1)=0;


nii2 = load_untouch_nii('C3_02.nii')
moving = nii2.img;

[selectedMovingPoints,selectedFixedPoints] = cpselect(mat2gray(moving),mat2gray(ref),'Wait',true);
t_concord = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'nonreflectivesimilarity');
disp(['scale: ' num2str(t_concord.T(1,1))])
disp(['rot: ' num2str(asind(t_concord.T(2,1)))])
       
ortho_ref = imref2d(size(ref)); %relate intrinsic and world coordinates
   
registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
    
% mask registered and ref with cylinder
% save

mask = load_untouch_nii_gz('CylinderMask.nii.gz')
mask = mask.img;

A = ref.*single(mask);

nii.img = A;
save_untouch_nii_gz(nii,'eddy_b0_masked.nii.gz')

B = registered.*int16(mask);
nii.img = B;
save_untouch_nii_gz(nii,'C3_02_2b0_masked.nii.gz')








