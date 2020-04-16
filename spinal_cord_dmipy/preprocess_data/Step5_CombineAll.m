% register all to #17
% load #17 and b0 from this (slice 4), get 4 points
% also get mask from manual made from 17
% for each, load eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs, dwmri_all.bval
% grab b0 from this (slice 4)
% get 4 points
% ICP to #17, apply to b0 and all dwi from forward PE (minus b=5)
% cat to image, cat bvals, cat bvecs and save as all
% in same loop, normalize to b0. b0 will have value of 1 everywhere. Mask with mask
% save

clc; clear all; close all;
addpath(genpath('/Volumes/schillkg/MATLAB/NIFTI/'))
%addpath(genpath('/Users/schilling/Downloads/ImageRegistrationApp/'))

%% load #17 and b0 from this (slice 4), get 4 points
nii = load_untouch_nii('diff_10/eddy_unwarped_images.nii.gz')
ref=nii.img;
b0=squeeze(ref(:,:,4,37));

nii.img = b0;
nii.hdr.dime.dim(4)=1;
nii.hdr.dime.dim(5)=1;
save_untouch_nii(nii,'eddy_b0.nii.gz')

CordMask_nii = load_untouch_nii('CordMask.nii.gz')
CordMask=CordMask_nii.img;
figure; imagesc(CordMask)

CyMask_nii = load_untouch_nii('CylinderMask.nii.gz')
CyMask=CyMask_nii.img;
figure; imagesc(CyMask)

%% also get mask from manual made from 17

WMGMCSFmask_nii = load_untouch_nii('WMGMCSFmask.nii.gz')
WMGMCSFmask=WMGMCSFmask_nii.img;
figure; imagesc(WMGMCSFmask)

%% for each, load eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs, dwmri_all.bval

% ref = b0;
% % get 4 points from ref
% %figure; imagesc(ref)
% %[x,y] = ginput(4); close all;
% 
% DWI = [];
% BVALS = [];
% BVECS = [];
% 
% DWInorm = [];
% DWInormcord = [];
% DWInormcy = [];
% 
% for i = 1:32
%     disp(i)
%     try
%     eddy = load_untouch_nii(['diff_' num2str(i) '/eddy_unwarped_images.nii.gz']);
%     eddy=eddy.img;
%     bvals = dlmread(['diff_' num2str(i) '/dwmri_all.bval']);
%     bvecs = dlmread(['diff_' num2str(i) '/dwmri_all.bvec']);
% 
%     % grab b0 from this (slice 4)
%     moving=squeeze(eddy(:,:,4,37));
% 
%     % get X points
%     %figure; imagesc(moving)
%     %[x2,y2] = ginput(4); close all;
% 
%     [selectedMovingPoints,selectedFixedPoints] = cpselect(mat2gray(moving),mat2gray(ref),'Wait',true)
%     t_concord = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'nonreflectivesimilarity');
%     disp(['scale: ' num2str(t_concord.T(1,1))])
%     disp(['rot: ' num2str(acos(t_concord.T(1,2)))])
%     % ImageRegistrationApp(moving_image,fixed_image); 
%     %ImageRegistrationApp(moving,ref); 
%     %[optimizer, metric] = imregconfig('monomodal')
%     %optimizer.MaximumIterations = 300;
%     %tform = imregtform(moving.*single(CordMask),ref.*single(CordMask),'rigid',optimizer, metric); 
%     
%     ortho_ref = imref2d(size(ref)); %relate intrinsic and world coordinates
%     b0_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
%     b0_registered2 = imwarp(moving,tform,'OutputView',ortho_ref);
%     %figure
%     %imshowpair(b0_registered,ref,'checkerboard')
%     
% 
%     % ICP to #17, apply to b0 and all dwi from forward PE (minus b=5)
%     sesDWI = [];
%     sesDWInorm = [];
%     sesDWInormcord = [];
%     sesDWInormcy = [];
%     sesBVALS=[];
%     sesBVECS=[];
%     
%     DWI = cat(3,DWI,b0_registered);
%     sesDWI = cat(3,sesDWI,b0_registered);
%     
%     BVALS = [BVALS bvals(37)];
%     BVECS = [BVECS bvecs(:,37)];
%     sesBVALS = [sesBVALS bvals(37)];
%     sesBVECS = [sesBVECS bvecs(:,37)];
%     
%     % in same loop, normalize to b0. b0 will have value of 1 everywhere. Mask with mask
%     b0_registered_norm = b0_registered; 
%     b0_registered_norm(WMGMCSFmask==1) = 1;
%     b0_registered_norm(WMGMCSFmask~=1) = 0;
%     DWInorm = cat(3,DWInorm,b0_registered_norm);
%     sesDWInorm = cat(3,sesDWInorm,b0_registered_norm);
%     
%     b0_registered_norm = b0_registered; 
%     b0_registered_norm(CordMask==1) = 1;
%     b0_registered_norm(CordMask~=1) = 0;
%     DWInormcord = cat(3,DWInormcord,b0_registered_norm);
%     sesDWInormcord = cat(3,sesDWInormcord,b0_registered_norm);
%     
%     b0_registered_norm = b0_registered; 
%     b0_registered_norm(CyMask==1) = 1;
%     b0_registered_norm(CyMask~=1) = 0;
%     DWInormcy = cat(3,DWInormcy,b0_registered_norm);
%     sesDWInormcy = cat(3,sesDWInormcy,b0_registered_norm);
%     
% 
%     % cat to image, cat bvals, cat bvecs and save as all
%     for j=1:18
%         moving=squeeze(eddy(:,:,4,j));
%         dwi_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
%         DWI = cat(3,DWI,dwi_registered);
%         sesDWI = cat(3,sesDWI,dwi_registered);
%         
%         BVALS = [BVALS bvals(j)];
%         BVECS = [BVECS bvecs(:,j)];
%         sesBVALS = [sesBVALS bvals(j)];
%         sesBVECS = [sesBVECS bvecs(:,j)];
%         
%         dwi_norm = dwi_registered./b0_registered;
%         dwi_norm(WMGMCSFmask~=1)=0;
%         DWInorm = cat(3,DWInorm,dwi_norm);
%         sesDWInorm = cat(3,sesDWInorm,dwi_norm);
%         
%         dwi_normcord = dwi_registered./b0_registered;
%         dwi_normcord(CordMask~=1)=0;
%         DWInormcord = cat(3,DWInormcord,dwi_normcord);
%         sesDWInormcord = cat(3,sesDWInormcord,dwi_normcord);
%         
%         dwi_normcy = dwi_registered./b0_registered;
%         dwi_normcy(CyMask~=1)=0;
%         DWInormcy = cat(3,DWInormcy,dwi_normcy);
%         sesDWInormcy = cat(3,sesDWInormcy,dwi_normcy);
%     end
%     
%     for j=20:36
%         moving=squeeze(eddy(:,:,4,j));
%         dwi_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
%         DWI = cat(3,DWI,dwi_registered);
%         sesDWI = cat(3,sesDWI,dwi_registered);
%         
%         BVALS = [BVALS bvals(j)];
%         BVECS = [BVECS bvecs(:,j)];
%         sesBVALS = [sesBVALS bvals(j)];
%         sesBVECS = [sesBVECS bvecs(:,j)];
%         
%         dwi_norm = dwi_registered./b0_registered;
%         dwi_norm(WMGMCSFmask~=1)=0;
%         DWInorm = cat(3,DWInorm,dwi_norm);
%         sesDWInorm = cat(3,sesDWInorm,dwi_norm);
%         
%         dwi_normcord = dwi_registered./b0_registered;
%         dwi_normcord(CordMask~=1)=0;
%         DWInormcord = cat(3,DWInormcord,dwi_normcord);
%         sesDWInormcord = cat(3,sesDWInormcord,dwi_normcord);
%         
%         dwi_normcy = dwi_registered./b0_registered;
%         dwi_normcy(CyMask~=1)=0;
%         DWInormcy = cat(3,DWInormcy,dwi_normcy);
%         sesDWInormcy = cat(3,sesDWInormcy,dwi_normcy);
%     end
%     
%     % save ses
%     sesDWI = reshape(sesDWI,[size(sesDWI,1) size(sesDWI,2) 1 size(sesDWI,3)]);
%     nii.img = sesDWI;
%     nii.hdr.dime.dim(2:5)=size(sesDWI);
%     save_untouch_nii(nii,['diff_' num2str(i) '/sesDWI.nii.gz']);
%     
%     sesDWInorm = reshape(sesDWInorm,[size(sesDWInorm,1) size(sesDWInorm,2) 1 size(sesDWInorm,3)]);
%     nii.img = sesDWInorm;
%     nii.hdr.dime.dim(2:5)=size(sesDWInorm);
%     save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInorm.nii.gz']);
%     
%     sesDWInormcord = reshape(sesDWInormcord,[size(sesDWInormcord,1) size(sesDWInormcord,2) 1 size(sesDWInormcord,3)]);
%     nii.img = sesDWInormcord;
%     nii.hdr.dime.dim(2:5)=size(sesDWInormcord);
%     save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInormcord.nii.gz']);
%     
%     sesDWInormcy = reshape(sesDWInormcy,[size(sesDWInormcy,1) size(sesDWInormcy,2) 1 size(sesDWInormcy,3)]);
%     nii.img = sesDWInormcy;
%     nii.hdr.dime.dim(2:5)=size(sesDWInormcy);
%     save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInormcy.nii.gz']);
% 
%     dlmwrite(['diff_' num2str(i) '/sesBVALS.bval'],sesBVALS,'delimiter',' ')
%     dlmwrite(['diff_' num2str(i) '/sesBVECS.bvec'],sesBVECS,'delimiter',' ')
% 
%     catch
%         disp('DNE')
%     end
% end
% 
% %% save all
% dlmwrite(['BVALS.bval'],BVALS,'delimiter',' ')
% dlmwrite(['BVECS.bvec'],BVECS,'delimiter',' ')

%% loop and save just b0's into stack
b0_all=[];
for i = 1:32
    disp(i)
    try
    nii = load_untouch_nii(['diff_' num2str(i) '/sesDWIcy.nii.gz']);
    b0=nii.img; b0=squeeze(b0(:,:,:,1));
    b0_all = cat(3,b0_all,b0);
    catch
        disp('DNE')
    end
end

b0s = reshape(b0_all,[size(b0_all,1) size(b0_all,2) 1 size(b0_all,3)]);
nii.img = b0s;
nii.hdr.dime.dim(2:5)=size(b0s);
save_untouch_nii(nii,['b0_all.nii.gz']);

%% loop and save just b0's into stack
b0_all=[];
for i = 1:32
    disp(i)
    try
    nii = load_untouch_nii(['diff_' num2str(i) '/sesDWIcy.nii.gz']);
    b0=nii.img; 
    b0(:,:,:,1)=[];
    b0=mean(b0,4);
    
    b0_all = cat(3,b0_all,b0);
    catch
        disp('DNE')
    end
end

b0s = reshape(b0_all,[size(b0_all,1) size(b0_all,2) 1 size(b0_all,3)]);
nii.img = b0s;
nii.hdr.dime.dim(2:5)=size(b0s);
save_untouch_nii(nii,['meanDWI_all.nii.gz']);
  
%% DWI
% for each, load eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs, dwmri_all.bval
nii = load_untouch_nii('diff_10/eddy_unwarped_images.nii.gz')
ref=nii.img; 
ref(:,:,:,[19 37 38 39 40])=[];
ref=squeeze(ref(:,:,4,:));
ref=mean(ref,3);

nii.img = ref;
nii.hdr.dime.dim(4)=1;
nii.hdr.dime.dim(5)=1;
save_untouch_nii(nii,'eddy_meanDWI.nii.gz')

DWI = [];
DWIcy = [];
BVALS = [];
BVECS = [];

DWInorm = [];
DWInormcord = [];
DWInormcy = [];


for i = 1:32
    disp(i)
    try
    eddy = load_untouch_nii(['diff_' num2str(i) '/eddy_unwarped_images.nii.gz']);
    eddy=eddy.img;
    bvals = dlmread(['diff_' num2str(i) '/dwmri_all.bval']);
    bvecs = dlmread(['diff_' num2str(i) '/dwmri_all.bvec']);

    dwis=eddy;
    dwis(:,:,:,[19 37 38 39 40])=[];
    
    dwis=squeeze(dwis(:,:,4,:));
    dwis=mean(dwis,3);

    %b0(:,:,:,1)=[];
    %b0=mean(b0,4);
    % grab b0 from this (slice 4)
    %moving=squeeze(eddy(:,:,4,37));
    moving = dwis;
    
    % get X points
    %figure; imagesc(moving)
    %[x2,y2] = ginput(4); close all;

    [selectedMovingPoints,selectedFixedPoints] = cpselect(mat2gray(moving),mat2gray(ref),'Wait',true);
    t_concord = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'nonreflectivesimilarity');
    disp(['scale: ' num2str(t_concord.T(1,1))])
    disp(['rot: ' num2str(asind(t_concord.T(2,1)))])
    % ImageRegistrationApp(moving_image,fixed_image); 
    %ImageRegistrationApp(moving,ref); 
    %[optimizer, metric] = imregconfig('monomodal')
    %optimizer.MaximumIterations = 300;
    %tform = imregtform(moving.*single(CordMask),ref.*single(CordMask),'rigid',optimizer, metric); 
    
    ortho_ref = imref2d(size(ref)); %relate intrinsic and world coordinates
    
    disp('b0')
    % grab b0 from this (slice 4)
    moving=squeeze(eddy(:,:,4,37));

    b0_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
    %b0_registered2 = imwarp(moving,tform,'OutputView',ortho_ref);
    %figure
    %imshowpair(b0_registered,ref,'checkerboard')
    

    % ICP to #17, apply to b0 and all dwi from forward PE (minus b=5)
    sesDWI = [];
    sesDWIcy = [];
    sesDWInorm = [];
    sesDWInormcord = [];
    sesDWInormcy = [];
    sesBVALS=[];
    sesBVECS=[];
    
    DWI = cat(3,DWI,b0_registered);
    sesDWI = cat(3,sesDWI,b0_registered);
    
    BVALS = [BVALS bvals(37)];
    BVECS = [BVECS bvecs(:,37)];
    sesBVALS = [sesBVALS bvals(37)];
    sesBVECS = [sesBVECS bvecs(:,37)];
    
    b0_registered_msk = b0_registered;
    b0_registered_msk(CyMask~=1) = 0;
    DWIcy = cat(3,DWIcy,b0_registered_msk);
    sesDWIcy = cat(3,sesDWIcy,b0_registered_msk);
    
    % in same loop, normalize to b0. b0 will have value of 1 everywhere. Mask with mask
    b0_registered_norm = b0_registered; 
    b0_registered_norm(WMGMCSFmask==1) = 1;
    b0_registered_norm(WMGMCSFmask~=1) = 0;
    DWInorm = cat(3,DWInorm,b0_registered_norm);
    sesDWInorm = cat(3,sesDWInorm,b0_registered_norm);
    
    b0_registered_norm = b0_registered; 
    b0_registered_norm(CordMask==1) = 1;
    b0_registered_norm(CordMask~=1) = 0;
    DWInormcord = cat(3,DWInormcord,b0_registered_norm);
    sesDWInormcord = cat(3,sesDWInormcord,b0_registered_norm);
    
    b0_registered_norm = b0_registered; 
    b0_registered_norm(CyMask==1) = 1;
    b0_registered_norm(CyMask~=1) = 0;
    DWInormcy = cat(3,DWInormcy,b0_registered_norm);
    sesDWInormcy = cat(3,sesDWInormcy,b0_registered_norm);
    
    disp('DWIs')
    % cat to image, cat bvals, cat bvecs and save as all
    for j=1:18
        moving=squeeze(eddy(:,:,4,j));
        dwi_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
        DWI = cat(3,DWI,dwi_registered);
        sesDWI = cat(3,sesDWI,dwi_registered);
        
        dwi_registered_msk = dwi_registered;
        dwi_registered_msk(CyMask~=1) = 0;
        DWIcy = cat(3,DWIcy,dwi_registered_msk);
        sesDWIcy = cat(3,sesDWIcy,dwi_registered_msk);
  
        BVALS = [BVALS bvals(j)];
        BVECS = [BVECS bvecs(:,j)];
        sesBVALS = [sesBVALS bvals(j)];
        sesBVECS = [sesBVECS bvecs(:,j)];
        
        dwi_norm = dwi_registered./b0_registered;
        dwi_norm(WMGMCSFmask~=1)=0;
        DWInorm = cat(3,DWInorm,dwi_norm);
        sesDWInorm = cat(3,sesDWInorm,dwi_norm);
        
        dwi_normcord = dwi_registered./b0_registered;
        dwi_normcord(CordMask~=1)=0;
        DWInormcord = cat(3,DWInormcord,dwi_normcord);
        sesDWInormcord = cat(3,sesDWInormcord,dwi_normcord);
        
        dwi_normcy = dwi_registered./b0_registered;
        dwi_normcy(CyMask~=1)=0;
        DWInormcy = cat(3,DWInormcy,dwi_normcy);
        sesDWInormcy = cat(3,sesDWInormcy,dwi_normcy);
    end
    
    for j=20:36
        moving=squeeze(eddy(:,:,4,j));
        dwi_registered = imwarp(moving,t_concord,'OutputView',ortho_ref);
        DWI = cat(3,DWI,dwi_registered);
        sesDWI = cat(3,sesDWI,dwi_registered);
        
        dwi_registered_msk = dwi_registered;
        dwi_registered_msk(CyMask~=1) = 0;
        DWIcy = cat(3,DWIcy,dwi_registered_msk);
        sesDWIcy = cat(3,sesDWIcy,dwi_registered_msk);
  
        
        BVALS = [BVALS bvals(j)];
        BVECS = [BVECS bvecs(:,j)];
        sesBVALS = [sesBVALS bvals(j)];
        sesBVECS = [sesBVECS bvecs(:,j)];
        
        dwi_norm = dwi_registered./b0_registered;
        dwi_norm(WMGMCSFmask~=1)=0;
        DWInorm = cat(3,DWInorm,dwi_norm);
        sesDWInorm = cat(3,sesDWInorm,dwi_norm);
        
        dwi_normcord = dwi_registered./b0_registered;
        dwi_normcord(CordMask~=1)=0;
        DWInormcord = cat(3,DWInormcord,dwi_normcord);
        sesDWInormcord = cat(3,sesDWInormcord,dwi_normcord);
        
        dwi_normcy = dwi_registered./b0_registered;
        dwi_normcy(CyMask~=1)=0;
        DWInormcy = cat(3,DWInormcy,dwi_normcy);
        sesDWInormcy = cat(3,sesDWInormcy,dwi_normcy);
    end
    
    disp('SAVING')
    % save ses
    sesDWI = reshape(sesDWI,[size(sesDWI,1) size(sesDWI,2) 1 size(sesDWI,3)]);
    nii.img = sesDWI;
    nii.hdr.dime.dim(2:5)=size(sesDWI);
    save_untouch_nii(nii,['diff_' num2str(i) '/sesDWI.nii.gz']);
    
    % save ses
    sesDWIcy = reshape(sesDWIcy,[size(sesDWIcy,1) size(sesDWIcy,2) 1 size(sesDWIcy,3)]);
    nii.img = sesDWIcy;
    nii.hdr.dime.dim(2:5)=size(sesDWIcy);
    save_untouch_nii(nii,['diff_' num2str(i) '/sesDWIcy.nii.gz']);
    
    sesDWInorm = reshape(sesDWInorm,[size(sesDWInorm,1) size(sesDWInorm,2) 1 size(sesDWInorm,3)]);
    nii.img = sesDWInorm;
    nii.hdr.dime.dim(2:5)=size(sesDWInorm);
    save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInorm.nii.gz']);
    
    sesDWInormcord = reshape(sesDWInormcord,[size(sesDWInormcord,1) size(sesDWInormcord,2) 1 size(sesDWInormcord,3)]);
    nii.img = sesDWInormcord;
    nii.hdr.dime.dim(2:5)=size(sesDWInormcord);
    save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInormcord.nii.gz']);
    
    sesDWInormcy = reshape(sesDWInormcy,[size(sesDWInormcy,1) size(sesDWInormcy,2) 1 size(sesDWInormcy,3)]);
    nii.img = sesDWInormcy;
    nii.hdr.dime.dim(2:5)=size(sesDWInormcy);
    save_untouch_nii(nii,['diff_' num2str(i) '/sesDWInormcy.nii.gz']);

    dlmwrite(['diff_' num2str(i) '/sesBVALS.bval'],sesBVALS,'delimiter',' ')
    dlmwrite(['diff_' num2str(i) '/sesBVECS.bvec'],sesBVECS,'delimiter',' ')
    
    disp('SUCCESS')
    close all
    catch
        disp('fail or DNE')
        close all
    end
end


%%

% quality check

% redo bad ones

% stack all and resave

DWI = [];
DWInorm = [];

for i = 1:32
    disp(i)
    try
        nii = load_untouch_nii(['diff_' num2str(i) '/sesDWIcy.nii.gz']);
        dwi=nii.img; dwi = squeeze(dwi);
        DWI = cat(3,DWI,dwi);
        
        nii2 = load_untouch_nii(['diff_' num2str(i) '/sesDWInormcy.nii.gz']);
        dwi2=nii2.img; dwi2 = squeeze(dwi2);
        DWInorm = cat(3,DWInorm,dwi2);
        
    catch
        disp('DNE')
    end
end

mask_nii = load_untouch_nii('meanDWI_CordMask.nii.gz')
mask=mask_nii.img;
figure; imagesc(mask)


 for i = 1:size(DWI,3)
     temp = DWI(:,:,i);
    
     blah(i) = mean(temp(mask==1));
 end
 
for i = 1:size(DWInorm,3)
     temp = DWInorm(:,:,i);
     blah(i) = mean(temp(mask==1));
 end
     

DWI = reshape(DWI,[size(DWI,1) size(DWI,2) 1 size(DWI,3)]);
nii.img = DWI;
nii.hdr.dime.dim(2:5)=size(DWI);
save_untouch_nii(nii,['dwmri.nii.gz']);
    
DWInorm = reshape(DWInorm,[size(DWInorm,1) size(DWInorm,2) 1 size(DWInorm,3)]);
nii.img = DWInorm;
nii.hdr.dime.dim(2:5)=size(DWInorm);
save_untouch_nii(nii,['dwmri_norm.nii.gz']);
    

%%
a=[17	6	24.11	61.79	218	43
2	6	28.51	61.67	260	48
21	6	40	61.74	374	59
12	6	47.18	61.76	445	66
5	6	55.62	61.84	530	75
24	6	65.82	61.86	631	85
11	6	71.38	61.86	686	91
23	6	84.67	61.81	816	104
25	10	24.57	61.57	577	48
26	10	28.62	61.58	687	52
10	10	40.3	61.72	1010	64
30	10	47.11	61.84	1202	70
18	10	56.01	61.79	1443	79
15	10	66.08	61.73	1718	89
27	10	71.71	61.82	1876	95
31	10	84.72	61.85	2233	108
3	15	26.05	61.35	1218	54
22	15	28.68	61.58	1449	57
1	15	40.48	61.89	2192	69
19	15	47.28	61.84	2608	76
6	15	56.2	61.89	3166	84
32	15	66.01	61.81	3768	94
9	15	71.94	61.89	4140	100
20	19.2	30.53	61.84	2375	63
14	20.42	31.8	61.87	2826	65
4	22	41.19	61.7	4470	76]

scheme = [];
for i=1:32
    try
       index = find(a(:,1) == i);
       vals = a(index,2:6);
       scheme = [scheme ; vals];
    catch
        disp('DNE')
    end
end

       
scheme = [];
for i=1:32
    try
       index = find(a(:,1) == i);
       vals = a(index,2:6);
       
       for j = 1:36
           if j==1
               val = vals;
               val(3) = 0;
               val(4) = 0;
           else
               val = vals;
           end
           
           scheme = [scheme ; val];
       end
       
    catch
        disp('DNE')
    end
end

dlmwrite('scheme.scheme',scheme,'delimiter',' ')               



