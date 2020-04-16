% plot signal as a function of orientation 
% colors designate delta
% symbols designate different b-valeus within delta

addpath(genpath('/Volumes/schillkg/MATLAB/anderson'))

scheme = dlmread('DWMRI_all/SingleVoxelSignals_norm/scheme.scheme')

dirs = dlmread('DWMRI_all/SingleVoxelSignals_norm/BVECS.bvec')

for i = 1:length(dirs)
    currdir=dirs(:,i);
    if sum(currdir)==0
        ang=pi/2;
    else
        ang = angleBetw2vect2(currdir(1),currdir(2),currdir(3),0,0,1,1);
    end
    direction(i) = cos(ang);
end

% colors = {'r','g','b','m'}
colors = [0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330;  0.9290 0.6940 0.1250; 0 0.4470 0.7410]
    
%symbols = {'o','+','*','x','s','d','^','p'}
symbols = {'o','<','h','>','s','d','^','p'}


%% LateralCST
signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/LateralCST.txt')

figure; hold on; ms=8;
Aind = 1;Bind = 1;Cind = 1;Dind = 1;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   if median(scheme(indices,1)) == 6
       c = colors(1,:); s = symbols{Aind}; Aind = Aind+1;
   elseif median(scheme(indices,1)) == 10
       c = colors(2,:); s = symbols{Bind}; Bind = Bind+1;
   elseif median(scheme(indices,1)) == 15
       c = colors(3,:); s = symbols{Cind}; Cind = Cind+1;    
   else
       c = colors(4,:); s = symbols{Dind}; Dind = Dind+1;
   end
   plot(direction(indices), signal(indices),'color',c,'marker',s,'LineStyle','none','MarkerSize',ms,'MarkerFaceColor',c,'MarkerEdgeColor','k')
end

xlabel('Direction'); ylabel('Signal'); set(gca,'FontSize',24); 
whitebg(gcf,[1 1 1]); set(gcf,'Color',[1,1,1]); box on; grid on;
fig = gcf; fig.Color = 'white'; fig.InvertHardcopy = 'off'; saveas(gcf,['Figures/F1_LateralCST_norm.png'])
pause(1);
close all;

%%  VentralHorn
signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/VentralHorn.txt')

figure; hold on; ms=8;
Aind = 1;Bind = 1;Cind = 1;Dind = 1;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   if median(scheme(indices,1)) == 6
       c = colors(1,:); s = symbols{Aind}; Aind = Aind+1;
   elseif median(scheme(indices,1)) == 10
       c = colors(2,:); s = symbols{Bind}; Bind = Bind+1;
   elseif median(scheme(indices,1)) == 15
       c = colors(3,:); s = symbols{Cind}; Cind = Cind+1;    
   else
       c = colors(4,:); s = symbols{Dind}; Dind = Dind+1;
   end
   plot(direction(indices), signal(indices),'color',c,'marker',s,'LineStyle','none','MarkerSize',ms,'MarkerFaceColor',c,'MarkerEdgeColor','k')
end

xlabel('Direction'); ylabel('Signal'); set(gca,'FontSize',24); 
whitebg(gcf,[1 1 1]); set(gcf,'Color',[1,1,1]); box on; grid on;
fig = gcf; fig.Color = 'white'; fig.InvertHardcopy = 'off'; saveas(gcf,['Figures/F1_VentralHorn_norm.png'])
pause(1);
close all;

%% VentralCST
signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/VentralCST.txt')

figure; hold on; ms=8;
Aind = 1;Bind = 1;Cind = 1;Dind = 1;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   if median(scheme(indices,1)) == 6
       c = colors(1,:); s = symbols{Aind}; Aind = Aind+1;
   elseif median(scheme(indices,1)) == 10
       c = colors(2,:); s = symbols{Bind}; Bind = Bind+1;
   elseif median(scheme(indices,1)) == 15
       c = colors(3,:); s = symbols{Cind}; Cind = Cind+1;    
   else
       c = colors(4,:); s = symbols{Dind}; Dind = Dind+1;
   end
   plot(direction(indices), signal(indices),'color',c,'marker',s,'LineStyle','none','MarkerSize',ms,'MarkerFaceColor',c,'MarkerEdgeColor','k')
end

xlabel('Direction'); ylabel('Signal'); set(gca,'FontSize',24); 
whitebg(gcf,[1 1 1]); set(gcf,'Color',[1,1,1]); box on; grid on;
fig = gcf; fig.Color = 'white'; fig.InvertHardcopy = 'off'; saveas(gcf,['Figures/F1_VentralCST_norm.png'])
pause(1);
close all;

%% SpinalLemniscus
signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/SpinalLemniscus.txt')

figure; hold on; ms=8;
Aind = 1;Bind = 1;Cind = 1;Dind = 1;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   if median(scheme(indices,1)) == 6
       c = colors(1,:); s = symbols{Aind}; Aind = Aind+1;
   elseif median(scheme(indices,1)) == 10
       c = colors(2,:); s = symbols{Bind}; Bind = Bind+1;
   elseif median(scheme(indices,1)) == 15
       c = colors(3,:); s = symbols{Cind}; Cind = Cind+1;    
   else
       c = colors(4,:); s = symbols{Dind}; Dind = Dind+1;
   end
   plot(direction(indices), signal(indices),'color',c,'marker',s,'LineStyle','none','MarkerSize',ms,'MarkerFaceColor',c,'MarkerEdgeColor','k')
end

xlabel('Direction'); ylabel('Signal'); set(gca,'FontSize',24); 
whitebg(gcf,[1 1 1]); set(gcf,'Color',[1,1,1]); box on; grid on;
fig = gcf; fig.Color = 'white'; fig.InvertHardcopy = 'off'; saveas(gcf,['Figures/F1_SpinalLemniscus_norm.png'])
pause(1);
close all;

%% FasciulusCuneatus
signal = dlmread('DWMRI_all/SingleVoxelSignals_norm/FasciulusCuneatus.txt')

figure; hold on; ms=8;
Aind = 1;Bind = 1;Cind = 1;Dind = 1;
for i = 1:26 % schemes
   indices = i*36-35:i*36;
   if median(scheme(indices,1)) == 6
       c = colors(1,:); s = symbols{Aind}; Aind = Aind+1;
   elseif median(scheme(indices,1)) == 10
       c = colors(2,:); s = symbols{Bind}; Bind = Bind+1;
   elseif median(scheme(indices,1)) == 15
       c = colors(3,:); s = symbols{Cind}; Cind = Cind+1;    
   else
       c = colors(4,:); s = symbols{Dind}; Dind = Dind+1;
   end
   plot(direction(indices), signal(indices),'color',c,'marker',s,'LineStyle','none','MarkerSize',ms,'MarkerFaceColor',c,'MarkerEdgeColor','k')
end

xlabel('Direction'); ylabel('Signal'); set(gca,'FontSize',24); 
whitebg(gcf,[1 1 1]); set(gcf,'Color',[1,1,1]); box on; grid on;
fig = gcf; fig.Color = 'white'; fig.InvertHardcopy = 'off'; saveas(gcf,['Figures/F1_FasciulusCuneatus_norm.png'])
pause(1);
close all;

   
       
       
   
   
   
   
   
   
   
   
   
   
    