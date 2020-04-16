% v2
% added BIC
% removed G3

C = {};

voxel = 'FC'

signal = 'SingleVoxelSignals/FasciulusCuneatus.txt'
sig = dlmread(signal);

%% intra_extra_rest
% two compartment
index = 1;
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        json_name = ['intra_extra_rest' filesep voxel filesep intra{1} '_' extra{1} '.json'];
        json = jsondecode(fileread(json_name));
        disp(json)
        K = length(fieldnames(json))-2+1;
        n = length(sig);
        txt_name = ['intra_extra_rest' filesep voxel filesep intra{1} '_' extra{1} '_signal.txt'];
        txt = dlmread(txt_name);
        C{index,1} = [intra{1} '_' extra{1}];
        error=txt(:)-sig(:);
        MSE = mean(error.^2);
        RSS = sum(error.^2);
        AICresearchgate = n*log(RSS/n) + 2*K;
        AICwiki = 2*K + n*log(RSS);
        AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
        BIC = n*log(RSS/n) + K*log(n);
        C{index,2} = MSE;
        C{index,3} = AICresearchgate;
        C{index,4} = AICwiki;
        C{index,5} = BIC;
        C{index,6} = K;
        C{index,7} = json.partial_volume_0;
        C{index,8} = json.partial_volume_1;
        C{index,9} = NaN;
        
           C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            
        switch intra{1}
            case 'C1'
                C{index,10} = json.C1Stick_1_mu(1);
                C{index,11} = json.C1Stick_1_mu(2);
                C{index,12} = json.C1Stick_1_lambda_par(1);
                C{index,13} = NaN;
            case 'C2'
                C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par(1);
                C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
            case 'C3'
                C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par(1);
                C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
            case 'C4'
                C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par(1);
                C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
        end
        
        switch extra{1}
            case 'G1'
                C{index,14} = json.G1Ball_1_lambda_iso;
                C{index,15} = NaN;
                C{index,16} = NaN;
            case 'G2'
                C{index,14} = json.G2Zeppelin_1_lambda_par;
                C{index,15} = json.G2Zeppelin_1_lambda_perp;
                C{index,16} = NaN;
            case 'G3'
                C{index,14} = json.G3TemporalZeppelin_1_lambda_par;
                C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                C{index,16} = json.G3TemporalZeppelin_1_A;
        end
        
        C{index,17} = NaN;
        C{index,18} = NaN;  
        
        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
        json_name = ['intra_extra_rest' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '.json'];
        json = jsondecode(fileread(json_name));
        disp(json)
        K = length(fieldnames(json))-2+1;
        n = length(sig);
        txt_name = ['intra_extra_rest' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_signal.txt'];
        txt = dlmread(txt_name);
        C{index,1} = [intra{1} '_' extra{1} '_' rest{1}];
        error=txt(:)-sig(:);
        MSE = mean(error.^2);
        RSS = sum(error.^2);
        AICresearchgate = n*log(RSS/n) + 2*K;
        AICwiki = 2*K + n*log(RSS);
        AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
        BIC = n*log(RSS/n) + K*log(n);
        C{index,2} = MSE;
        C{index,3} = AICresearchgate;
        C{index,4} = AICwiki;
        C{index,5} = BIC;
        C{index,6} = K;
        C{index,7} = json.partial_volume_0;
        C{index,8} = json.partial_volume_1;
        C{index,9} = json.partial_volume_2;
        
           C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            
        switch intra{1}
            case 'C1'
                C{index,10} = json.C1Stick_1_mu(1);
                C{index,11} = json.C1Stick_1_mu(2);
                C{index,12} = json.C1Stick_1_lambda_par(1);
                C{index,13} = NaN;
            case 'C2'
                C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par(1);
                C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
            case 'C3'
                C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par(1);
                C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
            case 'C4'
                C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par(1);
                C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
        end
        
        switch extra{1}
            case 'G1'
                C{index,14} = json.G1Ball_1_lambda_iso;
                C{index,15} = NaN;
                C{index,16} = NaN;
            case 'G2'
                C{index,14} = json.G2Zeppelin_1_lambda_par;
                C{index,15} = json.G2Zeppelin_1_lambda_perp;
                C{index,16} = NaN;
            case 'G3'
                C{index,14} = json.G3TemporalZeppelin_1_lambda_par;
                C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                C{index,16} = json.G3TemporalZeppelin_1_A;
        end
        
         switch rest{1}
            case 'S1'
                C{index,17} = NaN;
            case 'S2'
                C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
            case 'S4'
                C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
         end
              
         C{index,18} = NaN;  
         
        index = index+1;
        end % rest
    end % extra
end % intra

%% par_fixed
% two compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        json_name = ['intra_extra_rest_par_fixed' filesep voxel filesep intra{1} '_' extra{1} '_par_fix.json'];
        json = jsondecode(fileread(json_name));
        disp(json)
        K = length(fieldnames(json))-2+1;
        n = length(sig);
        txt_name = ['intra_extra_rest_par_fixed' filesep voxel filesep intra{1} '_' extra{1} '_par_fix_signal.txt'];
        txt = dlmread(txt_name);
        C{index,1} = [intra{1} '_' extra{1} '_fix'];
        error=txt(:)-sig(:);
        MSE = mean(error.^2);
        RSS = sum(error.^2);
        AICresearchgate = n*log(RSS/n) + 2*K;
        AICwiki = 2*K + n*log(RSS);
        AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
        BIC = n*log(RSS/n) + K*log(n);
        
        C{index,2} = MSE;
        C{index,3} = AICresearchgate;
        C{index,4} = AICwiki;
        C{index,5} = BIC;
        C{index,6} = K;
        C{index,7} = json.partial_volume_0;
        C{index,8} = json.partial_volume_1;
        C{index,9} = NaN;

           C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            
        switch intra{1}
            case 'C1'
                C{index,10} = json.C1Stick_1_mu(1);
                C{index,11} = json.C1Stick_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = NaN;
            case 'C2'
                C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
            case 'C3'
                C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
            case 'C4'
                C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
        end
        
        switch extra{1}
            case 'G1'
                C{index,14} = 1.7e-09;
                C{index,15} = NaN;
                C{index,16} = NaN;
            case 'G2'
                C{index,14} = 1.7e-09;
                C{index,15} = json.G2Zeppelin_1_lambda_perp;
                C{index,16} = NaN;
            case 'G3'
                C{index,14} = 1.7e-09;
                C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                C{index,16} = json.G3TemporalZeppelin_1_A;
        end
        
        C{index,17} = NaN;
        C{index,18} = NaN;    
        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
        json_name = ['intra_extra_rest_par_fixed' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_par_fix.json'];
        json = jsondecode(fileread(json_name));
        disp(json)
        K = length(fieldnames(json))-2+1;
        n = length(sig);
        txt_name = ['intra_extra_rest_par_fixed' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_par_fix_signal.txt'];
        txt = dlmread(txt_name);
        C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_fix'];
        error=txt(:)-sig(:);
        MSE = mean(error.^2);
        RSS = sum(error.^2);
        AICresearchgate = n*log(RSS/n) + 2*K;
        AICwiki = 2*K + n*log(RSS);
        AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
        BIC = n*log(RSS/n) + K*log(n);
        
        C{index,2} = MSE;
        C{index,3} = AICresearchgate;
        C{index,4} = AICwiki;
        C{index,5} = BIC;
        C{index,6} = K;
        C{index,7} = json.partial_volume_0;
        C{index,8} = json.partial_volume_1;
        C{index,9} = json.partial_volume_2;

           C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            
        switch intra{1}
            case 'C1'
                C{index,10} = json.C1Stick_1_mu(1);
                C{index,11} = json.C1Stick_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = NaN;
            case 'C2'
                C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
            case 'C3'
                C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
            case 'C4'
                C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                C{index,12} = 1.7e-09;
                C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
        end
        
        switch extra{1}
            case 'G1'
                C{index,14} = 1.7e-09;
                C{index,15} = NaN;
                C{index,16} = NaN;
            case 'G2'
                C{index,14} = 1.7e-09;
                C{index,15} = json.G2Zeppelin_1_lambda_perp;
                C{index,16} = NaN;
            case 'G3'
                C{index,14} = 1.7e-09;
                C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                C{index,16} = json.G3TemporalZeppelin_1_A;
        end
        
         switch rest{1}
            case 'S1'
                C{index,17} = NaN;
            case 'S2'
                C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
            case 'S4'
                C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
         end
        C{index,18} = NaN;  
        
        index = index+1;
        end % rest
    end % extra
end % intra

%% intraextra_equal
% two compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
            json_name = ['intra_extra_rest_intraextra_equal' filesep voxel filesep intra{1} '_' extra{1} 'intra_extra_eq.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_intraextra_equal' filesep voxel filesep intra{1} '_' extra{1} 'intra_extra_eq_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_equal'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = NaN;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = json.C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = C{index,12};
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = C{index,12};
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = C{index,12};
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
            C{index,17} = NaN;
           C{index,18} = NaN;  
           
        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
            json_name = ['intra_extra_rest_intraextra_equal' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} 'intra_extra_eq.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_intraextra_equal' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} 'intra_extra_eq_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_equal'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = json.partial_volume_2;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN;
            
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = json.C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = C{index,12};
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = C{index,12};
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = C{index,12};
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
            switch rest{1}
                case 'S1'
                    C{index,17} = NaN;
                case 'S2'
                    C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
                case 'S4'
                    C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
            end
        
            C{index,18} = NaN;  
            
        index = index+1;
        end % rest
    end % extra
end % intra


%% watson
% two compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
            json_name = ['intra_extra_rest_watson' filesep voxel filesep intra{1} '_' extra{1} '_watson.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-1+1; % changed because no partial volume given
            n = length(sig);
            txt_name = ['intra_extra_rest_watson' filesep voxel filesep intra{1} '_' extra{1} '_watson_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_watson'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.SD1WatsonDistributed_1_partial_volume_0;
            C{index,8} = 1-json.SD1WatsonDistributed_1_partial_volume_0;
            C{index,9} = NaN;
            
            C{index,17} = NaN; % size/diameter of restricted
            C{index,18} = json.SD1WatsonDistributed_1_SD1Watson_1_odi;  
            C{index,19} = NaN; 
            C{index,20} = NaN;
            switch intra{1}
                case 'C1'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C2CylinderStejskalTannerApproximation_1_;
                    C{index,13} = json.SD1WatsonDistributed_1_C2CylinderStejskalTannerApproximation__1;
                case 'C3'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C3CylinderCallaghanApproximation_1_lambd;
                    C{index,13} = json.SD1WatsonDistributed_1_C3CylinderCallaghanApproximation_1_diame;
                case 'C4'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C4CylinderGaussianPhaseApproximation_1_l;
                    C{index,13} = json.SD1WatsonDistributed_1_C4CylinderGaussianPhaseApproximation_1_d;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.SD1WatsonDistributed_1_G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.SD1WatsonDistributed_1_G2Zeppelin_1_lambda_par;
                    C{index,15} = json.SD1WatsonDistributed_1_G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_A;
            end
            
        index = index+1;
    end % extra
end % intra


% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
            json_name = ['intra_extra_rest_watson' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_watson.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_watson' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_watson_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_watson'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0 * json.SD1WatsonDistributed_1_partial_volume_0;
            C{index,8} = json.partial_volume_0 * (1-json.SD1WatsonDistributed_1_partial_volume_0);
            C{index,9} = json.partial_volume_1;
            
            C{index,18} = json.SD1WatsonDistributed_1_SD1Watson_1_odi;                           
            C{index,19} = NaN; 
            C{index,20} = NaN;
            switch intra{1}
                case 'C1'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C2CylinderStejskalTannerApproximation_1_;
                    C{index,13} = json.SD1WatsonDistributed_1_C2CylinderStejskalTannerApproximation__1;
                case 'C3'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C3CylinderCallaghanApproximation_1_lambd;
                    C{index,13} = json.SD1WatsonDistributed_1_C3CylinderCallaghanApproximation_1_diame;
                case 'C4'
                    C{index,10} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(1);
                    C{index,11} = json.SD1WatsonDistributed_1_SD1Watson_1_mu(2);
                    C{index,12} = json.SD1WatsonDistributed_1_C4CylinderGaussianPhaseApproximation_1_l;
                    C{index,13} = json.SD1WatsonDistributed_1_C4CylinderGaussianPhaseApproximation_1_d;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.SD1WatsonDistributed_1_G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.SD1WatsonDistributed_1_G2Zeppelin_1_lambda_par;
                    C{index,15} = json.SD1WatsonDistributed_1_G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.SD1WatsonDistributed_1_G3TemporalZeppelin_1_A;
            end
            
            switch rest{1}
                case 'S1'
                    C{index,17} = NaN;
                case 'S2'
                    C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
                case 'S4'
                    C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
            end
            
        index = index+1;
        end % rest
    end % extra
end % intra


%% bingham
% two compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
            json_name = ['intra_extra_rest_bingham' filesep voxel filesep intra{1} '_' extra{1} '_bingham.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-1+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_bingham' filesep voxel filesep intra{1} '_' extra{1} '_bingham_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_bingham'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
      BIC = n*log(RSS/n) + K*log(n);
      
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;

            C{index,7} = json.SD2BinghamDistributed_1_partial_volume_0;
            C{index,8} = 1-json.SD2BinghamDistributed_1_partial_volume_0;
            C{index,9} = NaN;
            
            C{index,17} = NaN; % size/diameter of restricted
            C{index,18} = json.SD2BinghamDistributed_1_SD2Bingham_1_odi; 
            C{index,19} = json.SD2BinghamDistributed_1_SD2Bingham_1_beta_fraction; 
            C{index,20} = json.SD2BinghamDistributed_1_SD2Bingham_1_psi; 
            
            switch intra{1}
                case 'C1'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C2CylinderStejskalTannerApproximation_1;
                    C{index,13} = json.SD2BinghamDistributed_1_C2CylinderStejskalTannerApproximation_2;
                case 'C3'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C3CylinderCallaghanApproximation_1_lamb;
                    C{index,13} = json.SD2BinghamDistributed_1_C3CylinderCallaghanApproximation_1_diam;
                case 'C4'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C4CylinderGaussianPhaseApproximation_1_;
                    C{index,13} = json.SD2BinghamDistributed_1_C4CylinderGaussianPhaseApproximation__1;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.SD2BinghamDistributed_1_G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.SD2BinghamDistributed_1_G2Zeppelin_1_lambda_par;
                    C{index,15} = json.SD2BinghamDistributed_1_G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_A;
            end
            
           
        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
            json_name = ['intra_extra_rest_bingham' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_bingham.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_bingham' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_bingham_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_bingham'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0 * json.SD2BinghamDistributed_1_partial_volume_0;
            C{index,8} = json.partial_volume_0 * (1-json.SD2BinghamDistributed_1_partial_volume_0);
            C{index,9} = json.partial_volume_1;
            
            C{index,18} = json.SD2BinghamDistributed_1_SD2Bingham_1_odi; 
            C{index,19} = json.SD2BinghamDistributed_1_SD2Bingham_1_beta_fraction; 
            C{index,20} = json.SD2BinghamDistributed_1_SD2Bingham_1_psi; 
            
            switch intra{1}
                case 'C1'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C2CylinderStejskalTannerApproximation_1;
                    C{index,13} = json.SD2BinghamDistributed_1_C2CylinderStejskalTannerApproximation_2;
                case 'C3'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C3CylinderCallaghanApproximation_1_lamb;
                    C{index,13} = json.SD2BinghamDistributed_1_C3CylinderCallaghanApproximation_1_diam;
                case 'C4'
                    C{index,10} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(1);
                    C{index,11} = json.SD2BinghamDistributed_1_SD2Bingham_1_mu(2);
                    C{index,12} = json.SD2BinghamDistributed_1_C4CylinderGaussianPhaseApproximation_1_;
                    C{index,13} = json.SD2BinghamDistributed_1_C4CylinderGaussianPhaseApproximation__1;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.SD2BinghamDistributed_1_G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.SD2BinghamDistributed_1_G2Zeppelin_1_lambda_par;
                    C{index,15} = json.SD2BinghamDistributed_1_G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.SD2BinghamDistributed_1_G3TemporalZeppelin_1_A;
            end
            
            switch rest{1}
                case 'S1'
                    C{index,17} = NaN;
                case 'S2'
                    C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
                case 'S4'
                    C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
            end
            
        index = index+1;
        end % rest
    end % extra
end % intra

%% intra
% two compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
            json_name = ['intra_extra_rest_intra_greater' filesep voxel filesep intra{1} '_' extra{1} '_intra_greater.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_intra_greater' filesep voxel filesep intra{1} '_' extra{1} '_intra_greater_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_intra'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
BIC = n*log(RSS/n) + K*log(n);

            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = NaN;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN; 
                     
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = json.C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.G1Ball_1_lambda_iso_fraction*C{index,12};
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.G2Zeppelin_1_lambda_par_fraction*C{index,12};
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.G3TemporalZeppelin_1_lambda_par_fraction*C{index,12};
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
            json_name = ['intra_extra_rest_intra_greater' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_intra_greater.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_intra_greater' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_intra_greater_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_intra'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = json.partial_volume_2;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN; 
                     
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = json.C1Stick_1_lambda_par;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = json.C4CylinderGaussianPhaseApproximation_1_lambda_par;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.G1Ball_1_lambda_iso_fraction*C{index,12};
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.G2Zeppelin_1_lambda_par_fraction*C{index,12};
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.G3TemporalZeppelin_1_lambda_par_fraction*C{index,12};
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
            switch rest{1}
                case 'S1'
                    C{index,17} = NaN;
                case 'S2'
                    C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
                case 'S4'
                    C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
            end
            
        index = index+1;
        end % rest
    end % extra
end % intra


%% extra
% tw0 compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
            json_name = ['intra_extra_rest_extra_greater' filesep voxel filesep intra{1} '_' extra{1} '_extra_greater.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_extra_greater' filesep voxel filesep intra{1} '_' extra{1} '_extra_greater_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_extra'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = NaN;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN; 
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.G2Zeppelin_1_lambda_par;
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = C{index,14}*json.C1Stick_1_lambda_par_fraction;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C2CylinderStejskalTannerApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C3CylinderCallaghanApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C4CylinderGaussianPhaseApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end

        index = index+1;
    end % extra
end % intra

% three compartment
for intra = {'C1', 'C2', 'C3', 'C4'}
    for extra = {'G1', 'G2'}
        for rest = {'S1', 'S2', 'S4'}
            json_name = ['intra_extra_rest_extra_greater' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_extra_greater.json'];
            json = jsondecode(fileread(json_name));
            disp(json)
            K = length(fieldnames(json))-2+1;
            n = length(sig);
            txt_name = ['intra_extra_rest_extra_greater' filesep voxel filesep intra{1} '_' extra{1} '_' rest{1} '_extra_greater_signal.txt'];
            txt = dlmread(txt_name);
            C{index,1} = [intra{1} '_' extra{1} '_' rest{1} '_extra'];
            error=txt(:)-sig(:);
            MSE = mean(error.^2);
            RSS = sum(error.^2);
            AICresearchgate = n*log(RSS/n) + 2*K;
            AICwiki = 2*K + n*log(RSS);
            AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
            BIC = n*log(RSS/n) + K*log(n);
            
            C{index,2} = MSE;
            C{index,3} = AICresearchgate;
            C{index,4} = AICwiki;
            C{index,5} = BIC;
            C{index,6} = K;
            
            C{index,7} = json.partial_volume_0;
            C{index,8} = json.partial_volume_1;
            C{index,9} = json.partial_volume_2;
            
            C{index,18} = NaN; 
            C{index,19} = NaN; 
            C{index,20} = NaN; 
            
            switch extra{1}
                case 'G1'
                    C{index,14} = json.G1Ball_1_lambda_iso;
                    C{index,15} = NaN;
                    C{index,16} = NaN;
                case 'G2'
                    C{index,14} = json.G2Zeppelin_1_lambda_par;
                    C{index,15} = json.G2Zeppelin_1_lambda_perp;
                    C{index,16} = NaN;
                case 'G3'
                    C{index,14} = json.G3TemporalZeppelin_1_lambda_par;
                    C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
                    C{index,16} = json.G3TemporalZeppelin_1_A;
            end
            
            switch intra{1}
                case 'C1'
                    C{index,10} = json.C1Stick_1_mu(1);
                    C{index,11} = json.C1Stick_1_mu(2);
                    C{index,12} = C{index,14}*json.C1Stick_1_lambda_par_fraction;
                    C{index,13} = NaN;
                case 'C2'
                    C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
                    C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C2CylinderStejskalTannerApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
                case 'C3'
                    C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
                    C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C3CylinderCallaghanApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
                case 'C4'
                    C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
                    C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
                    C{index,12} = C{index,14}*json.C4CylinderGaussianPhaseApproximation_1_lambda_par_fraction;
                    C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
            end
            
            switch rest{1}
                case 'S1'
                    C{index,17} = NaN;
                case 'S2'
                    C{index,17} = json.S2SphereStejskalTannerApproximation_1_diameter;
                case 'S4'
                    C{index,17} = json.S4SphereGaussianPhaseApproximation_1_diameter;
            end
            
        index = index+1;
        end % rest
    end % extra
end % intra

%% tort
% three compartment
% for intra = {'C1', 'C2', 'C3', 'C4'}
%     for extra = {'G1', 'G2'}
%        
%             json_name = ['intra_extra_rest_tort' filesep voxel filesep intra{1} '_' extra{1} '_tort.json'];
%             json = jsondecode(fileread(json_name));
%             disp(json)
%             K = length(fieldnames(json))-2+1;
%             n = length(sig);
%             txt_name = ['intra_extra_rest_tort' filesep voxel filesep intra{1} '_' extra{1} '_tort_signal.txt'];
%             txt = dlmread(txt_name);
%             C{index,1} = [intra{1} '_' extra{1} '_' rest{1}];
%             error=txt(:)-sig(:);
%             MSE = mean(error.^2);
%             RSS = sum(error.^2);
%             AICresearchgate = n*log(RSS/n) + 2*K;
%             AICwiki = 2*K + n*log(RSS);
%             AICanderson = sum((error).^2) ./ std(error).^2 +2*K;
%             BIC = n*log(RSS/n) + K*log(n);
%             C{index,2} = MSE;
%             C{index,3} = AICresearchgate;
%             C{index,4} = AICwiki;
%             C{index,5} = BIC;
%             C{index,6} = K;
%             
%             C{index,7} = json.partial_volume_0;
%             C{index,8} = json.partial_volume_1;
%             C{index,9} = NaN;
%             
%             C{index,18} = NaN; 
%             C{index,19} = NaN; 
%             C{index,20} = NaN; 
%             
%             switch intra{1}
%                 case 'C1'
%                     C{index,10} = json.C1Stick_1_mu(1);
%                     C{index,11} = json.C1Stick_1_mu(2);
%                     C{index,12} = json.C1Stick_1_lambda_par;
%                     C{index,13} = NaN;
%                 case 'C2'
%                     C{index,10} = json.C2CylinderStejskalTannerApproximation_1_mu(1);
%                     C{index,11} = json.C2CylinderStejskalTannerApproximation_1_mu(2);
%                     C{index,12} = json.C2CylinderStejskalTannerApproximation_1_lambda_par;
%                     C{index,13} = json.C2CylinderStejskalTannerApproximation_1_diameter;
%                 case 'C3'
%                     C{index,10} = json.C3CylinderCallaghanApproximation_1_mu(1);
%                     C{index,11} = json.C3CylinderCallaghanApproximation_1_mu(2);
%                     C{index,12} = json.C3CylinderCallaghanApproximation_1_lambda_par;
%                     C{index,13} = json.C3CylinderCallaghanApproximation_1_diameter;
%                 case 'C4'
%                     C{index,10} = json.C4CylinderGaussianPhaseApproximation_1_mu(1);
%                     C{index,11} = json.C4CylinderGaussianPhaseApproximation_1_mu(2);
%                     C{index,12} = json.C4CCylinderGaussianPhaseApproximation_1_lambda_par;
%                     C{index,13} = json.C4CylinderGaussianPhaseApproximation_1_diameter;
%             end
%             
%             switch extra{1}
%                 case 'G1'
%                     C{index,14} = json.G1Ball_1_lambda_iso;
%                     C{index,15} = NaN;
%                     C{index,16} = NaN;
%                 case 'G2'
%                     C{index,14} = json.G2Zeppelin_1_lambda_par;
%                     C{index,15} = NaN;
%                     C{index,16} = NaN;
%                 case 'G3'
%                     C{index,14} = json.G3TemporalZeppelin_1_lambda_par;
%                     C{index,15} = json.G3TemporalZeppelin_1_lambda_inf;
%                     C{index,16} = json.G3TemporalZeppelin_1_A;
%             end
%             
%            
%             
%         index = index+1;
%         end % rest
%     end % extra
% end % intra

%%
[new,I]=sort(cell2mat(C(:,2)')');
A=[]; A(:,1)=I;
[new,I]=sort(cell2mat(C(:,3)')'); A(:,2)=I;
[new,I]=sort(cell2mat(C(:,4)')');A(:,3)=I;
[new,I]=sort(cell2mat(C(:,5)')');A(:,4)=I;

tmp = cell2mat(C(:,13)')

tmp = cell2mat(C(:,7)')