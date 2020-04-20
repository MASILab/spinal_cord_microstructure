'''
This Script will cover the 48 Dmipy Model Combinations for Intra-Cellular, Extra-Cellular and Restricted Diffusion
C1: Stick; C2: Stejskal Cylinder; C3: Callaghan Cylinder; C4; GaussPhase Cylinder
G1: Ball; G2: Zeppelin; G3: TemporalZeppelin
S1: Dot; S2: SphereStejskalTannerApproximation; S4: SphereGaussianPhaseApproximation

Orientation Condition that will be applied is that if there are more than two sets of 'mu' then they will be set to equal
as the orientational measures from two different compartments should be equivalent

If there are any parallel diffusivity compartments then they will be set to hard-value of 1.7e-9

Outputs Per Method 2 Files are generated:
1.)The fitting will be reported in the form of the fitted parameters that will be saved as numpy text files.
2.)The second file will contain the name of the fitted parameters
'''

import os
import numpy as np
import csv
import json

# Dmipy Imports
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
from dmipy.signal_models import cylinder_models as cyn  # Intra Diffusion Modeling
from dmipy.signal_models import gaussian_models as gsn  # Extra Diffusion Modeling
from dmipy.signal_models import sphere_models as sph  # Restricted Diffusion Modeling
from dmipy.core.modeling_framework import MultiCompartmentModel


def main():
    # Define Base Data Paths here
    base_data_path = r'/nfs/masi/nathv/spinal_cord_data_2020/SingleVoxelSignals_norm/SingleVoxelSignals_norm'
    base_data_path = os.path.normpath(base_data_path)

    # Define Saving Paths here
    save_data_path = r'/nfs/masi/nathv/spinal_cord_data_2020/results_norm/intra_extra_rest_par_fixed_csf'
    save_data_path = os.path.normpath(save_data_path)
    if os.path.exists(save_data_path) == False:
        os.mkdir(save_data_path)

    # Scheme and The Directions file Paths
    scheme_path = os.path.join(base_data_path, 'scheme.scheme')
    bvecs_path = os.path.join(base_data_path, 'BVECS.bvec')
    bvals_path = os.path.join(base_data_path, 'BVALS.bval')

    # Voxel Paths
    voxel_fc_path = os.path.join(base_data_path, 'FasciulusCuneatus.txt')
    voxel_lc_path = os.path.join(base_data_path, 'LateralCST.txt')
    voxel_sl_path = os.path.join(base_data_path, 'SpinalLemniscus.txt')
    voxel_vc_path = os.path.join(base_data_path, 'VentralCST.txt')
    voxel_vh_path = os.path.join(base_data_path, 'VentralHorn.txt')

    # Reading the Scheme and the Directions
    scheme_data = np.loadtxt(scheme_path)
    bvecs_data = np.loadtxt(bvecs_path)
    bvals_data = np.loadtxt(bvals_path)

    # Read the voxel Data
    fc_data = []
    with open(voxel_fc_path) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=',')
        for row in readcsv:
            fc_data.append(row)
    csvfile.close()
    fc_data = np.asarray(fc_data, dtype='float32')
    print('FC Voxel Shape: {}'.format(fc_data.shape))

    lc_data = []
    with open(voxel_lc_path) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=',')
        for row in readcsv:
            lc_data.append(row)
    csvfile.close()
    lc_data = np.asarray(lc_data, dtype='float32')
    print('LC Voxel Shape: {}'.format(lc_data.shape))

    sl_data = []
    with open(voxel_sl_path) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=',')
        for row in readcsv:
            sl_data.append(row)
    csvfile.close()
    sl_data = np.asarray(sl_data, dtype='float32')
    print('SL Voxel Shape: {}'.format(sl_data.shape))

    vc_data = []
    with open(voxel_vc_path) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=',')
        for row in readcsv:
            vc_data.append(row)
    csvfile.close()
    vc_data = np.asarray(vc_data, dtype='float32')
    print('VC Voxel Shape: {}'.format(vc_data.shape))

    vh_data = []
    with open(voxel_vh_path) as csvfile:
        readcsv = csv.reader(csvfile, delimiter=',')
        for row in readcsv:
            vh_data.append(row)
    csvfile.close()
    vh_data = np.asarray(vh_data, dtype='float32')
    print('VH Voxel Shape: {}'.format(vh_data.shape))

    print('All Data Loaded ...')

    print('Constructing Acquisition Schemes')
    all_bvals = bvals_data * 1e6
    all_bvecs = np.transpose(bvecs_data)

    little_delta = scheme_data[:, 0]
    big_delta = scheme_data[:, 1]
    t_e = scheme_data[:, 4]

    Acq_Scheme = acquisition_scheme_from_bvalues(all_bvals,
                                                 all_bvecs,
                                                 delta=little_delta * 1e-3,
                                                 Delta=big_delta * 1e-3)

    cylinder_dict = {
        'C1': cyn.C1Stick,
        'C2': cyn.C2CylinderStejskalTannerApproximation,
        'C3': cyn.C3CylinderCallaghanApproximation,
        'C4': cyn.C4CylinderGaussianPhaseApproximation
    }

    gaussian_dict = {
        'G1': gsn.G1Ball,
        'G2': gsn.G2Zeppelin
    }

    sphere_dict = {
        'S1': sph.S1Dot,
        'S2': sph.S2SphereStejskalTannerApproximation,
        'S4': sph.S4SphereGaussianPhaseApproximation
    }

    # FC Saving path
    fc_save_path = os.path.join(save_data_path, 'FC')
    if os.path.exists(fc_save_path) == False:
        os.mkdir(fc_save_path)

    lc_save_path = os.path.join(save_data_path, 'LC')
    if os.path.exists(lc_save_path) == False:
        os.mkdir(lc_save_path)

    sl_save_path = os.path.join(save_data_path, 'SL')
    if os.path.exists(sl_save_path) == False:
        os.mkdir(sl_save_path)

    vc_save_path = os.path.join(save_data_path, 'VC')
    if os.path.exists(vc_save_path) == False:
        os.mkdir(vc_save_path)

    vh_save_path = os.path.join(save_data_path, 'VH')
    if os.path.exists(vh_save_path) == False:
        os.mkdir(vh_save_path)

    # TODO Double Combinations of Intra and Extra.
    for cyn_key, cyn_val in cylinder_dict.items():
        for gsn_key, gsn_val in gaussian_dict.items():

            # File name
            model_file_name = cyn_key + '_' + gsn_key + '_par_fix' + '_csf.json'
            signal_file_name = cyn_key + '_' + gsn_key + '_par_fix' + '_csf_signal.txt'

            cylinder = cyn_val()
            gaussian = gsn_val()

            # TODO Defining the compartment of a fixed CSF ball
            csf_ball = gsn.G1Ball()

            multi_compat_model = MultiCompartmentModel(models=[cylinder, gaussian, csf_ball])

            #TODO Always set the CSF Ball the first so the additional iso gets eliminated from the list
            #the edge case is that if the original gaussian is not a ball then the first isotropic does not exist
            #and hence we need an additional if condition to check if the first iso exists or not
            if gsn_key == 'G2':
                multi_compat_model.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)
            elif gsn_key == 'G1':
                multi_compat_model.set_fixed_parameter('G1Ball_2_lambda_iso', 3e-9)

            # TODO If more than two mu exist, implies multiple orientation based measures exist
            # Hence for them we will identify them and set them to be equal to each other.
            mu_list = []
            for each_para_name in multi_compat_model.parameter_names:
                # Last three characters of parameter
                mu_type = each_para_name[-2:]
                if mu_type == 'mu':
                    mu_list.append(each_para_name)

            if len(mu_list) == 2:
                multi_compat_model.set_equal_parameter(mu_list[0], mu_list[1])
            # End of mu conditions

            # TODO We will be setting all parallel diffusivity parameters to 1.7e-9
            for each_para_name in multi_compat_model.parameter_names:
                # Last three characters of parameter
                diffusivity_type = each_para_name[-3:]
                if diffusivity_type == 'par':
                    multi_compat_model.set_fixed_parameter(each_para_name, 1.7e-9)
            # End of Parallel Diffusivity Parameters

            # TODO We will be setting Isotropic diffusivity parameters to 1.7e-9
            for each_para_name in multi_compat_model.parameter_names:
                # Last three characters of parameter
                diffusivity_type = each_para_name[-3:]
                if diffusivity_type == 'iso':
                    multi_compat_model.set_fixed_parameter(each_para_name, 1.7e-9)
            # End of Isotropic Diffusivity Parameters

            print(multi_compat_model.parameter_names)

            ######## FC #########
            fc_model_fit = multi_compat_model.fit(Acq_Scheme, fc_data, use_parallel_processing=False, solver='mix')
            fc_fitted_params = fc_model_fit.fitted_parameters
            fc_model_signal = fc_model_fit.predict()

            ## Save FC Signal
            fc_model_signal = fc_model_signal[0, :].tolist()
            signal_save_path = os.path.join(fc_save_path, signal_file_name)
            np.savetxt(signal_save_path, fc_model_signal)
            #################

            ## Error Calculations
            fc_mse = fc_model_fit.mean_squared_error(fc_data)
            ##

            new_params = {}
            new_params['mse'] = fc_mse.tolist()

            for key, value in fc_fitted_params.items():
                new_params[key] = value.tolist()

            model_save_path = os.path.join(fc_save_path, model_file_name)

            with open(model_save_path, 'w') as json_file:
                json.dump(new_params, json_file)
            json_file.close()
            #####################

            ######## LC #########
            lc_model_fit = multi_compat_model.fit(Acq_Scheme, lc_data, use_parallel_processing=False, solver='mix')
            lc_fitted_params = lc_model_fit.fitted_parameters
            lc_model_signal = lc_model_fit.predict()

            ## Save LC Signal
            lc_model_signal = lc_model_signal[0, :].tolist()
            signal_save_path = os.path.join(lc_save_path, signal_file_name)
            np.savetxt(signal_save_path, lc_model_signal)
            #################

            ## Error Calculations
            lc_mse = lc_model_fit.mean_squared_error(lc_data)
            ##

            new_params = {}
            new_params['mse'] = lc_mse.tolist()

            for key, value in lc_fitted_params.items():
                new_params[key] = value.tolist()

            model_save_path = os.path.join(lc_save_path, model_file_name)

            with open(model_save_path, 'w') as json_file:
                json.dump(new_params, json_file)
            json_file.close()
            #####################

            ######## SL #########
            sl_model_fit = multi_compat_model.fit(Acq_Scheme, sl_data, use_parallel_processing=False, solver='mix')
            sl_fitted_params = sl_model_fit.fitted_parameters
            sl_model_signal = sl_model_fit.predict()

            ## Save SL Signal
            sl_model_signal = sl_model_signal[0, :].tolist()
            signal_save_path = os.path.join(sl_save_path, signal_file_name)
            np.savetxt(signal_save_path, sl_model_signal)
            #################

            ## Error Calculations
            sl_mse = sl_model_fit.mean_squared_error(sl_data)
            ##

            new_params = {}
            new_params['mse'] = sl_mse.tolist()

            for key, value in sl_fitted_params.items():
                new_params[key] = value.tolist()

            model_save_path = os.path.join(sl_save_path, model_file_name)

            with open(model_save_path, 'w') as json_file:
                json.dump(new_params, json_file)
            json_file.close()
            ######################

            ######## VC ##########
            vc_model_fit = multi_compat_model.fit(Acq_Scheme, vc_data, use_parallel_processing=False, solver='mix')
            vc_fitted_params = vc_model_fit.fitted_parameters
            vc_model_signal = vc_model_fit.predict()

            ## Save VC Signal
            vc_model_signal = vc_model_signal[0, :].tolist()
            signal_save_path = os.path.join(vc_save_path, signal_file_name)
            np.savetxt(signal_save_path, vc_model_signal)
            #################

            ## Error Calculations
            vc_mse = vc_model_fit.mean_squared_error(vc_data)
            ##

            new_params = {}
            new_params['mse'] = vc_mse.tolist()

            for key, value in vc_fitted_params.items():
                new_params[key] = value.tolist()

            model_save_path = os.path.join(vc_save_path, model_file_name)

            with open(model_save_path, 'w') as json_file:
                json.dump(new_params, json_file)
            json_file.close()
            ######################

            ######## VH #########
            vh_model_fit = multi_compat_model.fit(Acq_Scheme, vh_data, use_parallel_processing=False, solver='mix')
            vh_fitted_params = vh_model_fit.fitted_parameters
            vh_model_signal = vh_model_fit.predict()

            ## Save VH Signal
            vh_model_signal = vh_model_signal[0, :].tolist()
            signal_save_path = os.path.join(vh_save_path, signal_file_name)
            np.savetxt(signal_save_path, vh_model_signal)
            #################

            ## Error Calculations
            vh_mse = vh_model_fit.mean_squared_error(vh_data)
            ##

            new_params = {}
            new_params['mse'] = vh_mse.tolist()

            for key, value in vh_fitted_params.items():
                new_params[key] = value.tolist()

            model_save_path = os.path.join(vh_save_path, model_file_name)

            with open(model_save_path, 'w') as json_file:
                json.dump(new_params, json_file)
            json_file.close()
            ######################

            print('Model Completed wit Combination of {} and {}'.format(cyn_key, gsn_key))

    # TODO Triple Combinations of Intra, Extra and Rest
    for cyn_key, cyn_val in cylinder_dict.items():
        for gsn_key, gsn_val in gaussian_dict.items():
            for sph_key, sph_val in sphere_dict.items():

                cylinder = cyn_val()
                gaussian = gsn_val()
                sphere = sph_val()

                # TODO Defining the compartment of a fixed CSF ball
                csf_ball = gsn.G1Ball()

                multi_compat_model = MultiCompartmentModel(models=[cylinder, gaussian, sphere, csf_ball])

                print(multi_compat_model.parameter_names)
                # TODO Always set the CSF Ball the first so the additional iso gets eliminated from the list
                # the edge case is that if the original gaussian is not a ball then the first isotropic does not exist
                # and hence we need an additional if condition to check if the first iso exists or not
                if gsn_key == 'G2':
                    multi_compat_model.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)
                elif gsn_key == 'G1':
                    multi_compat_model.set_fixed_parameter('G1Ball_2_lambda_iso', 3e-9)

                #print(multi_compat_model.parameter_names)

                # TODO If more than two mu exist, implies multiple orientation based measures exist
                # Hence for them we will identify them and set them to be equal to each other.
                mu_list = []
                for each_para_name in multi_compat_model.parameter_names:
                    # Last three characters of parameter
                    mu_type = each_para_name[-2:]
                    if mu_type == 'mu':
                        mu_list.append(each_para_name)
                        # multi_compat_model.set_fixed_parameter(each_para_name, 1.7e-9)

                if len(mu_list) == 2:
                    multi_compat_model.set_equal_parameter(mu_list[0], mu_list[1])
                # End of mu conditions

                # TODO We will be setting all parallel diffusivity parameters to 1.7e-9
                for each_para_name in multi_compat_model.parameter_names:
                    # Last three characters of parameter
                    diffusivity_type = each_para_name[-3:]
                    if diffusivity_type == 'par':
                        multi_compat_model.set_fixed_parameter(each_para_name, 1.7e-9)
                # End of Parallel Diffusivity Parameters

                # TODO We will be setting Isotropic diffusivity parameters to 1.7e-9
                for each_para_name in multi_compat_model.parameter_names:
                    # Last three characters of parameter
                    diffusivity_type = each_para_name[-3:]
                    if diffusivity_type == 'iso':
                        multi_compat_model.set_fixed_parameter(each_para_name, 1.7e-9)
                # End of Isotropic Diffusivity Parameters

                # This file name is common to all voxels and describes the nomenclature
                # as the selection of models that were used based on the three components
                model_file_name = cyn_key + '_' + gsn_key + '_' + sph_key + '_par_fix' + '_csf.json'
                signal_file_name = cyn_key + '_' + gsn_key + '_' + sph_key + '_par_fix' + '_csf_signal.txt'

                ######## FC #########
                fc_model_fit = multi_compat_model.fit(Acq_Scheme, fc_data, use_parallel_processing=False, solver='mix')
                fc_fitted_params = fc_model_fit.fitted_parameters
                fc_model_signal = fc_model_fit.predict()

                ## Save FC Signal
                fc_model_signal = fc_model_signal[0, :].tolist()
                signal_save_path = os.path.join(fc_save_path, signal_file_name)
                np.savetxt(signal_save_path, fc_model_signal)
                #################

                ## Error Calculations
                fc_mse = fc_model_fit.mean_squared_error(fc_data)
                ##

                new_params = {}
                new_params['mse'] = fc_mse.tolist()

                for key, value in fc_fitted_params.items():
                    new_params[key] = value.tolist()

                model_save_path = os.path.join(fc_save_path, model_file_name)

                with open(model_save_path, 'w') as json_file:
                    json.dump(new_params, json_file)
                json_file.close()
                #####################

                ######## LC #########
                lc_model_fit = multi_compat_model.fit(Acq_Scheme, lc_data, use_parallel_processing=False, solver='mix')
                lc_fitted_params = lc_model_fit.fitted_parameters
                lc_model_signal = lc_model_fit.predict()

                ## Save LC Signal
                lc_model_signal = lc_model_signal[0, :].tolist()
                signal_save_path = os.path.join(lc_save_path, signal_file_name)
                np.savetxt(signal_save_path, lc_model_signal)
                #################

                ## Error Calculations
                lc_mse = lc_model_fit.mean_squared_error(lc_data)
                ##

                new_params = {}
                new_params['mse'] = lc_mse.tolist()

                for key, value in lc_fitted_params.items():
                    new_params[key] = value.tolist()

                model_save_path = os.path.join(lc_save_path, model_file_name)

                with open(model_save_path, 'w') as json_file:
                    json.dump(new_params, json_file)
                json_file.close()
                #####################

                ######## SL #########
                sl_model_fit = multi_compat_model.fit(Acq_Scheme, sl_data, use_parallel_processing=False, solver='mix')
                sl_fitted_params = sl_model_fit.fitted_parameters
                sl_model_signal = sl_model_fit.predict()

                ## Save SL Signal
                sl_model_signal = sl_model_signal[0, :].tolist()
                signal_save_path = os.path.join(sl_save_path, signal_file_name)
                np.savetxt(signal_save_path, sl_model_signal)
                #################

                ## Error Calculations
                sl_mse = sl_model_fit.mean_squared_error(sl_data)
                ##

                new_params = {}
                new_params['mse'] = sl_mse.tolist()

                for key, value in sl_fitted_params.items():
                    new_params[key] = value.tolist()

                model_save_path = os.path.join(sl_save_path, model_file_name)

                with open(model_save_path, 'w') as json_file:
                    json.dump(new_params, json_file)
                json_file.close()
                ######################

                ######## VC ##########
                vc_model_fit = multi_compat_model.fit(Acq_Scheme, vc_data, use_parallel_processing=False, solver='mix')
                vc_fitted_params = vc_model_fit.fitted_parameters
                vc_model_signal = vc_model_fit.predict()

                ## Save VC Signal
                vc_model_signal = vc_model_signal[0, :].tolist()
                signal_save_path = os.path.join(vc_save_path, signal_file_name)
                np.savetxt(signal_save_path, vc_model_signal)
                #################

                ## Error Calculations
                vc_mse = vc_model_fit.mean_squared_error(vc_data)
                ##

                new_params = {}
                new_params['mse'] = vc_mse.tolist()

                for key, value in vc_fitted_params.items():
                    new_params[key] = value.tolist()

                model_save_path = os.path.join(vc_save_path, model_file_name)

                with open(model_save_path, 'w') as json_file:
                    json.dump(new_params, json_file)
                json_file.close()
                ######################

                ######## VH #########
                vh_model_fit = multi_compat_model.fit(Acq_Scheme, vh_data, use_parallel_processing=False, solver='mix')
                vh_fitted_params = vh_model_fit.fitted_parameters
                vh_model_signal = vh_model_fit.predict()

                ## Save VH Signal
                vh_model_signal = vh_model_signal[0, :].tolist()
                signal_save_path = os.path.join(vh_save_path, signal_file_name)
                np.savetxt(signal_save_path, vh_model_signal)
                #################

                ## Error Calculations
                vh_mse = vh_model_fit.mean_squared_error(vh_data)
                ##

                new_params = {}
                new_params['mse'] = vh_mse.tolist()

                for key, value in vh_fitted_params.items():
                    new_params[key] = value.tolist()

                model_save_path = os.path.join(vh_save_path, model_file_name)

                with open(model_save_path, 'w') as json_file:
                    json.dump(new_params, json_file)
                json_file.close()
                ######################

                print('Model Completed with Combination of {} and {} and {}'.format(cyn_key, gsn_key, sph_key))

    print('All Done')


if __name__ == "__main__":
    main()
