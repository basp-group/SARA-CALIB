% *************************************************************************
% Matlab code for joint imaging and DDE calibration in RI
% *************************************************************************
% Reconstruction considering image of M31 and band limited DDE
% 
% Generate random antenna distribution
% DDE generated randomly in Fourier
% 
% Related paper:
% A non-convex optimization algorithm for joint DDE calibration and imaging
% in radio interferometry
% Audrey Repetti, Jasleen Birdy, Arwa Dabbech, Yves Wiaux
% Arxiv preprint: 1701.03689
%
% Contact: a.repetti@hw.ac.uk
% *************************************************************************

close all
clc
clear all

addpath('./Tools')
addpath('./Create_data')



%% Generate inverse problem

% ========================================== %
Generate_parameters_M31
% ========================================== %

[param_data, param_dde, param_die, param_im] = Generate_random_RI_data...
    (param_data, param_dde, param_die, param_im) ;



%% Initialization: Joint reconstruction of the image and the DIEs

disp(' ')
disp('==========================================')
disp('DIEs estimation...')
disp('==========================================')

% Initialization for DIEs
[param_algo, param_die, param_im] = Init_param_for_DIE_init...
    (param_algo, param_die, param_im, param_data) ;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
results = Joint_imaging_dde_calibration_reconstruction...
    (param_im, param_die, param_algo, param_data) ;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




%% Main algorithm: Joint reconstruction of the image and the DDEs

disp(' ')
disp('==========================================')
disp('DDEs estimation...')
disp('==========================================')

% Initialization for DDEs
[param_algo, param_dde, param_im] = Init_param_for_DDE_rec...
    (param_algo, param_dde, param_die, param_im, param_data, results) ;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
results_dde = Joint_imaging_dde_calibration_reconstruction...
    (param_im, param_dde, param_algo, param_data) ;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





