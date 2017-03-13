% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Definition of the different parameters 
% -------------------------------------------------------------------------
% param_dde   : parameters for the DDEs
% param_die   : parameters for the DIEs
%               structure containing
%               Sd: dimension of DDE/DIE = Sd x Sd
%               var: std. dev. to generate DDE/DIE
%               [...]
%
% param_im    : parameters for the image
%               structure containing
%               Ni: dimension of the image
%               Nt: dimension in Fourier
%               TF/TF_adj: discrete 2D Fourier transform
%               [...]
%
% param_data  : parameters for the visibilities
%               structure containing
%               na: number of antennas
%               T: measurements/antenna pair
%               ant_pair: number of antenna pairs
%               M: number of total measurements (dimension of visibility
%               vector)
%               [...]
%
% param_algo  : parameters for the reconstruction algorithm
%               structure containing
%               Psi/Psit: Sparsity operator
%               [...]
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% DDES parameters
% -------------------------------------------------------------------------
param_dde.Sd = 5 ; % choose an odd number, DDE of size (Sd x Sd)
param_die.Sd = 1 ;
param_dde.var = 0.05 ;
param_die.var = param_dde.var ;
% ------------------------



% -------------------------------------------------------------------------
% Image parameters
% -------------------------------------------------------------------------
param_im.im_choice = 'rand_im' ;  % random image
param_im.Ni = 128 * [1,1] ;
% -----------------------------------------
% If random image
if strcmp(param_im.im_choice,'rand_im')
param_im.n_lev = 3 ; % number of intensity levels
param_im.I_lev = [10, 10, 10] ;  % number of sources for each level
param_im.En_lev = [10, 0.1, 1e-5] ; % energy of each level
param_im.rho = 0.1 ; % percentage for std deviation: std_lev = rho * moy_lev
end
% -----------------------------------------
param_im.min = 0 ;
param_im.max = +Inf ;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% FFT definition 
% -------------------------------------------------------------------------
[param_im.TF, param_im.TF_adj,param_im.Nt,param_im.Nadd] = ...
TF_operators2(param_dde.Sd,param_im.Ni) ; 
param_im.Index_tshift = reshape(1:prod(param_im.Nt), param_im.Nt);
param_im.Index_tshift = fftshift(param_im.Index_tshift) ;
param_im.Index_tshift = param_im.Index_tshift(:) ;
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Coverage parameters
% -------------------------------------------------------------------------
param_data.T = 1 ;
param_data.na = 200 ;
param_data.M = param_data.T*param_data.na*(param_data.na-1)/2 ;
param_data.ant_pair = param_data.na*(param_data.na-1)/2 ;

disp('--------------------------------------------')
disp(['nombre de mesures : ',num2str(param_data.M)])
disp(['nombre d antennes : ',num2str(param_data.na)])
disp(['nombre de paires d antennes : ',num2str(param_data.ant_pair)])
disp(['nombre de pts / paires d antennes : ',num2str(param_data.T)])
disp('--------------------------------------------')
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Noise level noisy measurements
% -------------------------------------------------------------------------
param_data.input_snr = 150 ; % no noise considered 
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% regularisation parameters
% -------------------------------------------------------------------------
param_algo.eta = 1e-4 ; % l1 image
param_algo.Jeps = 1000 ;
param_algo.tol_x = 1e-5 ;
param_algo.tol_crit = 2e-2 ;
% parameters for DIE estimation
param_die.max_it =  10 ;
param_die.JUtot = 100 ;
param_die.JU1 = 2 ;
param_die.JU2 = 2 ;
param_die.nu = 2000 ;
param_die.tol_norm = 1e-6 ;
% parameters for DDE estimation
param_dde.max_it =  100 ;
param_dde.JUtot = 10 ;
param_dde.JU1 = 5 ;
param_dde.JU2 = 5 ;
param_dde.nu = 1000 ;
param_dde.tol_norm = 1e-5 ;
% -------------------------------------------------------------------------
% percentage for the threshold of the image after reconstruction with DIEs
% -------------------------------------------------------------------------
param_algo.Tau = 0.01 ; 
param_algo.p_min = 0.6 ; 
param_algo.p_max = 0.9 ; 
% -------------------------------------------------------------------------
% Regularisation 
% -------------------------------------------------------------------------
param_algo.nlevel = 4 ;
param_algo.opt_dict = 'Dirac' ;
[param_algo.Psi,param_algo.Psit] = SARA_sparse_operator(randn(param_im.Ni-1), param_algo.nlevel,param_algo.opt_dict) ;
% -------------------------------------------------------------------------
