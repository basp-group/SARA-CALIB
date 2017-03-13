function [param_algo, param_dde, param_im] = Init_param_for_DDE_rec...
    (param_algo, param_dde, param_die, param_im, param_data, results)

% *************************************************************************
% Initialize variables to estimate the DDEs and the image, taking into
% account the estimations obtained during the -initialization process-
% consisting of estimating the DIEs
% see Init_param_for_DIE_init for the initialization of the DIEs
% *************************************************************************
% 
% param_dde: structure containing
%            theta_min: lower bound on DDEs
%            theta_max: upper bound on DDEs
%            U1_init/U2_init: random initialization of DDEs (taking into
%            account the DIEs estimated in the initialization process)
%
% param_im: structure containing
%            min_x: lower bound for the image
%            max_x: upper bound for the image (this bound takes into account
%            the assumption that bright sources are known exactly)
%            im_rec_die: initialization of the image (taking into account
%            the image estimated in the initialization process)
% *************************************************************************


%% ------------------------------------------------------------------------
% parameter for algo dealing with ddes
% -------------------------------------------------------------------------
param_algo.initialization = 0 ; 
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Initialize parameters for DDEs
% -------------------------------------------------------------------------
D_kernel = cell(param_data.T,1) ;
for s = 1:param_data.T
D_kernel{s} = cell2mat(param_dde.U_true{s}) ;
end
dkernel = cell2mat(D_kernel) ;
die_b = dkernel(:,floor(param_dde.Sd^2/2)+1) ;
dde_b = dkernel ;
dde_b(:,floor(param_dde.Sd^2/2)+1) = 0 ;
umax = 1.5 * max(max(abs(real(dde_b(:))),abs(imag(dde_b(:)))))* ones(param_dde.Sd^2,1) ;
umax(floor(param_dde.Sd^2/2)+1) = 1.5 * max(max(abs(real(die_b(:))), abs(imag(die_b(:)))));
umin = - umax ;

param_dde.theta_min = umin ;
param_dde.theta_max = umax ;
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Initialize DDE variables
% -------------------------------------------------------------------------
D0 = cell(param_data.T,1) ;
Dk1 = cell(param_data.T,1) ;
Dk2 = cell(param_data.T,1) ;
c = floor(param_dde.Sd/2) *(param_dde.Sd+1) +1 ;
for s = 1:param_data.T
D0{s} = cell(param_data.na,1) ;
Dk1{s} = cell(param_data.na,1) ;
Dk2{s} = cell(param_data.na,1) ;
for alpha = 1:param_data.na
% initialisation des noyaux
D0{s}{alpha} = (0.3*umax') .* ( 0.5* randn(1, param_dde.Sd^2) + 0.5* 1i * randn(1, param_dde.Sd^2) ) ;
D0{s}{alpha}(c) = D0{s}{alpha}(c) + sign(randn) ;
D0r = max( min( real(D0{s}{alpha}), umax' ), umin' ) ;
D0i = max( min( imag(D0{s}{alpha}), umax' ), umin' ) ;
D0{s}{alpha} = D0r + 1i * D0i ;
% noyaux D1
Dk1{s}{alpha} = D0{s}{alpha} ; 
Dk1{s}{alpha}(floor(param_dde.Sd^2/2)+1) = results.U1_rec{s}{alpha} ; 
% noyaux D2
Dk2{s}{alpha} = D0{s}{alpha} ; 
Dk2{s}{alpha}(floor(param_dde.Sd^2/2)+1) = results.U2_rec{s}{alpha} ; 
end
end

param_dde.U1_init = Dk1 ;
param_dde.U2_init = Dk2 ;
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Initialize image from previous estimation
% -------------------------------------------------------------------------
param_die_tmp = param_die ;
param_die_tmp.U1 = results.U1_rec ;
param_die_tmp.U2 = results.U2_rec ;
G = create_matrix_G_D1_D2_toolbox_test(param_data, param_die_tmp, param_im)  ;
Bstar_die   = @(x)   G * param_im.TF(x);  
x = results.im_rec ;
% Shrink the reconstructed taking into account information in Fourier
% we have Y = YIO + YID + err
YIO = Bstar_die(param_im.xo) ;
YID = Bstar_die(x - param_im.xo) ;
% if YID<<YIO --> do not trust reconstructed image with DIEs
if norm(YID) < param_algo.p_min * norm(YIO)
x_die_T = param_im.xo ;
tau = 1 ;
% else --> we can take information from the reconstructed image
else
im_tmp = x ;
tau = param_algo.Tau ;
while norm(YID) > param_algo.p_max * norm(YIO) && tau < 1
im_tmp(im_tmp<tau*param_im.min_xo) = 0 ;
YID = Bstar_die(im_tmp - param_im.xo) ;
tau = tau + param_algo.Tau ;
end
x_die_T = im_tmp ;
end

param_im.im_rec_die = x_die_T ;
% -------------------------------------------------------------------------

end