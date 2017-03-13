function [param_algo, param_die, param_im] = Init_param_for_DIE_init...
    (param_algo, param_die, param_im, param_data)


% *************************************************************************
% Initialize variables to estimate the DIEs and have a first approximation
% of the image, starting from known brightest sources.
% *************************************************************************
% 
% param_die: structure containing
%            theta_min: lower bound on DIEs
%            theta_max: upper bound on DIEs
%            U1_init/U2_init: random initialization of DIEs
%
% param_im: structure containing
%            min_x: lower bound for the image
%            max_x: upper bound for the image (this bound takes into account
%            the assumption that bright sources are known exactly)
% *************************************************************************


%% ------------------------------------------------------------------------
% parameter for algo dealing with dies
% -------------------------------------------------------------------------
param_algo.initialization = 1 ; 
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Initialize parameters for DIEs
% -------------------------------------------------------------------------
D_kernel = cell(param_data.T,1) ;
parfor s = 1:param_data.T
D_kernel{s} = cell2mat(param_die.U_true{s}) ;
end
dkernel = cell2mat(D_kernel) ;
die_b = dkernel(:,floor(param_die.Sd^2/2)+1) ;
dde_b = dkernel ;
dde_b(:,floor(param_die.Sd^2/2)+1) = 0 ;
umax = 1.5 * max(max(abs(real(dde_b(:))),abs(imag(dde_b(:)))))* ones(param_die.Sd^2,1) ;
umax(floor(param_die.Sd^2/2)+1) = 1.5 * max(max(abs(real(die_b(:))), abs(imag(die_b(:)))));
umin = - umax ;

param_die.theta_min = umin ;
param_die.theta_max = umax ;
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Initialize parameters for image
% -------------------------------------------------------------------------
param_im.min_x = zeros(size(param_im.xo)) ;
param_im.min_x(param_im.xo>0) = 0 ;
param_im.max_x = 0 * ones(size(param_im.min_x(param_im.xo>0))) ;
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Initialize DIE variables
% -------------------------------------------------------------------------
D0 = cell(param_data.T,1) ;
Dk1 = cell(param_data.T,1) ;
Dk2 = cell(param_data.T,1) ;
for s = 1:param_data.T
D0{s} = cell(param_data.na,1) ;
Dk1{s} = cell(param_data.na,1) ;
Dk2{s} = cell(param_data.na,1) ;
for alpha = 1:param_data.na
% initialisation des noyaux
D0{s}{alpha} = ...
    sign(randn(1,param_die.Sd^2)) + 0.5*0.3* randn(1,param_die.Sd^2) ...
    +1i * sign(randn(1,param_die.Sd^2)) + 0.5*0.3* randn(1,param_die.Sd^2) ;
D0r = max( min( real(D0{s}{alpha}), umax' ), umin' ) ;
D0i = max( min( imag(D0{s}{alpha}), umax' ), umin' ) ;
D0{s}{alpha} = D0r + 1i * D0i ;
% noyaux D1
Dk1{s}{alpha} = D0{s}{alpha} ; 
% noyaux D2
Dk2{s}{alpha} = D0{s}{alpha} ; 
end
end

param_die.U1_init = Dk1 ;
param_die.U2_init = Dk2 ;
% -------------------------------------------------------------------------


end