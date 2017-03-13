function [param_data, param_dde, param_die, param_im] = Generate_random_RI_data...
    (param_data, param_dde, param_die, param_im)

% *************************************************************************
% Generate_random_RI_data for joint calibration & imaging
%
% Main output variables:
%
% param_im : structure containing
%            im_true: original simulated image
%            eps_true: original unknown sources of im_true
%            xo: known sources of im_true
%
% param_dde/param_die : structure containing
%            U_true: original simulated Fourier kernels of ddes/dies
%            
% param_data : structure containing
%            ud/vd: discrete u-v antenna positions
%            Uabd: discrete corresponding u-v coverage
%            y: visibilities as a vector
%            Y: vibilities as a cube
% *************************************************************************

%% Generate the image

if strcmp(param_im.im_choice,'rand_im')
param_im = create_random_image(param_im) ;
elseif strcmp(param_im.im_choice,'M31')
param_im = create_M31_image(param_im) ;
end



%% U-V coverage + antenna positions


[param_data, param_dde, param_die] = ...
    Generate_discrete_uv_coverage ...
    (param_data, param_dde, param_die, param_im) ;



%% Generate DDEs and useful variables for convolutions


[param_dde, param_die, param_data] = ...
            Generate_ddes(param_dde, param_die, param_data, param_im) ;



%% Generate visibilities


% Construire visibilites vecteur a partir de Y
y_ind_mask = zeros(param_data.na) ;
for alpha = 1:param_data.na-1
for beta = alpha+1:param_data.na
y_ind_mask(alpha, beta) = 1 ;
end
end
y_ind = find(transpose(y_ind_mask)) ;



param_dde.U1 = param_dde.U_true ;
param_dde.U2 = param_dde.U_true ;
G = create_matrix_G_D1_D2_toolbox(param_data, param_dde, param_im) ;
B_   = @(x)   G * param_im.TF(x);  

y0=B_(param_im.im_true);
sigma_noise = sum(abs(y0).^2)/ size(y0,1) * 10^(-param_data.input_snr/20) ;
noise = (randn(size(y0)) + 1i*randn(size(y0)))*sigma_noise/sqrt(2);
% noisy data vector -------------------------------------------------------
y=y0 +noise;

Y = cell(param_data.T,1) ; % noisy data cube
for s = 1:param_data.T
Ytmp = zeros(param_data.na^2,1) ;
Ytmp(y_ind) = y((s-1)*param_data.ant_pair+1 : s*param_data.ant_pair) ;
Y{s} = reshape(conj(Ytmp),param_data.na,param_data.na) + transpose(reshape(Ytmp,param_data.na,param_data.na)) ;
end




param_data.Y = Y ;
param_data.y = y ;










end