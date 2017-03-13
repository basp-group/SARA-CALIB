function param_im = create_M31_image(param_im)

SNR =@(x, xtrue) 20 * log10( sqrt( sum( xtrue(:).^2 ) / sum( (x(:)-xtrue(:)).^2 ) ) ) ;


%% Load image

im=fitsread('M31.fits'); 
im = imresize(im,param_im.Ni);
im=(im+abs(im))./2;
im=im./max(im(:));



%% Generate approximated image from the true image


% energy of the images image
E_im_2 = sum(im(:).^2) ;
E_im = sqrt(E_im_2) ;

E_im_app = sqrt(1/(1+param_im.kappa^2)) * E_im ;
E_im_b = param_im.kappa * E_im_app ;


% generate approximation
tau_min = 0 ; tau_max = 1 ;
E_max = E_im ; E_min = 0 ;
diff_sup = E_im_app - E_max ;
diff_inf = E_im_app - E_min ;
if diff_sup > diff_inf
tau = tau_max ; xo = 0*im ;
else
tau = tau_min ; xo = im ;
end

cond_stop = 1 ;
k = 0 ;
k_max = 100 ;
tol = 1e-4 ;

disp(['Find approximated image with theoretical energy = ',num2str(E_im_app)])
disp(['Global energy of the image : ',num2str(E_im)])
disp(['Theoretical energy of the background image : ',num2str(E_im_b)])
disp(' ')

while cond_stop>tol && k<k_max
k = k+1 ;

tau_tmp = (tau_min + tau_max)/2 ;

im_app_tmp = im ;
im_app_tmp(im_app_tmp<tau_tmp) = 0 ;

E_tmp = sqrt(sum(im_app_tmp(:).^2)) ;


if E_im_app > E_tmp
tau_max = tau_tmp ;
E_min = E_tmp ;
cond = E_im_app - E_min ;
diff_inf = cond ;
if cond < diff_sup
xo = im_app_tmp ;
tau = tau_max ;
end
else
tau_min = tau_tmp ;
E_max = E_tmp ;
cond = E_max - E_im_app ;
diff_sup = cond ;
if cond < diff_inf
xo = im_app_tmp ;
tau = tau_min ;
end
end

cond_stop = abs(E_im_app - E_tmp) ;

end


eps_true = im - xo ;

E1 = sqrt( sum(xo(:).^2) ) ;
E2 = sqrt( sum(eps_true(:).^2) ) ;
disp('********************************************')
disp(['Threshold parameter tau = ',num2str(tau)])
disp(['Energy of the approximated image = ',num2str(sqrt( sum(xo(:).^2) ))])
disp(['Energy of the background image = ',num2str(sqrt( sum(eps_true(:).^2) ))])
disp(['SNR approximated image = ',num2str(SNR(xo,im))])
disp('********************************************')



param_im.tau = tau ; % save threshold operator
param_im.E1 = E1 ;
param_im.E2 = E2 ;



%% Truncate image to have an odd image size
if mod(param_im.Ni(1),2)==0 % probleme de symetrie convolution
im = im(1:end-1,1:end-1) ;
xo = xo(1:end-1,1:end-1) ;
end

%% Save image
param_im.im_true = im ;
param_im.xo = xo ;
param_im.min_xo = min(param_im.xo(:)>0) ;
param_im.eps_true = eps_true ;

param_im.min = 0 ;
param_im.max = Inf ;
param_im.mask_xo = (xo>0) ;


end