function [param_im] = create_random_image(param_im)

% for the convolution to have Gaussian sources instead of just points
size_ker_h = [3,3] ; 
var_ker_h = 0.5 ;
h = fspecial('gaussian',size_ker_h,var_ker_h);


%% Generate positions


c = floor(param_im.Ni(1)/2)+1 ;
c_inf = c- floor(param_im.Ni(1)/(3)) ;
c_sup = c+ floor(param_im.Ni(1)/(3)) ;
min_ind = 3 ;
ind_tot = [] ;

for l = 1:param_im.n_lev

if l <=2 
window_sources = [c_inf c_sup] ;
else
window_sources = [ 2+min_ind , param_im.Ni(1)-min_ind ] ;
end

ind = randi(window_sources, param_im.I_lev(l),2) ;
ind = unique(ind,'rows') ;
if l>1
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
end
size_ind = max(size(ind)) ;
while size_ind<param_im.I_lev(l)
ind_tmp = randi(window_sources, param_im.I_lev(l)-size_ind,2) ;
ind = unique([ind ; ind_tmp],'rows') ;
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
size_ind = max(size(ind)) ;
end


if l == 1
ind_lev{l} = ind ;
x_tmp = zeros(param_im.Ni) ;
for i =1:param_im.I_lev(l)
x_tmp(ind(i,1), ind(i,2)) = 1 ;
end
h_tmp = fspecial('gaussian',size_ker_h+2,var_ker_h);
x_tmp = conv2(x_tmp,h_tmp,'same') ;
x_tmp = x_tmp(:) ;
ind_tmp = find(x_tmp>0) ;
ind_tot = convert_1D_2D( ind_tmp, param_im.Ni(1)) + (param_im.Ni(1)/2+1) ;
ind_tot = [ind_tot(:,2), ind_tot(:,1)] ;
else
ind_lev{l} = ind ;
ind_tot = [ind_tot ; ind] ;
end
    

end



%% Generate amplitudes


mean_lev = zeros(param_im.n_lev,1) ; 
std_lev = zeros(param_im.n_lev,1) ;
amp_lev = cell(1,param_im.n_lev) ;

for l = 1:param_im.n_lev
mean_lev(l) = param_im.En_lev(l) / sqrt( param_im.I_lev(l) * (param_im.rho^2+1) ) ;
std_lev(l) = param_im.rho * mean_lev(l) ;

amp_lev{l}(1:param_im.I_lev(l)) = std_lev(l) * randn(param_im.I_lev(l),1) + mean_lev(l) ;
end



%% Create images


im = zeros(param_im.Ni) ;
x_level = cell(length(param_im.n_lev),1) ;

for l = 1:param_im.n_lev

x_level{l} = zeros(param_im.Ni) ;
for i = 1:param_im.I_lev(l)
x_level{l}(ind_lev{l}(i,1), ind_lev{l}(i,2)) = amp_lev{l}(i) ;
end
x_level{l} = conv2(x_level{l},h,'same') ;
energy_temp = sqrt(sum(x_level{l}(:).^2)) ;
alpha = param_im.En_lev(l) / energy_temp ;
x_level{l} = alpha*x_level{l} ;
im = im+x_level{l} ;
end

% save the different levels
param_im.x_level = x_level ;
% save the source positions
param_im.ind_lev = ind_lev ;

% define image containing bright sources
xo = x_level{1} ;

if mod(param_im.Ni(1),2)==0 % to have an odd image
im = im(1:end-1,1:end-1) ;
xo = xo(1:end-1,1:end-1) ;
end

param_im.im_true = im ;
param_im.x = im ;
param_im.xo = xo ;
param_im.min_xo = min(param_im.xo(:)>0) ;
param_im.eps_true = param_im.im_true-param_im.xo ;

param_im.min = 0 ;
param_im.max = Inf ;
param_im.mask_xo = (xo>0) ;

end
