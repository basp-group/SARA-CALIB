function [param_data, param_dde, param_die] = ...
    Generate_discrete_uv_coverage(param_data, param_dde, param_die, param_im)



% ----------------------------------------------------- %
% random cont. antenna positions                        %
% ----------------------------------------------------- %
uv(1,:) = rand(1,2) ;                                   %
for alpha = 2:param_data.na                             %
uv_ = rand(1,2) ;                                       %
while ismember(uv_,uv,'rows')                           %
uv_ = rand(1,2) ;                                       %
end                                                     %
uv(alpha,:) = uv_ ;                                     %
end                                                     %
Pos_ant = 1e06*[uv, zeros(param_data.na,1)] ;           %
% ----------------------------------------------------- %



% -------------------------------------------------------------------------
h = linspace(-5, 5, param_data.T)*pi/12;% hour angle range of +/-  hours
dec = (pi/180)*(40);  % Cas A is 56.4
lat = (pi/180)*(38. + 59./60. + 48./3600.);  % College Park
% position de reference x0 dans le plan x-y
x0 = [mean(Pos_ant(:,1)), mean(Pos_ant(:,2))] ;
% -------------------------------------------------------------------------
[u,v,w] = generate_uv_cov_antennas(Pos_ant,x0,h,lat,dec, param_data.T) ;



% -------------------------------------------------------------------------
% discrete antenna positions
% -------------------------------------------------------------------------
bmax=4*max(sqrt(cell2mat(u).^2+cell2mat(v).^2));
du=floor(bmax/(param_im.Ni(1)));
bmax=du*(param_im.Ni(1)-1) ;
% -------------------------------------------
ud = cell(param_data.na,1) ; % useful for antenna positions
vd = cell(param_data.na,1) ;
for alpha = 1:param_data.na
ud{alpha} = floor(u{alpha}./du).*du ;
vd{alpha} = floor(v{alpha}./du).*du ;

ud{alpha} = ud{alpha}/bmax + 1/2 ;
vd{alpha} = vd{alpha}/bmax + 1/2 ;
ud{alpha} = ceil(ud{alpha}*(param_im.Ni(1)-1))+1 ;
vd{alpha} = ceil(vd{alpha}*(param_im.Ni(1)-1))+1 ;
end
param_data.ud = ud ;
param_data.vd = vd ;
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Discrete uv-cov
% -------------------------------------------------------------------------
[Uab_d, Vab_d, Ant] = uaub2uab_disc(param_data.na,ud,vd,param_data.T,param_im.Ni(1)) ;
param_data.Uabd = cell2mat(Uab_d) ; param_data.Vabd = cell2mat(Vab_d) ; % useful to create Fourier selection
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% antennes associees au uv-coverage
% -------------------------------------------------------------------------
AntF = [] ;
for s = 1:param_data.T
ant = [Ant{s}{1} Ant{s}{2}] ;
AntF = [AntF ; ant ] ;
end
if param_data.M ~= max(size(AntF))
disp('**************** probleme ****************')
end
clear AntB
for i = 1:param_data.M
AntB(i,:)=[AntF(param_data.M-i+1,2) ; AntF(param_data.M-i+1,1)] ;
end
param_data.antennas=[AntF;AntB];
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
% Support ddes for convolution
% -------------------------------------------------------------------------
Support_Dalpha = cell(param_data.T,1) ;
Support_Dalpha_ = cell(param_data.T,1) ;
Support_tot = [] ;
for s = 1:param_data.T
Support_Dalpha{s} = cell(param_data.na,1) ;
Support_Dalpha_{s} = cell(param_data.na,1) ;
for alpha = 1:param_data.na
temp = zeros(param_im.Ni) ;
temp(vd{alpha}(s)-floor(param_dde.Sd/2)+1:vd{alpha}(s)+floor(param_dde.Sd/2)+1, ...
    ud{alpha}(s)-floor(param_dde.Sd/2)+1:ud{alpha}(s)+floor(param_dde.Sd/2)+1) = 1 ;
temp_z = zeros(param_im.Nt) ;
temp_z(param_im.Nadd/2+1:end-param_im.Nadd/2, param_im.Nadd/2+1:end-param_im.Nadd/2) = temp ;
temp2 = temp_z(2:end, 2:end) ;
temp2 = flipud(fliplr(temp2)) ;
temp2_z = zeros(param_im.Nt) ;
temp2_z(2:end, 2:end) = temp2 ;
%-----------------
temp = temp_z(:) ;
temp2 = temp2_z(:) ;
Support_Dalpha{s}{alpha} = find(temp) ;
Support_Dalpha_{s}{alpha} = find(temp2) ;
Support_tot = [Support_tot ; Support_Dalpha{s}{alpha} ; Support_Dalpha_{s}{alpha}] ;
Support_tot = unique(Support_tot) ;
end
end
param_dde.Support_Dalpha1 = Support_Dalpha ;
param_dde.Support_Dalpha2 = Support_Dalpha_ ;
param_dde.Support_tot = Support_tot ;
param_die.Support_tot = Support_tot ;

Support_Dalpha1_tmp = cell(param_data.T,1) ;
Support_Dalpha2_tmp = cell(param_data.T,1) ;
parfor s = 1:param_data.T
for alpha = 1:param_data.na
Support_Dalpha1_tmp{s}{alpha} = Support_Dalpha{s}{alpha}(floor(param_dde.Sd^2/2)+1) ;
Support_Dalpha2_tmp{s}{alpha} = Support_Dalpha_{s}{alpha}(floor(param_dde.Sd^2/2)+1) ;
end
end
param_die.Support_Dalpha1 = Support_Dalpha1_tmp ;
param_die.Support_Dalpha2 = Support_Dalpha2_tmp ;

% -------------------------------------------------------------------------






% -------------------------------------------------------------------------
% Generate convolution matrices
% -------------------------------------------------------------------------
G0=sparse(param_data.M,prod(param_im.Nt));
for p = 1:param_data.M
mask = zeros(param_im.Ni) ;
mask(param_data.Vabd(p)+1,param_data.Uabd(p)+1) = 1 ;
maskfzero = zeros(param_im.Nt) ;
maskfzero(param_im.Nadd/2+1:end-param_im.Nadd/2, param_im.Nadd/2+1:end-param_im.Nadd/2) = mask ;
maskf = fftshift(maskfzero) ;
r = find(maskf(:)) ;
if max(size(r))>1
    disp('Erreur echantillonnage')
end
G0(p,r)=1;
end
G0 = sparse(G0) ;
    
dummyG0=transpose(G0);%centering G_gridd and transpose
for i=1:param_data.M
      dummyG0(:,i)=reshape(fftshift((reshape(dummyG0(:,i),param_im.Nt))), prod(param_im.Nt),1);
end
dummyG0=sparse(dummyG0);


param_data.G0 = G0 ;
param_data.SFourier = dummyG0 ;
% -------------------------------------------------------------------------




end