function [param_dde, param_die, param_data] = ...
            Generate_ddes(param_dde, param_die, param_data, param_im)


center = floor(param_dde.Sd/2)*param_dde.Sd + floor(param_dde.Sd/2) + 1 ;
D_kernel = cell(param_data.T,1);
DIE_kernel = D_kernel ;

for s = 1:param_data.T
D_kernel{s} = param_dde.var*randn(param_data.na,param_dde.Sd^2) ...
    + 1i * param_dde.var*randn(param_data.na,param_dde.Sd^2) ;
D_kernel{s}(:,center) = ...
              sign(randn(param_data.na,1)) + param_dde.var*randn(param_data.na,1) + ...
              1i * (sign(randn(param_data.na,1)) + param_dde.var*randn(param_data.na,1)) ;
DIE_kernel{s} = D_kernel{s}(:,center) ;

D_kernel{s} = mat2cell(D_kernel{s},ones(param_data.na,1),param_dde.Sd^2) ;
DIE_kernel{s} = num2cell(DIE_kernel{s}) ;
end

param_dde.U_true = D_kernel ;
param_die.U_true = DIE_kernel ;


i0 = param_im.Nt/2+1 ;
i_min = i0-(param_dde.Sd-1)/2 ;
i_max = i0 + (param_dde.Sd-1)/2 ;

temp = zeros(param_im.Nt) ;
temp(i_min:i_max, i_min:i_max) = 1 ;
param_dde.Support_D = find(temp(:)~=0) ;

temp = zeros(param_im.Nt) ;
temp(i0, i0) = 1 ;
param_die.Support_D = find(temp(:)~=0) ;




param_data.Perm_antennas = ...
    [param_data.antennas(1:param_data.na*(param_data.na-1)/2,:) ; ...
    param_data.antennas(end-( param_data.na*(param_data.na-1)/2 ) +1:end, :)] ;















% -------------------------------------------------------------------
% Recuperation des indices pour les convolutions
% -------------------------------------------------------------------
Ntot = prod(param_im.Nt) ;
Elem_D = zeros(param_dde.Sd^2, Ntot) ;
Elem_D_convol = cell(1,Ntot) ;
Elem_D_convol_die = cell(1,Ntot) ;
convol_pos_D = cell(1,Ntot) ;
Sd2 = floor(param_dde.Sd/2) ;
% -----------------------------------
t = 1 ;
for j = 1:param_im.Nt(2)
for i = 1:param_im.Nt(1)
temp = zeros(param_dde.Sd,param_dde.Sd) ;
if Sd2>=i
temp(1:Sd2+1-i,:) = 1 ;
end
if param_im.Nt(1)-i<Sd2
    temp(Sd2+2+param_im.Nt(1)-i:end,:) = 1 ; 
end
if Sd2>=j
  temp(:,1:Sd2+1-j) = 1 ;
end
if param_im.Nt(2)-j<Sd2
    temp(:,Sd2+2+param_im.Nt(2)-j:end) = 1 ; 
end
convol_pos_D{t} = param_dde.Support_D(temp(:)==0) ;
Elem_D(:,t) = temp(:)==0 ;
Elem_D_convol{t} = find(temp(:)==0) ;
Elem_D_convol_die{t} = 1 ;
t = t+1 ;
end
end
Elem_D = sparse(Elem_D) ;

param_dde.Elem_D_convol = Elem_D_convol ;
param_die.Elem_D_convol = Elem_D_convol_die ;



% -----------------------------------
% Indices convolution
% -----------------------------------
convol_pos = cell(1,Ntot) ;
convol_pos_die = cell(1,Ntot) ;
temp = zeros(Ntot,1) ;
temp(param_dde.Support_D) = 1 ;
temp = sparse(temp) ;
first = min(param_dde.Support_D) ;
last = max(param_dde.Support_D) ;
z0=(param_im.Nt(1)+1)*param_im.Nt(1)/2+1;
for t = 1 : Ntot
kernels = shift_kernels(temp,t,z0,param_im.Nt(1),first,last,Ntot) ;
convol_pos{t} = find(kernels) ;
convol_pos_die{t} = t ;
end

param_dde.convol_pos = convol_pos ;
param_die.convol_pos = convol_pos_die ;


% -----------------------------------
% Indices pour le conjugue de d_alpha
% -----------------------------------
Ind_conj = 1:param_dde.Sd^2 ;
Ind_conj = fliplr(Ind_conj) ;

param_dde.Ind_conj = Ind_conj ;
param_die.Ind_conj = 1 ;
% -------------------------------------------------------------------
% -------------------------------------------------------------------







end