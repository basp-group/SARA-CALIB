function G = create_matrix_G_D1_D2_toolbox(param_data, param_dde, param_im) 

% *************************************************************************
% Generate convolution matrix for the image
% G = Mask * D1 * conj(D2)
% *************************************************************************
% -- Code slow with Matlab --
% *************************************************************************

% -------------------------------------------------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% d1 and d2 cells
Perm_antennas1 = param_data.Perm_antennas(1:end/2,1) ;
Perm_antennas2 = param_data.Perm_antennas(1:end/2,2) ;
d1 = cell(param_data.T,1) ;
d2 = cell(param_data.T,1) ;
U1 = param_dde.U1 ;
U2 = param_dde.U2 ;
Ind_conj = param_dde.Ind_conj ;
parfor s = 1:param_data.T
D1 = cell2mat(U1{s}) ;
d1{s} = D1(Perm_antennas1,:) ;
D2 = cell2mat(U2{s}) ;
d2{s} = conj(D2(Perm_antennas2,Ind_conj)) ;
end
d1 = cell2mat(d1) ;
d2conj = transpose(cell2mat(d2)) ;
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------

% -------------------------------------------------------------
% convolution Gchirp et D2_conj -------------------------------
% -------------------------------------------------------------
G1 = cell(1,prod(param_im.Nt));
parfor t = 1:prod(param_im.Nt)
G1{t} = transpose(sum(param_data.SFourier(param_dde.convol_pos{t},:).*d2conj(param_dde.Elem_D_convol{t},:),1)) ;
end 
% -------------------------------------------------------------
% -------------------------------------------------------------

% -------------------------------------------------------------
% Creation matrice G ------------------------------------------
% -------------------------------------------------------------
G= cell(1,prod(param_im.Nt));
parfor t = 1:prod(param_im.Nt)
G{t} = (sum(cell2mat({G1{param_dde.convol_pos{t}}}).*d1(:, param_dde.Elem_D_convol{t}),2)) ;
end 
G = cell2mat(G) ;
G = G(:,param_im.Index_tshift) ;
G = sparse(G) ;
% -------------------------------------------------------------

end



