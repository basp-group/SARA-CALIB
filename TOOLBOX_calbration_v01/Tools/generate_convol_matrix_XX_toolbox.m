function XX = generate_convol_matrix_XX_toolbox(param_im, param_dde)

xhat = reshape(fftshift(reshape(param_im.TF(param_im.x),param_im.Nt)), prod(param_im.Nt),1) ;
XX = sparse(prod(param_im.Nt),prod(param_im.Nt)) ;

Xtemp = cell(length(param_dde.Support_tot)) ;
Supp_tot_2D = convert_1D_2D(param_dde.Support_tot, param_im.Nt(1)) ;
parfor t = 1:length(param_dde.Support_tot)

ind_bar = Supp_tot_2D(t,:) ;

Supp_temp = [ Supp_tot_2D(:,1) + ind_bar(1), Supp_tot_2D(:,2) + ind_bar(2) ] ;
    
ind_select = (Supp_temp(:,1) >= - param_im.Nt(1)/2) ...
            &(Supp_temp(:,2) >= - param_im.Nt(1)/2) ...
            &(Supp_temp(:,1) <= param_im.Nt(1)/2-1) ...
            &(Supp_temp(:,2) <= param_im.Nt(1)/2-1) ;

Supp_star = Supp_temp(ind_select,:) ;
Supp_bar = Supp_tot_2D(ind_select,:) ;
Supp_star_vect = convert_2D_1D(Supp_star,param_im.Nt(1)) ;
Supp_bar_vect = convert_2D_1D(Supp_bar,param_im.Nt(1)) ;

Xtemp{t} = zeros(1,prod(param_im.Nt)) ;
Xtemp{t}(Supp_bar_vect) = xhat(Supp_star_vect) ;
end
Xtemp = cell2mat(Xtemp) ;
XX(param_dde.Support_tot,:) = sparse(Xtemp) ;
end