function c = convert_1D_2D(n,N)

% n : position 1D dans {1,...,N^2}
% c : position 2D dans {-N/2, ..., N/2-1}^2
% N : taille de l image = NxN

d = floor(n/N) ;
c = [ d- (N/2) , n- (d*N) - (N/2) - 1 ] ;

end