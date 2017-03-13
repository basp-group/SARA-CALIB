function n = convert_2D_1D(c, N)

% c : position 2D dans {-N/2, ..., N/2-1}^2
% n : position 1D dans {1,...,N^2}
% N : taille de l image = NxN

n = ( c(:,1) + N/2 )*N + c(:,2) + N/2 +1 ;

end