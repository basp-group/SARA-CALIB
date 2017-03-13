function x = fft2_adj_new(y,Nt,zeropadd)
%%
Y = reshape(y,Nt) ;
X = ifft2(Y) ;
x = X(zeropadd/2+2:end-zeropadd/2,zeropadd/2+2:end-zeropadd/2) ;


end
