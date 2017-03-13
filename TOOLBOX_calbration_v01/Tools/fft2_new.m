function y= fft2_new(x, No,zeropadd)
%%

% Add zeros in image domain
X = zeros(No) ;
X(zeropadd/2+2:end-zeropadd/2,zeropadd/2+2:end-zeropadd/2) = x ;
% 2D FT
Y = fft2(X) ;

y = Y(:) ;
y(abs(y)<1e-16) = 0 ;
end