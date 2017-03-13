function [TF, TF_adj,Nt,Nadd] = TF_operators2(Sd,Ni)


Nadd= (2*(2*Sd-1)) ; 

Nu=Ni+Nadd;
Nt=Nu;
disp(['Zero Padding of Fx with ',num2str((Nadd)), ' pixels.'])
disp ' '


TF = @(x)fft2_new(x, Nt,Nadd) ;
TF_adj = @(x) fft2_adj_new(x,Nt,Nadd) ;



end