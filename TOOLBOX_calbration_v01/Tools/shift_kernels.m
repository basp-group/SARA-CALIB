function kernels = shift_arwa(G1,t,z0,nx,first,last,Ntot)
% G1: matrice a decaler
% t: position actuelle dans l'image
% z0: milieu de l'image
% nx^2: taille de la TF
% first: premier indice non nul
% last: dernier indice non nul
% pix^2: taille image originale x 4

kernels=G1;
pos=t-z0;
sg=sign(pos);
jpos=mod(t,nx);% column index

if sg==-1 && ~isempty(nonzeros(first+pos<1))
    if first+pos<1
    kernels(1:-pos,:)=0;
    kernels=sparse(kernels);
    end
elseif sg==1 && ~isempty(nonzeros(last+pos>Ntot)) 
    if last+pos>Ntot
    kernels(Ntot-pos+1:Ntot,:)=0; 
    kernels=sparse(kernels);   
    end
end

% now circular shifting  of  sparse kernels 
kernels=(circshift(kernels,[pos 0]));

if jpos == 0
    jpos = nx ;
end

temp=zeros(nx,nx);
temp(1:jpos-nx/2, : ) = 1 ;
temp(jpos+nx/2 : end, : ) = 1 ;
kernels(temp(:)>0,:) = 0 ;

    
end