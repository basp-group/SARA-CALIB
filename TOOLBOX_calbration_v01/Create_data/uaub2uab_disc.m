function [uab, vab, Ant] = uaub2uab_disc(na,u,v,nstep, Nx) 

Ant = cell(nstep,2) ;
uab = cell(nstep,1) ;
vab = cell(nstep,1) ;


for s = 1:nstep ;
    
Ant{s}{1} = [] ;
Ant{s}{2} = [] ;
uab{s} = [] ;
vab{s} = [] ;

for alpha = 1:na-1
Ant{s}{1} = [Ant{s}{1} ; alpha * ones((na-alpha),1)] ; % antenne alpha
Ant{s}{2} = [Ant{s}{2} ; (alpha+1:na)'] ; % antenne beta
for beta = alpha+1:na
uab{s} = [uab{s} ; u{alpha}(s) - u{beta}(s)] ; % u_alpha - u_beta
vab{s} = [vab{s} ; v{alpha}(s) - v{beta}(s)] ; % v_alpha - v_beta
end
end
uab{s} = uab{s} + Nx/2 ; % recentrage
vab{s} = vab{s} + Nx/2 ; % recentrage
end





% % % Uab = cell2mat(uab) ;
% % % Vab = cell2mat(vab) ;
% % % figure
% % % plot(Uab, Vab, '.')

end