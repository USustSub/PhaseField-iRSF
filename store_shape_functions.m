 % ---------------------
    % Loop on Gauss points   
    % ---------------------
[W,Q] = quadrature_v2(n_order,'GAUSS',2);
% P_order = 38;
N_p = [];
for kk = 1 : size(W,1)
    if rem(kk,200)==0
        kk/size(W,1)
    end
        pt = Q(kk,:);                             % quadrature_v2 point
        
        % Shape functions and its derivatives
%             [N,dNdxi] = lobatto_basis(elemType,pt,1,P_order);  % element shape functions
        [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
        N_p{kk}{1} = N ; 
        N_p{kk}{2} = dNdxi ; 
        
end 


% save data_p5_nq40.mat 'N_p'
eval(['save data_p' num2str(P_order) '_nq' num2str(n_order) '.mat N_p'])
