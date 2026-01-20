all = zeros(numelem*4,17); 
% all = [] ;
cc = 1 ;
for iel = 1 : numelem 
    if rem(iel,1000)==0
        iel/numelem 
    end
    sctr = element(iel,:); % element connectivity
    sctr_app = element_app(iel,:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;
    nn   = length(sctr_app);   % number of nodes per element

    % scatter vector for element assembly
    sctr_phi = sctr_app ; 

    B = zeros(3,2*nn);
    B_ = zeros(2,1*nn);
    
%     if ismember(iel,desir_elems)
    if length(sctr_app)>4
        % Gauss quadrature_v2
        order = n_order ;
    else
        order = 2 ;
    end
%         [W,Q] = quadrature_v2(order,'GAUSS',2);
        [W,Q] = quadrature_v2( order, 'GAUSS', 2 );
    
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    for kk = 1 : size(W,1)
        if rem(kk,200)==0
            kk/size(W,1)
        end
        pt = Q(kk,:);                             % quadrature_v2 point
        
        % Shape functions and its derivatives
%         if iel == 1 
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('Q9',pt);  % element shape functions
            J0 = node_app(sctr_app,:)'*dNdxi(:,:);                 % element Jacobian matrix
        elseif strcmp(method,'Lobatto')
%             [N,dNdxi] = lobatto_basis(elemType,pt,1,P_order);  % element shape functions
%             if order == 2 
%                 [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,1);
%             else
%                 [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
                N = N_p{kk}{1};
                dNdxi = N_p{kk}{2};
%             end
            
            J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
            N=N(idx);
            dNdxi=dNdxi(idx,:); 
        end
%         end
        xgp = N(1:4)'*node(sctr,:) ;

        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
         all(cc,:) = [dNdx(:)'  N' W(kk)*det(J0) sctr_phi];
         cc = cc + 1 ;

    end                  % end of looping on GPs
end                      % end of looping on elements
