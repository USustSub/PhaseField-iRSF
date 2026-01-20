if strcmp(Const_method,'Penalty') || strcmp(Const_method,'Lagrange')

%% add Penalty stiffness for imposing theta = 1
KW = sparse(numnode,numnode);
RW = zeros(numnode,1);
K_A =sparse(size(crack_nodes,1),1*numnode);
R_A =sparse(size(crack_nodes,1),1);
O = sparse(size(crack_nodes,1),size(crack_nodes,1));
normal = null(seg)'; 
for ip = 1 : size(crack_nodes,1)
    if rem(ip,100)==0
        ip/size(crack_nodes,1)
    end
    
    xx = crack_nodes(ip,:);
    
    
%% find parent element and local coordinate
    local = crack_local(ip,:);
    sctr = element(crack_elems_(ip),:) ; 
    sctr_app = element_app(crack_elems_(ip),:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;

%     [N,~] = lobatto_basis(elemType,local,1,P_order);  % element shape functions
    [N,dNdxi] = lobatto_basis_v2(elemType,local,1,P_order);  % element shape functions
            
    J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
    N=N(idx);
    dNdxi=dNdxi(idx,:); 
    invJ0 = inv(J0);
    dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
    
    if norm(xx-N(1:4)'*node(sctr,:)) > 1e-10
        error('errr')
    end
    
%     Penalty constraint
    kw = N*N';
    rw = N; 
    KW(sctr_app,sctr_app) = KW(sctr_app,sctr_app) + kw'*w_fac ;
    RW(sctr_app) = RW(sctr_app) - rw*w_fac ;
    
%     kw = dN(:,1)*dN(:,1)';
%     rw = dN(:,1); 
%     KW(sctr_app,sctr_app) = KW(sctr_app,sctr_app) + kw'*w_fac ;
%     RW(sctr_app) = RW(sctr_app) - rw*w_fac ;
    
%     kw = dN(:,2)*dN(:,2)';
%     rw = dN(:,2); 
%     KW(sctr_app,sctr_app) = KW(sctr_app,sctr_app) + kw'*w_fac ;
%     RW(sctr_app) = RW(sctr_app) - rw*w_fac ;
    
% %     K_A(ip,0*numnode+sctr_app) = K_A(ip,0*numnode+sctr_app)  + N' ;  
% %     R_A(ip) = 1 ; 
    
    M = dNdx*normal';
    K_A(ip,0*numnode+sctr_app) = K_A(ip,0*numnode+sctr_app)  + M' ;  
    R_A(ip) = 0 ; 

    
end

end




