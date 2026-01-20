%% Phase field solve
vec_ = 0 ;
Const_method = 'Penalty' ; 
n_Cr = 5000 ; % number of penalty points along crack 
w_fac = 1e7 ; % penatly factor for imposing phase field = 1 along crack
P_order = 1;  

G_ = 2.7 ; 

NR.iter      = 30;
NR.tolerance = 1e-10;

seg   = xCr(2,:) - xCr(1,:);   % tip segment

numelem = size(element,1); 
element_app = element ; 
method =  'Lobatto' ;
elemType = 'Q4'  ; 

crack_nodes = []  ;
for i_ = 1 : n_Cr
    crack_nodes(i_,:) = xCr(1,:)+seg*i_/n_Cr;
end
crack_nodes(crack_nodes(:,1)>L,:) = []  ; 

desir_elems = split_elem ; 

numnode = size(node,1) ; 

%% find local coordinate of each node
crack_elems_ = crack_nodes(:,1)*0 ; 
crack_local = crack_nodes*0 ; 
for ip = 1 : length(crack_nodes)
    xx = crack_nodes(ip,:);
    
    for ii = 1 : length(desir_elems)
            sc = node(element(desir_elems(ii),:),:);
            min_x = min(sc(:,1));
            max_x = max(sc(:,1));
            min_y = min(sc(:,2));
            max_y = max(sc(:,2));
            GG.type = 'Q4';
            if xx(1)>=min_x && xx(1)<=max_x && xx(2)>=min_y && xx(2)<=max_y
                [local] = element_coordinate (xx , element(desir_elems(ii),:) , node , GG);
                crack_elems_(ip) = desir_elems(ii);
                crack_local(ip,:) = local ; 
            end
    end
    
end


%%
% compute initial phase field 
Compute_Penalty  ;


nDofs = numnode;
u   = zeros( nDofs + strcmp(Const_method,'Lagrange')*n_Cr , 1 );
du  = zeros( nDofs + strcmp(Const_method,'Lagrange')*n_Cr, 1         );
ddu = zeros( nDofs + strcmp(Const_method,'Lagrange')*n_Cr, 1         );

for iNR = 1 : 2

% Initialize Global Tangent Matrix K and Global RHS
    K   = sparse( nDofs , nDofs);
    RHS = sparse( nDofs, 1 );   

%% ASSEMBLE
if vec_ == 0 
volume_ = 0 ; imat = 0 ; 
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
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('Q9',pt);  % element shape functions
            J0 = node_app(sctr_app,:)'*dNdxi(:,:);                 % element Jacobian matrix
        elseif strcmp(method,'Lobatto')
%             [N,dNdxi] = lobatto_basis(elemType,pt,1,P_order);  % element shape functions
            if order == 2 
                [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,1);
            else
%                 [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
                N = N_p{kk}{1};
                dNdxi = N_p{kk}{2};
            end
            
            J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
            N=N(idx);
            dNdxi=dNdxi(idx,:); 
        end
        
        xgp = N(1:4)'*node(sctr,:) ;

        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
         
        B_(1,1:1:1*nn)      = dNdx(:,1)';
        B_(2,1:1:1*nn) = dNdx(:,2)';
        
        
%         get phase field parameter at this gp
%         phi = N'*u(sctr_phi) ;
%         phi_grad = B_*u(sctr_phi) ;
        phi = N'*u(sctr) ;
        dphi = dNdx'*u(sctr);
        gamma_ = (1/2/l_)*phi*phi + l_/2*(dphi'*dphi) ;
        imat = imat + gamma_ * W(kk)*det(J0) ;

        sai_0 = 0 ; 
        Kphi_phi = ((2*sai_0 + G_/l_)*N*N' + G_*l_*B_'*B_)*W(kk)*det(J0);
        fphi = ( B_'*B_*G_*l_*u(sctr_phi) + N*(phi*(G_/l_+2*sai_0)-2*sai_0) )*W(kk)*det(J0);

        K(sctr_phi,sctr_phi) = K(sctr_phi,sctr_phi) + Kphi_phi ;
        RHS(sctr_phi) = RHS(sctr_phi) + fphi ; 
        
        volume_ = volume_ +  W(kk)*det(J0) ; 
    end                  % end of looping on GPs
end                      % end of looping on elements


end
%         fphi_P = Kphi_phi*u

%         Assemble_staggered  ; 
        KK = K ; 
        
%% External Boundary Conditions
if strcmp(Const_method,'Penalty')
    K = K + KW ;
    RHS = RHS + RW ;
    A = K ; 
    r = RHS+KW*u ; 
elseif strcmp(Const_method,'Lagrange')
    K = K + KW ;
    RHS = RHS + RW ;

    A = [K K_A' ; K_A sparse(n_Cr,n_Cr)];     
    r = [RHS ; -(iNR==1)*R_A ] ;
else
    A = K ; r = RHS ;
    for  ui = 1 : length(crack_nodes) 
        cur_node1 = crack_nodes ( ui )  ;
        [A,r] = boundary_1point(A,r,cur_node1,-(iNR==1)*1);
    end
    
    % zero out extra degrees of freedom 
    if P_order>1
        crackNodes_app = find(node_app(:,2)==0.5 & node_app(:,1)>=0.25 & node_app(:,1)<=0.75) ;
        for  ui = 1 : length(crackNodes_app) 
    %             ui
            cur_node2 = crackNodes_app ( ui )  ;
    %         [A,r] = boundary_1point(A,r,cur_node2,0);
        end
    end
    
    
end

        ddu = A\r;
        
% update solution
        du = du - ddu ;


%% Check convergence?!
        % First iteration setup
        if iNR == 1
            ddu1 = ddu;
        end

        % Residual for generic iteration
        resD = norm( ddu ) / norm( ddu1 );
        resR = norm( r ) ;
        
%         resU = norm( ddu(1:2*numnode ) );
        resPhi = norm( ddu ); 
        
        fprintf(' step %i/%i, itr %i, res displ %i, res rhs %i ,  resPhi %i\n',...
                    1, 1 , iNR, resD, resR, resPhi);

        % Update increment of Nodal Variables
%         du = du + ddu;

        % Update Nodal variables
        u( :, 1 ) = du;
     
         
        % Tolerance check
        if  resD < NR.tolerance            
            break
        elseif iNR == NR.iter
            warning('no convergence')
            return
        end

        if vec_ == 1 
            break ;
        end
max(du)


pause(0.1)

end


PHI = u ; 
% return
% Plot_phase_field ; 
%% 
% %{
ccc = 1 ; imat = 0 ; 
stress = zeros(numelem*4,6) ; 
QQQ = zeros(numelem*4,2) ;
for iel = 1 : numelem 
    if rem(iel,1000)==0
    iel/numelem 
    end

    sctr = element(iel,:); % element connectivity
%     sctr_app = element_app(iel,:); % element connectivity
%     nn   = length(sctr_app);   % number of nodes per element

    sctr_app = element_app(iel,:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;
    nn   = length(sctr_app);   % number of nodes per element

    
    if length(sctr_app)>4
        % Gauss quadrature
        order = n_order ;
    else
        order = 2 ;
    end    
        [W,Q] = quadrature_v2(order,'GAUSS',2);
    
    % scatter vector for element assembly
%     sctrB = [sctr_app sctr_app+numnode];
    sctrB = zeros(1,length(sctr_app)*2) ;
    sctrB(1:2:end) = 2*sctr_app-1 ; 
    sctrB(2:2:end) = 2*sctr_app ;
    sctr_phi = 0*numnode+sctr_app ; 

%     U = zeros(2*length(sctr_app),1) ; 
%     U(1:2:end) = u(2*sctr_app-1) ; 
%     U(2:2:end) = u(2*sctr_app) ; 
        
    
    B = zeros(3,2*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    KK = 0 ; 
    

    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('Q9',pt);  % element shape functions
            J0 = node_app(sctr_app,:)'*dNdxi(:,:);                 % element Jacobian matrix
        elseif strcmp(method,'Lobatto')
            if order == 2 
                [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,1);
            else
%                 [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
                N = N_p{kk}{1};
                dNdxi = N_p{kk}{2};
            end
            J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
            N=N(idx);
            dNdxi=dNdxi(idx,:); 
        end
        xgp = N(1:4)'*node(sctr,:) ; 
        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
%         % B matrix
%         B(1,1:2:2*nn)      = dNdx(:,1)';
%         B(2,2:2:2*nn) = dNdx(:,2)';
%         B(3,1:2:2*nn)      = dNdx(:,2)';
%         B(3,2:2:2*nn) = dNdx(:,1)';
%         
        strain = [0;0;0];
        
%         get phase field parameter at this gp
        phi = N'*u(sctr_phi) ;
        dphi = dNdx'*u(sctr_phi);

        gamma_ = (1/2/l_)*phi*phi + l_/2*(dphi'*dphi) ;
        imat = imat + gamma_ * W(kk)*det(J0) ;

%         dispx = N'* u(2*sctr_app-1) ;
%         dispy = N'* u(2*sctr_app-0) ;
        % Then stress from strain via constitutive equation
        stress(ccc,:) = [ [D*strain]' 0 gamma_ phi] ;
        QQQ(ccc,:)  = xgp ;
        ccc = ccc + 1 ; 
    end                  % end of looping on GPs
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK;

end                      % end of looping on elements

facc_ = imat/L

% volume_
% for ip = 1 : size(crack_nodes,1)
%     if rem(ip,1000)==0
%         ip/size(crack_nodes,1)
%     end
%     
%     xx = crack_nodes(ip,:);
%     
%     
% %% find parent element and local coordinate
%     local = crack_local(ip,:);
%     sctr = element(crack_elems_(ip),:) ; 
%     sctr_app = element_app(crack_elems_(ip),:); % element connectivity
%     idx = find(sctr_app~=0);
%     sctr_app(sctr_app==0)=[] ;
%         
%     xgp = N(1:4)'*node(sctr,:) ; 
% 
% %     [N,~] = lobatto_basis(elemType,local,1,P_order);  % element shape functions
%     [N,~] = lobatto_basis_v2(elemType,local,1,P_order);  % element shape functions
%     
%     if norm(xx-N(1:4)'*node(sctr,:)) > 1e-10
%         error('errr')
%     end
%     
% %     Penalty constraint
%     stress = [ stress;  [D*strain]' 0 0 N'*u(sctr_app)] ;
%     QQQ = [QQQ  ; xgp] ;
% 
%         
%     
% end

tri = delaunay(QQQ(:,1),QQQ(:,2));
% tri = tricheck(QQQ(:,1:2),tri) ;
%%
% ! rm -r out/
% ! mkdir out/

xx = QQQ ;    
if size(xx,2)==2
    xx(:,3) = 0;
end 
yy = tri ;        
VTU_ (xx,yy,stress,'out/all.vtu'); 

xx = node ;    
if size(xx,2)==2
    xx(:,3) = 0;
end 

VTU_2 (xx,element,u,'out/phi_node.vtu'); 


!paraview out/all.vtu&

 %}





























