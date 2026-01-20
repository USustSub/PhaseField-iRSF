% Inputs
% ---------------------------------------
% Mohsen Goudarzi (July 2021 - Utrecht ) 
% ---------------------------------------

addpath ('PFEM_functions')

clear all ; close all ; 
clc
tic;     
disp ('P-Finite element method for phase field ' )

profile_ = 0  ; % if =1, matlab profile on 

if profile_ == 1
    ! rm -r profile_results/
    profile on
end

Inputs ; 

Initialize_ ;

% compute initial phase field 
Compute_Penalty  ;


if vec_ == 1
    Pre_compute_values ;
end
%% Phase field solve
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
volume_ = 0 ; 
for iel = 1 : numelem 
    if rem(iel,10)==0
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
        phi = N'*u(sctr_phi) ;
        phi_grad = B_*u(sctr_phi) ;
        
        sai_0 = 0 ; 
        Kphi_phi = ((2*sai_0 + G/l_)*N*N' + G*l_*B_'*B_)*W(kk)*det(J0);
        fphi = ( B_'*B_*G*l_*u(sctr_phi) + N*(phi*(G/l_+2*sai_0)-2*sai_0) )*W(kk)*det(J0);

        K(sctr_phi,sctr_phi) = K(sctr_phi,sctr_phi) + Kphi_phi ;
        RHS(sctr_phi) = RHS(sctr_phi) + fphi ; 
        
        volume_ = volume_ +  W(kk)*det(J0) ; 
    end                  % end of looping on GPs
end                      % end of looping on elements



        
% or in vectorized format
%%
else
    
%{
    K   = sparse( nDofs , nDofs);
    RHS = sparse( nDofs, 1 );   
Cell = [] ; 
volume_ = 0 ; 
for iel = 1 : numelem 
    if rem(iel,10)==0
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
    Kphi_phi = 0  ;
    fphi  =  0 ; 
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
        phi = N'*u(sctr_phi) ;
        phi_grad = B_*u(sctr_phi) ;
        
        sai_0 = 0 ; 
        Kphi_phi = Kphi_phi + ((2*sai_0 + G/l_)*N*N' + G*l_*B_'*B_)*W(kk)*det(J0);
        fphi = fphi + ( B_'*B_*G*l_*u(sctr_phi) + N*(phi*(G/l_+2*sai_0)-2*sai_0) )*W(kk)*det(J0);

        K(sctr_phi,sctr_phi) = K(sctr_phi,sctr_phi) + ((2*sai_0 + G/l_)*N*N' + G*l_*B_'*B_)*W(kk)*det(J0);
        
        
    end                  % end of looping on GPs
    
    

%         K(sctr_phi,sctr_phi) = K(sctr_phi,sctr_phi) + Kphi_phi ;
        RHS(sctr_phi) = RHS(sctr_phi) + fphi ; 
        
        volume_ = volume_ +  W(kk)*det(J0) ; 
        
        [i j s] = find(Kphi_phi) ; 
%         Cell = [Cell ; [sctr_phi(i)' ,sctr_phi(j)', s]];

        Cell{iel} = [sctr_phi(i)' ,sctr_phi(j)', s];
    
end                      % end of looping on elements

            IJV = cell2mat( Cell' );
    K_2 = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,2));
%}
    dNdx1 = all(:,1);     dNdx2 = all(:,2);     dNdx3 = all(:,3);     dNdx4 = all(:,4); 
    dNdy1 = all(:,5);     dNdy2 = all(:,6);     dNdy3 = all(:,7);     dNdy4 = all(:,8); 
    N1 = all(:,9);     N2 = all(:,10);     N3 = all(:,11);     N4 = all(:,12); 
    wq =     all(:,13) ;     sctrB = all(:,14:17);
    kpp = [] ; 
    kpp(:,1) = wq.*((G*N1.^2)/l_ + G*l_.*dNdx1.^2 + G*l_*dNdy1.^2);
    kpp(:,2) =wq.*((G.*N1.*N2)/l_ + G.*dNdx1.*dNdx2.*l_ + G.*dNdy1.*dNdy2.*l_);
    kpp(:,3) =wq.*((G.*N1.*N3)/l_ + G.*dNdx1.*dNdx3.*l_ + G.*dNdy1.*dNdy3.*l_);
    kpp(:,4) =wq.*((G.*N1.*N4)/l_ + G.*dNdx1.*dNdx4.*l_ + G.*dNdy1.*dNdy4.*l_);
    kpp(:,5) =wq.*((G.*N1.*N2)/l_ + G.*dNdx1.*dNdx2.*l_ + G.*dNdy1.*dNdy2.*l_);
    kpp(:,6) =         wq.*((G.*N2.^2)/l_ + G.*l_.*dNdx2.^2 + G.*l_.*dNdy2.^2);
    kpp(:,7) =wq.*((G.*N2.*N3)/l_ + G.*dNdx2.*dNdx3.*l_ + G.*dNdy2.*dNdy3.*l_);
    kpp(:,8) = wq.*((G.*N2.*N4)/l_ + G.*dNdx2.*dNdx4.*l_ + G.*dNdy2.*dNdy4.*l_);
    kpp(:,9) =  wq.*((G.*N1.*N3)/l_ + G.*dNdx1.*dNdx3.*l_ + G.*dNdy1.*dNdy3.*l_);
    kpp(:,10) = wq.*((G.*N2.*N3)/l_ + G.*dNdx2.*dNdx3.*l_ + G.*dNdy2.*dNdy3.*l_);
    kpp(:,11) =          wq.*((G.*N3.^2)/l_ + G.*l_.*dNdx3.^2 + G.*l_.*dNdy3.^2);
    kpp(:,12) = wq.*((G.*N3.*N4)/l_ + G.*dNdx3.*dNdx4.*l_ + G.*dNdy3.*dNdy4.*l_);
    kpp(:,13) = wq.*((G.*N1.*N4)/l_ + G.*dNdx1.*dNdx4.*l_ + G.*dNdy1.*dNdy4.*l_);
    kpp(:,14) = wq.*((G.*N2.*N4)/l_ + G.*dNdx2.*dNdx4.*l_ + G.*dNdy2.*dNdy4.*l_);
    kpp(:,15) = wq.*((G.*N3.*N4)/l_ + G.*dNdx3.*dNdx4.*l_ + G.*dNdy3.*dNdy4.*l_);
    kpp(:,16) =          wq.*((G.*N4.^2)/l_ + G.*l_.*dNdx4.^2 + G.*l_.*dNdy4.^2);

    S_all = zeros(numelem*4*16,3);
    S_all(1:16:end,1) = sctrB(:,1);
    S_all(1:16:end,2) = sctrB(:,1);
    S_all(1:16:end,3) = kpp(:,1) ;
    
    S_all(2:16:end,1) = sctrB(:,2);
    S_all(2:16:end,2) = sctrB(:,1);
    S_all(2:16:end,3) = kpp(:,2) ;

    S_all(3:16:end,1) = sctrB(:,3);
    S_all(3:16:end,2) = sctrB(:,1);
    S_all(3:16:end,3) = kpp(:,3) ;

    S_all(4:16:end,1) = sctrB(:,4);
    S_all(4:16:end,2) = sctrB(:,1);
    S_all(4:16:end,3) = kpp(:,4) ;
    
    S_all(5:16:end,1) = sctrB(:,1);
    S_all(5:16:end,2) = sctrB(:,2);
    S_all(5:16:end,3) = kpp(:,5) ;
    
    S_all(6:16:end,1) = sctrB(:,2);
    S_all(6:16:end,2) = sctrB(:,2);
    S_all(6:16:end,3) = kpp(:,6) ;

    S_all(7:16:end,1) = sctrB(:,3);
    S_all(7:16:end,2) = sctrB(:,2);
    S_all(7:16:end,3) = kpp(:,7) ;

    S_all(8:16:end,1) = sctrB(:,4);
    S_all(8:16:end,2) = sctrB(:,2);
    S_all(8:16:end,3) = kpp(:,8) ;
    
    S_all(9:16:end,1) = sctrB(:,1);
    S_all(9:16:end,2) = sctrB(:,3);
    S_all(9:16:end,3) = kpp(:,9) ;
    
    S_all(10:16:end,1) = sctrB(:,2);
    S_all(10:16:end,2) = sctrB(:,3);
    S_all(10:16:end,3) = kpp(:,10) ;

    S_all(11:16:end,1) = sctrB(:,3);
    S_all(11:16:end,2) = sctrB(:,3);
    S_all(11:16:end,3) = kpp(:,11) ;

    S_all(12:16:end,1) = sctrB(:,4);
    S_all(12:16:end,2) = sctrB(:,3);
    S_all(12:16:end,3) = kpp(:,12) ;
    
    S_all(13:16:end,1) = sctrB(:,1);
    S_all(13:16:end,2) = sctrB(:,4);
    S_all(13:16:end,3) = kpp(:,13) ;
    
    S_all(14:16:end,1) = sctrB(:,2);
    S_all(14:16:end,2) = sctrB(:,4);
    S_all(14:16:end,3) = kpp(:,14) ;

    S_all(15:16:end,1) = sctrB(:,3);
    S_all(15:16:end,2) = sctrB(:,4);
    S_all(15:16:end,3) = kpp(:,15) ;

    S_all(16:16:end,1) = sctrB(:,4);
    S_all(16:16:end,2) = sctrB(:,4);
    S_all(16:16:end,3) = kpp(:,16) ;
    
%     K_3 = sparse(S_all(:,1),S_all(:,2),S_all(:,3),size(K,1),size(K,2));
    
% divide it into blocks (less memory)
b = 200000; % block size
n = numelem.*4.*16;
c_ = diff([0:b:n-1,n]) ; 
K_3 = 0*K ; 
for i_step = 1 : size(c_,2)
    i_step/size(c_,2)
    
    s0 = sum(c_(1:i_step-1));
    s_step = c_(i_step);
    s_des = [s0+1:s0+s_step] ;
    K_3 = K_3 + sparse(S_all(s_des,1),S_all(s_des,2),S_all(s_des,3),size(K,1),size(K,2));
end
    
    
    K = K_3 ; 

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

Plot_phase_field ; 
% return

%% Momentum solve 
nDofs = 2*numnode;
u   = zeros( nDofs , 1 );
du  = zeros( nDofs , 1 );
ddu = zeros( nDofs , 1 );

for iNR = 1 : 2

% Initialize Global Tangent Matrix K and Global RHS
    K   = sparse( nDofs , nDofs);
    RHS = sparse( nDofs, 1 );   

%% ASSEMBLE
gg = [] ; 
volume_ = 0 ; 
if vec_ == 0 
    for iel = 1 : numelem 
    sctr = element(iel,:); % element connectivity
    sctr_app = element_app(iel,:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;
    nn   = length(sctr_app);   % number of nodes per element

%    
    sctrB = zeros(1,length(sctr_app)*2) ;
    sctrB(1:2:end) = 2*sctr_app-1 ; 
    sctrB(2:2:end) = 2*sctr_app ;
        
    sctr_phi = sctr_app ; 


    
%     if ismember(iel,desir_elems)
    if length(sctr_app)>4
        % Gauss quadrature_v2
        order = n_order ;
    else
        order = 2 ;
    end
%         [W,Q] = quadrature_v2(order,'GAUSS',2);
        [W,Q] = quadrature_v2( order, 'GAUSS', 2 );
        
    
    B = zeros(3,2*nn);
    B_ = zeros(2,1*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    for kk = 1 : size(W,1)
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
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
         
        B_(1,1:1:1*nn)      = dNdx(:,1)';
        B_(2,1:1:1*nn) = dNdx(:,2)';
        
        
%         get phase field parameter at this gp
        phi = N'*PHI(sctr_phi) ;
        phi_grad = B_*PHI(sctr_phi) ;
        gg  = [ gg ; phi ] ; 
%         phi = 0  ; 
        Ku_u = ( (1-phi)^2 + kkk )* B'*D*B*W(kk)*det(J0) ; 
        fu   = ( (1-phi)^2 + kkk )* B'*D*B*W(kk)*det(J0)*u(sctrB) ; 
        
        eps = B*u(sctrB);
        sai_0 = 1/2*eps'*D*eps*0; 
        
        
% get phase field at current node        
        K(sctrB,sctrB) = K(sctrB,sctrB) + Ku_u;
        RHS(sctrB) = RHS(sctrB) + fu ; 
                
        volume_ = volume_ +  W(kk)*det(J0) ; 
    end                  % end of looping on GPs
end                      % end of looping on elements
else
    dNdx1 = all(:,1);     dNdx2 = all(:,2);     dNdx3 = all(:,3);     dNdx4 = all(:,4); 
    dNdy1 = all(:,5);     dNdy2 = all(:,6);     dNdy3 = all(:,7);     dNdy4 = all(:,8); 
    N1 = all(:,9);     N2 = all(:,10);     N3 = all(:,11);     N4 = all(:,12); 
    wq =     all(:,13) ;     
    sctrB_ = all(:,14:17);
    sctrB = [] ;
    sctrB(:,[1 3 5 7]) = sctrB_*2-1 ; 
    sctrB(:,[2 4 6 8]) = sctrB_*2-0 ; 
    
    phi = N1.*PHI(sctrB_(:,1)) + N2.*PHI(sctrB_(:,2)) + N3.*PHI(sctrB_(:,3)) + N4.*PHI(sctrB_(:,4));
    
    D11 = D(1,1);
    D12 = D(1,2);
    D13 = D(1,3);
    D21 = D(2,1);
    D22 = D(2,2);
    D23 = D(2,3);
    D31 = D(3,1);
    D32 = D(3,2);
    D33 = D(3,3);
    
    kpp = [] ; 
    
    kpp(:,1) = wq.*(dNdx1.*(D11.*dNdx1.*(kkk + (phi - 1).^2) + D31.*dNdy1.*(kkk + (phi - 1).^2)) + dNdy1.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));;
    kpp(:,2) = wq.*(dNdx1.*(D21.*dNdy1.*(kkk + (phi - 1).^2) + D31.*dNdx1.*(kkk + (phi - 1).^2)) + dNdy1.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,3) = wq.*(dNdx1.*(D11.*dNdx2.*(kkk + (phi - 1).^2) + D31.*dNdy2.*(kkk + (phi - 1).^2)) + dNdy1.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,4) = wq.*(dNdx1.*(D21.*dNdy2.*(kkk + (phi - 1).^2) + D31.*dNdx2.*(kkk + (phi - 1).^2)) + dNdy1.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,5) = wq.*(dNdx1.*(D11.*dNdx3.*(kkk + (phi - 1).^2) + D31.*dNdy3.*(kkk + (phi - 1).^2)) + dNdy1.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,6) = wq.*(dNdx1.*(D21.*dNdy3.*(kkk + (phi - 1).^2) + D31.*dNdx3.*(kkk + (phi - 1).^2)) + dNdy1.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,7) = wq.*(dNdx1.*(D11.*dNdx4.*(kkk + (phi - 1).^2) + D31.*dNdy4.*(kkk + (phi - 1).^2)) + dNdy1.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,8) = wq.*(dNdx1.*(D21.*dNdy4.*(kkk + (phi - 1).^2) + D31.*dNdx4.*(kkk + (phi - 1).^2)) + dNdy1.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,9) = wq.*(dNdy1.*(D12.*dNdx1.*(kkk + (phi - 1).^2) + D32.*dNdy1.*(kkk + (phi - 1).^2)) + dNdx1.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,10) = wq.*(dNdy1.*(D22.*dNdy1.*(kkk + (phi - 1).^2) + D32.*dNdx1.*(kkk + (phi - 1).^2)) + dNdx1.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,11) = wq.*(dNdy1.*(D12.*dNdx2.*(kkk + (phi - 1).^2) + D32.*dNdy2.*(kkk + (phi - 1).^2)) + dNdx1.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,12) = wq.*(dNdy1.*(D22.*dNdy2.*(kkk + (phi - 1).^2) + D32.*dNdx2.*(kkk + (phi - 1).^2)) + dNdx1.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,13) = wq.*(dNdy1.*(D12.*dNdx3.*(kkk + (phi - 1).^2) + D32.*dNdy3.*(kkk + (phi - 1).^2)) + dNdx1.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,14) = wq.*(dNdy1.*(D22.*dNdy3.*(kkk + (phi - 1).^2) + D32.*dNdx3.*(kkk + (phi - 1).^2)) + dNdx1.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,15) = wq.*(dNdy1.*(D12.*dNdx4.*(kkk + (phi - 1).^2) + D32.*dNdy4.*(kkk + (phi - 1).^2)) + dNdx1.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,16) = wq.*(dNdy1.*(D22.*dNdy4.*(kkk + (phi - 1).^2) + D32.*dNdx4.*(kkk + (phi - 1).^2)) + dNdx1.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,17) = wq.*(dNdx2.*(D11.*dNdx1.*(kkk + (phi - 1).^2) + D31.*dNdy1.*(kkk + (phi - 1).^2)) + dNdy2.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,18) = wq.*(dNdx2.*(D21.*dNdy1.*(kkk + (phi - 1).^2) + D31.*dNdx1.*(kkk + (phi - 1).^2)) + dNdy2.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,19) = wq.*(dNdx2.*(D11.*dNdx2.*(kkk + (phi - 1).^2) + D31.*dNdy2.*(kkk + (phi - 1).^2)) + dNdy2.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,20) = wq.*(dNdx2.*(D21.*dNdy2.*(kkk + (phi - 1).^2) + D31.*dNdx2.*(kkk + (phi - 1).^2)) + dNdy2.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,21) = wq.*(dNdx2.*(D11.*dNdx3.*(kkk + (phi - 1).^2) + D31.*dNdy3.*(kkk + (phi - 1).^2)) + dNdy2.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,22) = wq.*(dNdx2.*(D21.*dNdy3.*(kkk + (phi - 1).^2) + D31.*dNdx3.*(kkk + (phi - 1).^2)) + dNdy2.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,23) = wq.*(dNdx2.*(D11.*dNdx4.*(kkk + (phi - 1).^2) + D31.*dNdy4.*(kkk + (phi - 1).^2)) + dNdy2.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,24) = wq.*(dNdx2.*(D21.*dNdy4.*(kkk + (phi - 1).^2) + D31.*dNdx4.*(kkk + (phi - 1).^2)) + dNdy2.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,25) = wq.*(dNdy2.*(D12.*dNdx1.*(kkk + (phi - 1).^2) + D32.*dNdy1.*(kkk + (phi - 1).^2)) + dNdx2.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,26) = wq.*(dNdy2.*(D22.*dNdy1.*(kkk + (phi - 1).^2) + D32.*dNdx1.*(kkk + (phi - 1).^2)) + dNdx2.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,27) = wq.*(dNdy2.*(D12.*dNdx2.*(kkk + (phi - 1).^2) + D32.*dNdy2.*(kkk + (phi - 1).^2)) + dNdx2.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,28) = wq.*(dNdy2.*(D22.*dNdy2.*(kkk + (phi - 1).^2) + D32.*dNdx2.*(kkk + (phi - 1).^2)) + dNdx2.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,29) = wq.*(dNdy2.*(D12.*dNdx3.*(kkk + (phi - 1).^2) + D32.*dNdy3.*(kkk + (phi - 1).^2)) + dNdx2.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,30) = wq.*(dNdy2.*(D22.*dNdy3.*(kkk + (phi - 1).^2) + D32.*dNdx3.*(kkk + (phi - 1).^2)) + dNdx2.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,31) = wq.*(dNdy2.*(D12.*dNdx4.*(kkk + (phi - 1).^2) + D32.*dNdy4.*(kkk + (phi - 1).^2)) + dNdx2.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,32) = wq.*(dNdy2.*(D22.*dNdy4.*(kkk + (phi - 1).^2) + D32.*dNdx4.*(kkk + (phi - 1).^2)) + dNdx2.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,33) = wq.*(dNdx3.*(D11.*dNdx1.*(kkk + (phi - 1).^2) + D31.*dNdy1.*(kkk + (phi - 1).^2)) + dNdy3.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,34) = wq.*(dNdx3.*(D21.*dNdy1.*(kkk + (phi - 1).^2) + D31.*dNdx1.*(kkk + (phi - 1).^2)) + dNdy3.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,35) = wq.*(dNdx3.*(D11.*dNdx2.*(kkk + (phi - 1).^2) + D31.*dNdy2.*(kkk + (phi - 1).^2)) + dNdy3.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,36) = wq.*(dNdx3.*(D21.*dNdy2.*(kkk + (phi - 1).^2) + D31.*dNdx2.*(kkk + (phi - 1).^2)) + dNdy3.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,37) = wq.*(dNdx3.*(D11.*dNdx3.*(kkk + (phi - 1).^2) + D31.*dNdy3.*(kkk + (phi - 1).^2)) + dNdy3.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,38) = wq.*(dNdx3.*(D21.*dNdy3.*(kkk + (phi - 1).^2) + D31.*dNdx3.*(kkk + (phi - 1).^2)) + dNdy3.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,39) = wq.*(dNdx3.*(D11.*dNdx4.*(kkk + (phi - 1).^2) + D31.*dNdy4.*(kkk + (phi - 1).^2)) + dNdy3.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,40) = wq.*(dNdx3.*(D21.*dNdy4.*(kkk + (phi - 1).^2) + D31.*dNdx4.*(kkk + (phi - 1).^2)) + dNdy3.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,41) = wq.*(dNdy3.*(D12.*dNdx1.*(kkk + (phi - 1).^2) + D32.*dNdy1.*(kkk + (phi - 1).^2)) + dNdx3.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,42) = wq.*(dNdy3.*(D22.*dNdy1.*(kkk + (phi - 1).^2) + D32.*dNdx1.*(kkk + (phi - 1).^2)) + dNdx3.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,43) = wq.*(dNdy3.*(D12.*dNdx2.*(kkk + (phi - 1).^2) + D32.*dNdy2.*(kkk + (phi - 1).^2)) + dNdx3.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,44) = wq.*(dNdy3.*(D22.*dNdy2.*(kkk + (phi - 1).^2) + D32.*dNdx2.*(kkk + (phi - 1).^2)) + dNdx3.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,45) = wq.*(dNdy3.*(D12.*dNdx3.*(kkk + (phi - 1).^2) + D32.*dNdy3.*(kkk + (phi - 1).^2)) + dNdx3.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,46) = wq.*(dNdy3.*(D22.*dNdy3.*(kkk + (phi - 1).^2) + D32.*dNdx3.*(kkk + (phi - 1).^2)) + dNdx3.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,47) = wq.*(dNdy3.*(D12.*dNdx4.*(kkk + (phi - 1).^2) + D32.*dNdy4.*(kkk + (phi - 1).^2)) + dNdx3.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,48) = wq.*(dNdy3.*(D22.*dNdy4.*(kkk + (phi - 1).^2) + D32.*dNdx4.*(kkk + (phi - 1).^2)) + dNdx3.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,49) = wq.*(dNdx4.*(D11.*dNdx1.*(kkk + (phi - 1).^2) + D31.*dNdy1.*(kkk + (phi - 1).^2)) + dNdy4.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,50) = wq.*(dNdx4.*(D21.*dNdy1.*(kkk + (phi - 1).^2) + D31.*dNdx1.*(kkk + (phi - 1).^2)) + dNdy4.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,51) = wq.*(dNdx4.*(D11.*dNdx2.*(kkk + (phi - 1).^2) + D31.*dNdy2.*(kkk + (phi - 1).^2)) + dNdy4.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,52) = wq.*(dNdx4.*(D21.*dNdy2.*(kkk + (phi - 1).^2) + D31.*dNdx2.*(kkk + (phi - 1).^2)) + dNdy4.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,53) = wq.*(dNdx4.*(D11.*dNdx3.*(kkk + (phi - 1).^2) + D31.*dNdy3.*(kkk + (phi - 1).^2)) + dNdy4.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,54) = wq.*(dNdx4.*(D21.*dNdy3.*(kkk + (phi - 1).^2) + D31.*dNdx3.*(kkk + (phi - 1).^2)) + dNdy4.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,55) = wq.*(dNdx4.*(D11.*dNdx4.*(kkk + (phi - 1).^2) + D31.*dNdy4.*(kkk + (phi - 1).^2)) + dNdy4.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));
    kpp(:,56) = wq.*(dNdx4.*(D21.*dNdy4.*(kkk + (phi - 1).^2) + D31.*dNdx4.*(kkk + (phi - 1).^2)) + dNdy4.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));
    kpp(:,57) = wq.*(dNdy4.*(D12.*dNdx1.*(kkk + (phi - 1).^2) + D32.*dNdy1.*(kkk + (phi - 1).^2)) + dNdx4.*(D13.*dNdx1.*(kkk + (phi - 1).^2) + D33.*dNdy1.*(kkk + (phi - 1).^2)));
    kpp(:,58) = wq.*(dNdy4.*(D22.*dNdy1.*(kkk + (phi - 1).^2) + D32.*dNdx1.*(kkk + (phi - 1).^2)) + dNdx4.*(D23.*dNdy1.*(kkk + (phi - 1).^2) + D33.*dNdx1.*(kkk + (phi - 1).^2)));
    kpp(:,59) = wq.*(dNdy4.*(D12.*dNdx2.*(kkk + (phi - 1).^2) + D32.*dNdy2.*(kkk + (phi - 1).^2)) + dNdx4.*(D13.*dNdx2.*(kkk + (phi - 1).^2) + D33.*dNdy2.*(kkk + (phi - 1).^2)));
    kpp(:,60) = wq.*(dNdy4.*(D22.*dNdy2.*(kkk + (phi - 1).^2) + D32.*dNdx2.*(kkk + (phi - 1).^2)) + dNdx4.*(D23.*dNdy2.*(kkk + (phi - 1).^2) + D33.*dNdx2.*(kkk + (phi - 1).^2)));
    kpp(:,61) = wq.*(dNdy4.*(D12.*dNdx3.*(kkk + (phi - 1).^2) + D32.*dNdy3.*(kkk + (phi - 1).^2)) + dNdx4.*(D13.*dNdx3.*(kkk + (phi - 1).^2) + D33.*dNdy3.*(kkk + (phi - 1).^2)));
    kpp(:,62) = wq.*(dNdy4.*(D22.*dNdy3.*(kkk + (phi - 1).^2) + D32.*dNdx3.*(kkk + (phi - 1).^2)) + dNdx4.*(D23.*dNdy3.*(kkk + (phi - 1).^2) + D33.*dNdx3.*(kkk + (phi - 1).^2)));
    kpp(:,63) = wq.*(dNdy4.*(D12.*dNdx4.*(kkk + (phi - 1).^2) + D32.*dNdy4.*(kkk + (phi - 1).^2)) + dNdx4.*(D13.*dNdx4.*(kkk + (phi - 1).^2) + D33.*dNdy4.*(kkk + (phi - 1).^2)));;
    kpp(:,64) = wq.*(dNdy4.*(D22.*dNdy4.*(kkk + (phi - 1).^2) + D32.*dNdx4.*(kkk + (phi - 1).^2)) + dNdx4.*(D23.*dNdy4.*(kkk + (phi - 1).^2) + D33.*dNdx4.*(kkk + (phi - 1).^2)));;    
    
    S_all = zeros(numelem.*4.*64,3);

 load ff.mat 'ff'
    for ik = 1 : 64
        i_ = ff(ik,1); 
        j_ = ff(ik,2); 
        S_all(ik:64:end,1) = sctrB(:,i_);
        S_all(ik:64:end,2) = sctrB(:,j_);
        S_all(ik:64:end,3) = kpp(:,ik) ;
    end
%     all at once (memory issue)
%     K_3 = sparse(S_all(:,1),S_all(:,2),S_all(:,3),size(K,1),size(K,2));
    
% divide it into blocks (less memory)
b = 2000000; % block size
n = numelem.*4.*64;
c_ = diff([0:b:n-1,n]) ; 
K_3 = 0*K ; 
for i_step = 1 : size(c_,2)
    i_step/size(c_,2)
    
    s0 = sum(c_(1:i_step-1));
    s_step = c_(i_step);
    s_des = [s0+1:s0+s_step] ;
    K_3 = K_3 + sparse(S_all(s_des,1),S_all(s_des,2),S_all(s_des,3),size(K,1),size(K,2));
end

    K = K_3 ; 
    
end
%         Assemble_staggered  ; 
        KK = K ; 

        
%% External Boundary Conditions

    A = K ; r = RHS ;
for  ui = 1 : length(topNodes) 
    cur_node1 = topNodes ( ui )  ;
    cur_node2 = botNodes  ( ui )  ;
    [A,r] = boundary_1point(A,r,2*cur_node1,-(iNR==1)*u_pres);
    [A,r] = boundary_1point(A,r,2*cur_node2,0);
%     [A,r] = boundary_1point(A,r,2*cur_node2-1,0);
end
    [A,r] = boundary_1point(A,r,2*lrn-1,0);
    

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

pause(0.1)

end


Plot_stress ;


if profile_
    profile viewer
    profsave
end


return

%% symbolic
clear all 
syms dNdx1 dNdx2 dNdx3 dNdx4
syms dNdy1 dNdy2 dNdy3 dNdy4
syms N1 N2 N3 N4
syms sctrB1 sctrB2 sctrB3 sctrB4
sctrB = [sctrB1 sctrB2 sctrB3 sctrB4];
assume([N1 N2 N3 N4]>0)
assume([sctrB1 sctrB2 sctrB3 sctrB4]>0)
syms G l_ wq  phi
syms u1 u2 u3 u4 
sai_0 = 0 ; 
u_sctr_phi =[u1 ; u2; u3 ; u4] 
assume(dNdx1>0)
assume(dNdx2>0)
assume(dNdx3>0)
assume(dNdx4>0)
assume(dNdy1>0)
assume(dNdy2>0)
assume(dNdy3>0)
assume(dNdy4>0)

dNdx = [dNdx1 dNdy1 ; dNdx2 dNdy2 ; dNdx3 dNdy3 ; dNdx4 dNdy4 ];
N = [N1 ; N2 ; N3 ; N4 ] ;
B_ = [ dNdx(:,1)  dNdx(:,2)]';
Kphi_phi = (G/l_*N*N' + G*l_*B_'*B_)*wq; 
fphi = ( B_'*B_*G*l_*u_sctr_phi + N*(phi*(G/l_+2*sai_0)-2*sai_0) )*wq;



[i j s] = find(Kphi_phi); 
Cell = [sctrB(i)' ,sctrB(j)', Kphi_phi(:)];



%% Q4 - Momentum
syms kkk
syms sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 sctrB7 sctrB8
syms D11 D12 D13 D21 D22 D23 D31 D32 D33 
D = [ D11 D12 D13 ; D21 D22 D23 ; D31 D32 D33 ] ;
sctrB = [sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 sctrB7 sctrB8];
assume([sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 sctrB7 sctrB8]>0)
nn = 4  ;
B = sym(zeros(3,8));
B(1,1:2:2*nn)      = dNdx(:,1)';
B(2,2:2:2*nn) = dNdx(:,2)';
B(3,1:2:2*nn)      = dNdx(:,2)';
B(3,2:2:2*nn) = dNdx(:,1)';

Ku_u = ( (1-phi)^2 + kkk )* B'*D*B*wq ; 
% fu   = ( (1-phi)^2 + kkk )* B'*D*B*wq*u(sctrB) ; 
[i j s] = find(Ku_u); 
Cell = [sctrB(i)' ,sctrB(j)', Ku_u(:)];

[sctrB(i)' ,sctrB(j)']
%% T3 - Momentum
clear all 
syms dNdx1 dNdx2 dNdx3 
syms dNdy1 dNdy2 dNdy3 
syms N1 N2 N3 
syms sctrB1 sctrB2 sctrB3 
sctrB = [sctrB1 sctrB2 sctrB3 ];
assume([N1 N2 N3 ]>0)
assume([sctrB1 sctrB2 sctrB3 ]>0)
syms G l_ wq  phi
syms u1 u2 u3  
sai_0 = 0 ; 
u_sctr_phi =[u1 ; u2; u3 ] 
assume(dNdx1>0)
assume(dNdx2>0)
assume(dNdx3>0)
assume(dNdy1>0)
assume(dNdy2>0)
assume(dNdy3>0)

dNdx = [dNdx1 dNdy1 ; dNdx2 dNdy2 ; dNdx3 dNdy3 ];
N = [N1 ; N2 ; N3  ] ;
B_ = [ dNdx(:,1)  dNdx(:,2)]';

syms kkk
syms sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 
syms D11 D12 D13 D21 D22 D23 D31 D32 D33 
D = [ D11 D12 D13 ; D21 D22 D23 ; D31 D32 D33 ] ;
sctrB = [sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 ];
assume([sctrB1 sctrB2 sctrB3 sctrB4 sctrB5 sctrB6 ]>0)
nn = 3  ;
B = sym(zeros(3,6));
B(1,1:2:2*nn)      = dNdx(:,1)';
B(2,2:2:2*nn) = dNdx(:,2)';
B(3,1:2:2*nn)      = dNdx(:,2)';
B(3,2:2:2*nn) = dNdx(:,1)';

Ku_u = ( (1-phi)^2 + kkk )* B'*D*B*wq ; 
% fu   = ( (1-phi)^2 + kkk )* B'*D*B*wq*u(sctrB) ; 
[i j s] = find(Ku_u); 
Cell = [sctrB(i)' ,sctrB(j)', Ku_u(:)];

[sctrB(i)' ,sctrB(j)']

%% generate FEM mesh
! rm -r out_FEM/
! mkdir out_FEM/
close all ; 
for alpha_ = -[0 0.1:0.1:pi/2 pi/2]
    
clearvars -except  alpha_

L = 1 ; 

cd ../gmsh 
eval(['!sed -i ',char(39),'2s/.*/alpha_ = ',num2str(alpha_),...
    '; /',char(39),' fault2d_2.geo']);

gmsh_2 ;
cd ../phase_field

element_app = element  ;
node_app = node ; 

% figure
% plot(node(:,1),node(:,2),'r.')
% axis equal

T3_solve ;


% get crack opening
g = [conint(:,[1 4]) ; conint(end,[2 3])]
xc = zeros(size(g,1),5);
for ii = 1 : size(g,1)
    sc = g(ii,[1 2]) ;
    dx = diff(u(2*sc-1));
    dy = diff(u(2*sc-0));
    dc = sqrt(dx^2+dy^2);
    xx = (node(sc(1),1));
    yy = (node(sc(1),2));
    xc(ii,:) = [xx yy dc dx dy] ;
end
figure
plot(xc(:,1),xc(:,3),'bsq-')

name = ['xc_' num2str(-alpha_) '.mat'];
eval(['save out_FEM/' name]);
pause(1)

end