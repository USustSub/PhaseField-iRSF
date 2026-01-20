%% Compute stress at gps
%{
stress = zeros(numelem*4,6) ; 
QQQ = zeros(numelem*4,2) ;
for iel = 1 : numelem 
    if rem(iel,100)==0
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
%         dispx = N'* u(2*sctr_app-1) ;
%         dispy = N'* u(2*sctr_app-0) ;
        % Then stress from strain via constitutive equation
        stress = [ stress;  [D*strain]' 0 0 phi] ;
        QQQ = [QQQ  ; xgp] ;
    end                  % end of looping on GPs
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK;

end                      % end of looping on elements

% volume_
for ip = 1 : size(crack_nodes,1)
    if rem(ip,1)==0
        ip/size(crack_nodes,1)
    end
    
    xx = crack_nodes(ip,:);
    
    
%% find parent element and local coordinate
    local = crack_local(ip,:);
    sctr = element(crack_elems_(ip),:) ; 
    sctr_app = element_app(crack_elems_(ip),:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;
        
    xgp = N(1:4)'*node(sctr,:) ; 

%     [N,~] = lobatto_basis(elemType,local,1,P_order);  % element shape functions
    [N,~] = lobatto_basis_v2(elemType,local,1,P_order);  % element shape functions
    
    if norm(xx-N(1:4)'*node(sctr,:)) > 1e-10
        error('errr')
    end
    
%     Penalty constraint
    stress = [ stress;  [D*strain]' 0 0 N'*u(sctr_app)] ;
    QQQ = [QQQ  ; xgp] ;

        
    
end

tri = delaunay(QQQ(:,1),QQQ(:,2));
tri = tricheck(QQQ(:,1:2),tri) ;
%%
! rm -r out/
! mkdir out/

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


% !paraview out/all.vtu&

 %}

%%

nq_ = 100 ; 
Q_ = [-1:2/nq_:1];
Q_ = [Q_'*0 Q_'];
QQQ_ = zeros(length(elem_center)*(nq_+1),3) ;
iii = 1 ; 
for iel_ = 1 : length(elem_center) 
    if rem(iel,10)==0
        iel_/length(elem_center) 
    end
    iel = elem_center(iel_);
    sctr = element(iel,:); % element connectivity

    sctr_app = element_app(iel,:); % element connectivity
    idx = find(sctr_app~=0);
    sctr_app(sctr_app==0)=[] ;
    nn   = length(sctr_app);   % number of nodes per element

    
    
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
    

    for kk = 1 : size(Q_,1)
        pt = Q_(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('Q9',pt);  % element shape functions
            J0 = node_app(sctr_app,:)'*dNdxi(:,:);                 % element Jacobian matrix
        elseif strcmp(method,'Lobatto')
%             if order == 2 
%                 [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,1);
%             else
                [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
%                 N = N_p{kk}{1};
%                 dNdxi = N_p{kk}{2};
%             end
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
        
%         get phase field parameter at this gp
        phi = N'*u(sctr_phi) ;
%         dispx = N'* u(2*sctr_app-1) ;
%         dispy = N'* u(2*sctr_app-0) ;
        % Then stress from strain via constitutive equation
        QQQ_(iii,:) = [xgp phi] ;
        iii = iii + 1; 
    end                  % end of looping on GPs
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK;

end

figure
plot(QQQ_(:,2),QQQ_(:,3),'.-')