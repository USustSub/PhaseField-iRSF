method = 'Lagrange'; % 'Lobatto' or 'Lagrange'
Const_method = 'Penalty' ; 
Penalty_ = 1 ; 



% Dimension of the domain 
D = 1 ; 
L = 1 ;
% Four corner points
pt1 = [0 0] ;
pt2 = [D 0] ;
pt3 = [D L] ;
pt4 = [0 L] ;

elemType = 'T3' ; 

% Number of nodes along two directions
% nnx = 50 ;
% nny = 50 ;

% Material properties
E  = 210000 ;
nu = 0.3 ;
compliance_ ;
G = 2.7 ; 
u_pres = 0.01 ; 
NR.iter      = 30;
NR.tolerance = 1e-10;

numnode = size(node_app,1);
numnode0 = size(node,1);
numelem = size(element,1);


theta = alpha ; 
a = L/4; % crack length

% crack geometry
xCr = [ L/2-a*cos(theta) L/2-a*sin(theta) ; L/2+a*cos(theta) L/2+a*sin(theta)];
xCr(:,2) = xCr(:,2)  ; % + .0012 ; 

% define essential boundaries
botNodes = find(node(:,2) == 0); 
topNodes = find(node(:,2) == L); 
rightNodes = find(node(:,1) == L); 
leftNodes = find(node(:,1) == 0);

lrn = find(node(:,2)==0 & node(:,1)==L);

dispNodes = [leftNodes];
dispNodes = unique(dispNodes);

figure
hold on
plot(node(:,1),node(:,2),'r.')
plot(xCr(:,1),xCr(:,2),'b-','LineWidth',2)
plot(node(topNodes,1),node(topNodes,2),'bsq')
plot(node(botNodes,1),node(botNodes,2),'bsq')
plot(node(lrn,1),node(lrn,2),'csq')
axis equal

%%
    nDofs = 2*numnode  ; 
    
    K   = sparse( nDofs , nDofs);
    RHS = sparse( nDofs, 1 );   
    volume_ = 0 ; 
for iel = 1 : size(element,1) 
    if rem(iel,2000)==0
        iel/size(element,1) 
    end
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
        [W,Q] = quadrature_v2( 2, 'TRIANGULAR', 2 );
        
    
    B = zeros(3,2*nn);
    B_ = zeros(2,1*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature_v2 point
        
        % Shape functions and its derivatives
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('T3',pt);  % element shape functions
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
        
        xgp = N(1:3)'*node(sctr,:) ;  

        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
         
        B_(1,1:1:1*nn)      = dNdx(:,1)';
        B_(2,1:1:1*nn) = dNdx(:,2)';
        
        
        Ku_u =  B'*D*B*W(kk)*det(J0) ; 
%         fu   =  B'*D*B*W(kk)*det(J0)*u(sctrB) ; 
        
        
        
% get phase field at current node        
        K(sctrB,sctrB) = K(sctrB,sctrB) + Ku_u;
%         RHS(sctrB) = RHS(sctrB) + fu ; 
                
        volume_ = volume_ +  W(kk)*det(J0) ; 
    end                  % end of looping on GPs
end                      % end of looping on elements



    A = K ; r = RHS ;
for  ui = 1 : length(topNodes) 
    cur_node1 = topNodes ( ui )  ;
    cur_node2 = botNodes  ( ui )  ;
    [A,r] = boundary_1point(A,r,2*cur_node1,-1*u_pres);
    [A,r] = boundary_1point(A,r,2*cur_node2,0);
%     [A,r] = boundary_1point(A,r,2*cur_node2-1,0);
end
    [A,r] = boundary_1point(A,r,2*lrn-1,0);
    

        u = -A\r;
                

%% Compute stress at gps

stress = zeros(numelem*length(W),6) ; 
QQQ = [] ;
cc =  1; 
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

    % scatter vector for element assembly
%     sctrB = [sctr_app sctr_app+numnode];
    sctrB = zeros(1,length(sctr_app)*2) ;
    sctrB(1:2:end) = 2*sctr_app-1 ; 
    sctrB(2:2:end) = 2*sctr_app ;
    sctr_phi = 2*numnode+sctr_app ; 

    U = zeros(2*length(sctr_app),1) ; 
    U(1:2:end) = u(2*sctr_app-1) ; 
    U(2:2:end) = u(2*sctr_app) ; 
        
    
    B = zeros(3,2*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    KK = 0 ; 
    

    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        if strcmp(method,'Lagrange')
            [N,dNdxi] = lagrange_basis('T3',pt);  % element shape functions
            J0 = node_app(sctr_app,:)'*dNdxi(:,:);                 % element Jacobian matrix
        elseif strcmp(method,'Lobatto')
            [N,dNdxi] = lobatto_basis(elemType,pt,1,P_order);  % element shape functions
            J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
            N=N(idx);
            dNdxi=dNdxi(idx,:); 
        end
        xgp = N(1:3)'*node(sctr,:) ; 
        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
        
        strain = B*U;
        
%         get phase field parameter at this gp
%         phi = N'*u(sctr_phi) ;
        dispx = N'* u(2*sctr_app-1) ;
        dispy = N'* u(2*sctr_app-0) ;
        % Then stress from strain via constitutive equation
        stress(cc,:) = [[D*strain]' dispx dispy 0] ;
        QQQ = [QQQ  ; xgp] ;
        cc = cc + 1 ; 
    end                  % end of looping on GPs
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK;

end                      % end of looping on elements

% volume_


tri = delaunay(QQQ(:,1),QQQ(:,2));
% tri = tricheck(QQQ(:,1:2),tri) ; 



% VTKPostProcess(QQQ,tri,1,'Tri3','stress',stress(:,1:3),1*[stress(:,4:6) ])
output = ['out_FEM/FEM_ref' num2str(theta) '.vtu'] ; 


xx = QQQ ;    
if size(xx,2)==2
    xx(:,3) = 0;
end 
yy = tri ;        
VTU_ (xx,yy,stress,output); 

