

%{
    kglob = sparse ( length(node)*2 , length(node)*2 ) ; 


for ielm = 1:nem
    [coordloc,uloc] = locvalues(coord,con,u_elastic,ndf,nne,ielm,lin_quad,0);    % read values (coordinates and displacements) of the finite element nodes
    sctr = con(ielm,:);
    coordloc = node(sctr,:) ;
    
    [bmat,xyIP_b,detJ] = bmatrix(rsIP_b,coordloc,lin_quad);                 % compatibility matrix (strain - displacements) for each integration point
    strain = multiprod(bmat,uloc);                                          % strain field for each integration point
    stress = multiprod(cmat_b,strain);                                      % stress field for each integration point
    [kmat] = kmatrix(bmat,cmat_b,detJ,b_b);                                 % local stiffness matrix

    % positions of local stiffness matrix in vector of MATLAB's sparse matrix format
    vectorID = ((ielm-1)*(nne*ndf)^2) + 1;
    kglob_value(vectorID:vectorID+(ndf*nne)^2-1,1) = kmat(:);

    sctrr = con(ielm,1:end); 
    dof = [sctrr(1)*2-1 sctrr(1)*2 sctrr(2)*2-1 sctrr(2)*2   sctrr(3)*2-1 sctrr(3)*2   sctrr(4)*2-1 sctrr(4)*2  ];
    kglob(dof,dof ) = kglob(dof,dof) + kmat;

    % internal force field (out-of-balance)
    [fint] = fmatrix(bmat,stress,detJ,b_b);             
    [fintglob] = assemblefintglob(fintglob,fint,ndf,con,ielm,lin_quad,0);

    [sss] = fmatrix_strain(bmat,uloc,detJ,cmat_b); 

end
        



%%
kglob = sparse ( length(node)*2 , length(node)*2 ) ; 

% Bulk assemble
volume_ = 0 ; 
for iel = 1 : nem 
%     iel/numelem 
    sctr = con(iel,:); % element connectivity
%     sctr_app = element_app(iel,:); % element connectivity
%     nn   = length(sctr_app);   % number of nodes per element

%     sctr_app = element_app(iel,:); % element connectivity
%     idx = find(sctr_app~=0);
%     sctr_app(sctr_app==0)=[] ;
    nn   = length(sctr);   % number of nodes per element

    % scatter vector for element assembly
%     sctrB = [sctr_app sctr_app+numnode];
    sctrB = zeros(1,length(sctr)*2) ;
    sctrB(1:2:end) = 2*sctr-1 ; 
    sctrB(2:2:end) = 2*sctr ;
    
    B = zeros(3,2*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    KK_ = 0 ; 
    
        [W,Q] = quadrature(2,'GAUSS',2);

    
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        [N,dNdxi] = lagrange_basis('Q4',pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
        
        % Stiffness matrix = B^T * D * B * detJ * w
        KK_ = KK_ + B'*cmat_b*B*W(kk)*det(J0);
        
        volume_ = volume_ +  W(kk)*det(J0) ; 
    end                  % end of looping on GPs
    kglob(sctrB,sctrB) = kglob(sctrB,sctrB) + KK_;

end                      % end of looping on elements

volume_

%}

%%
% tic

kglob = sparse ( length(node)*2 , length(node)*2 ) ; 
mglob = sparse ( length(node)*2 , length(node)*2 ) ; 
C = cell(nem,1);
C2 = cell(nem,1);
stepF = 1000 ;
Fb = sparse ( length(node)*2 , 1 ) ; 
g = 1*9.806; 
% Bulk assemble
volume_ = 0 ; 
for iel = 1 : nem 
    if rem(iel,1000)==0
    iel/nem
    end
    
    sctr = con(iel,:); % element connectivity

    nn   = length(sctr);   % number of nodes per element

    % scatter vector for element assembly
%     sctrB = [sctr_app sctr_app+numnode];
    sctrB = zeros(1,length(sctr)*2) ;
    sctrB(1:2:end) = 2*sctr-1 ; 
    sctrB(2:2:end) = 2*sctr ;
    
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    B = zeros(3,2*nn);
    Nu = zeros(2,2*nn);
    KK_ = 0 ; NN = 0 ; Nb = 0; 
    [W,Q] = quadrature(4,'GAUSS',2);    
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        [N,dNdxi] = Lag_basis('Q4',pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
        Nu(1,1:2:2*nn) = N; 
        Nu(2,2:2:2*nn) = N; 
        
        
%         get phase field parameter at this gp
        phi = N'*PHI(sctr) ;
%         Ku_u = ( (1-phi)^2 + kkk )* B'*cmat_b*B*W(kk)*det(J0) ; 
%         fu   = ( (1-phi)^2 + kkk )* B'*D*B*W(kk)*det(J0)*u(sctrB) ; 

        
        
        % Stiffness matrix = B^T * D * B * detJ * w
        KK_ = KK_ + ( (1-phi)^2 + kkk )*B'*cmat_b*B*W(kk)*det(J0);
        NN = NN + Nu'*rho*Nu*W(kk)*det(J0);
        Nb = Nb + Nu'*rho*[0;g]*W(kk)*det(J0);

        volume_ = volume_ +  W(kk)*det(J0) ; 
    end  
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK_;
    Fb(sctrB,1) = Fb(sctrB,1) + Nb;

        [i j s] = find(KK_); 
        C{iel} = [sctrB(i)' ,sctrB(j)', s];
    
        [i2 j2 s2] = find(NN); 
        C2{iel} = [sctrB(i2)' ,sctrB(j2)', s2];
        
if rem(iel,stepF)==0
IJV = cell2mat( C );
A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(kglob,1),size(kglob,1));
kglob = kglob + A;
C = cell(stepF,1);
c = 0 ; 

IJV = cell2mat( C2 );
A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(mglob,1),size(mglob,1));
mglob = mglob + A;
C2 = cell(stepF,1);
c2 = 0 ; 
end
    
end                      % end of looping on elements

IJV = cell2mat( C );
kglob = kglob + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(kglob,1),size(kglob,1));

IJV = cell2mat( C2 );
mglob = mglob + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(mglob,1),size(mglob,1));

% toc
% KKKKKKK = kglob ; 