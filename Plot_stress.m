%% Compute stress at gps
% %{
if P_order > 1 
stress = zeros(size(W,1)*numelem,6) ; 
QQQ = zeros(size(W,1)*numelem,2) ;
cc = 1 ; 
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
    sctr_phi = sctr_app ; 

    U = zeros(2*length(sctr_app),1) ; 
    U(1:2:end) = u(2*sctr_app-1) ; 
    U(2:2:end) = u(2*sctr_app) ; 
        
    
    B = zeros(3,2*nn);
    % ---------------------
    % Loop on Gauss points   
    % ---------------------
    KK = 0 ; 
    
    
%     if ismember(iel,desir_elems)
    if length(sctr_app)>4
        % Gauss quadrature
        order = n_order ;
    else
        order = 2 ;
    end    
        [W,Q] = quadrature_v2(order,'GAUSS',2);
        

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
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
        
        strain = B*U;
        
%         get phase field parameter at this gp
        phi = N'*PHI(sctr_phi) ;
        dispx = N'* u(2*sctr_app-1) ;
        dispy = N'* u(2*sctr_app-0) ;
        % Then stress from strain via constitutive equation
        stress(cc,:) = [ [D*strain]' dispx dispy phi] ;
        QQQ(cc,:) = [xgp] ;
        cc= cc + 1 ; 
    end                  % end of looping on GPs
%     K(sctrB,sctrB) = K(sctrB,sctrB) + KK;

end                      % end of looping on elements

% volume_


tri = delaunay(QQQ(:,1),QQQ(:,2));
% tri = tricheck(QQQ(:,1:2),tri) ; 

% VTKPostProcess(QQQ,tri,1,'Tri3','stress',stress(:,1:3),1*[stress(:,4:6) ])
% !paraview stress.vtu&


fac = 10;
QQQ_d = QQQ + fac*[stress(:,4) stress(:,5) ];
figure
plot(QQQ_d(:,1),QQQ_d(:,2),'r.')
axis equal


%%
! rm -r out/
! mkdir out/

xx = QQQ ;    
if size(xx,2)==2
    xx(:,3) = 0;
end 
yy = tri ;        
VTU_ (xx,yy,stress,'out/all.vtu'); 
!paraview out/all.vtu&


end
%}

%% get values along a line
% ! rm -r out/
! mkdir out/

nq_ = 50 ; 
Q_0 = [-1:2/nq_:1];
% Q_ = [Q_'*0 Q_'];
% Q_ = [Q_'*0 Q_'];
QQQ_ = zeros(length(elem_center)*(nq_+1),5) ;
iii = 1 ; 
figure
hold on

for lj = 0% : 1 : nnx/2
    
for kj = 0%[-1:0.2:1]

Q_ = [Q_0'*0+kj Q_0'];
QQQ_ = [] ; 
    for iel_ = 1 : length(elem_center) 
        if rem(iel,10)==0
            iel_/length(elem_center) 
        end
        iel = elem_center(iel_)-lj;
        sctr = element(iel,:); % element connectivity

        sctr_app = element_app(iel,:); % element connectivity
        idx = find(sctr_app~=0);
        sctr_app(sctr_app==0)=[] ;
        nn   = length(sctr_app);   % number of nodes per element

        % scatter vector for element assembly
        sctr_phi = 0*numnode+sctr_app ; 

        B = zeros(3,2*nn);
        % ---------------------
        % Loop on Gauss points   
        % ---------------------
        KK = 0 ; 


        for kk = 1 : size(Q_,1)
            pt = Q_(kk,:);                             % quadrature point


            % Shape functions and its derivatives
            [N,dNdxi] = lobatto_basis_v2(elemType,pt,1,P_order);
            J0 = node(sctr,:)'*dNdxi(1:4,:);                 % element Jacobian matrix
            N=N(idx);
            dNdxi=dNdxi(idx,:); 
            xgp = N(1:4)'*node(sctr,:) ; 

    %         invJ0 = inv(J0);
    %         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY

    %         

    %         get phase field parameter at this gp
            phi = N'*PHI(sctr_phi) ;
            dispx = N'* u(2*sctr_app-1) ;
            dispy = N'* u(2*sctr_app-0) ;
            % Then stress from strain via constitutive equation
            QQQ_(iii,:) = [xgp phi dispx dispy] ;
            iii = iii + 1; 
        end                  % end of looping on GPs

    end

plot(QQQ_(:,2),QQQ_(:,3),'.-')
pause(0.01)
end
% pause

end
figure
plot(QQQ_(:,2),QQQ_(:,3),'.-')


figure
plot(QQQ_(:,2),QQQ_(:,5),'b.-')
hold on
plot(QQQ_(:,2),QQQ_(:,3)/100,'r.-')
plot(QQQ_(:,2),QQQ_(:,4),'.-')




%%

load gg.mat 'ff'
plot(ff(:,10),ff(:,1),'g-')
plot(ff(:,10),ff(:,2),'g-')

if vec_ == 1 
       xx = N1.*node(sctrB_(:,1),1) + N2.*node(sctrB_(:,2),1) + N3.*node(sctrB_(:,3),1) + N4.*node(sctrB_(:,4),1);
       yy = N1.*node(sctrB_(:,1),2) + N2.*node(sctrB_(:,2),2) + N3.*node(sctrB_(:,3),2) + N4.*node(sctrB_(:,4),2);
       ux = N1.*u(2*sctrB_(:,1)-1) + N2.*u(2*sctrB_(:,2)-1) + N3.*u(2*sctrB_(:,3)-1) + N4.*u(2*sctrB_(:,4)-1);
       uy = N1.*u(2*sctrB_(:,1)-0) + N2.*u(2*sctrB_(:,2)-0) + N3.*u(2*sctrB_(:,3)-0) + N4.*u(2*sctrB_(:,4)-0);
       phi = N1.*PHI(sctrB_(:,1)) + N2.*PHI(sctrB_(:,2)) + N3.*PHI(sctrB_(:,3)) + N4.*PHI(sctrB_(:,4));
       QQQ = [xx yy 0*xx];
       tri = delaunay(QQQ(:,1),QQQ(:,2));
%         tri = tricheck(QQQ(:,1:2),tri) ; 
%        stress = [
        yy = tri ;        
        VTU_2 (QQQ,yy,ux,'out/ux.vtu'); 
        VTU_2 (QQQ,yy,uy,'out/uy.vtu'); 
        VTU_2 (QQQ,yy,phi,'out/phi.vtu'); 
% !paraview out/all.vtu&

end
%%
fac = 10;
QQQ_d = QQQ + fac*[ux uy 0*uy];
figure
plot(QQQ_d(:,1),QQQ_d(:,2),'r.')
axis equal
%% compute crack opening
% method 1
% %{
d_crack = zeros(size(crack_nodes2,1),4);
for ip = 1 : size(crack_nodes2,1)
    if rem(ip,100)==0
        ip/size(crack_nodes2,1)
    end
    
    xx = crack_nodes2(ip,:);
    
    
%% find parent element and local coordinate
    local = crack_local2(ip,:);
    sctr = element(crack_elems_2(ip),:) ; 
    sctr_app = element_app(crack_elems_2(ip),:); % element connectivity
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
    
    ux = N'*u(sctr_app*2-1);
    uy = N'*u(sctr_app*2-0);
    d_crack(ip,:) = [xx ux uy];
end

opening_x = d_crack(1:n_Cr,3)-d_crack(n_Cr+1:end,3);

opening_y = d_crack(1:n_Cr,4)-d_crack(n_Cr+1:end,4);

opening = sqrt(opening_x.^2+opening_y.^2);


figure
hold on
plot(d_crack(1:n_Cr,1),opening,'c.')

%}

% or 
d_crack = zeros(size(crack_nodes,1),8);
for ip = 1 : size(crack_nodes,1)
    if rem(ip,100)==0
        ip/size(crack_nodes,1)
    end
    
    xx = crack_nodes(ip,:);
    
    
%% find parent element and local coordinate
    local = crack_local(ip,:);
    local_t = [local(1) 1];
    local_b = [local(1) -1];
%     local_t = [1 local(2) ];
%     local_b = [-1 local(2) ];

    sctr = element(crack_elems_(ip),:) ;     sctr_app = element_app(crack_elems_(ip),:);     idx = find(sctr_app~=0);    sctr_app(sctr_app==0)=[] ;
    
% top
    [N,dNdxi] = lobatto_basis_v2(elemType,local_t,1,P_order);  % element shape functions
    J0 = node(sctr,:)'*dNdxi(1:4,:);                  N=N(idx);    dNdxi=dNdxi(idx,:);     invJ0 = inv(J0);  dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
    xx_t = N(1:4)'*node(sctr,:);
    ux_t = N'*u(sctr_app*2-1);
    uy_t = N'*u(sctr_app*2-0);

% bot
    [N,dNdxi] = lobatto_basis_v2(elemType,local_b,1,P_order);  % element shape functions
    J0 = node(sctr,:)'*dNdxi(1:4,:);                  N=N(idx);    dNdxi=dNdxi(idx,:);     invJ0 = inv(J0);  dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
    xx_b = N(1:4)'*node(sctr,:);
    ux_b = N'*u(sctr_app*2-1);
    uy_b = N'*u(sctr_app*2-0);

    
    
    d_crack(ip,:) = [xx_t ux_t uy_t xx_b ux_b uy_b ];
end

opening_x = d_crack(1:n_Cr,3)-d_crack(1:n_Cr:end,7);
opening_y = d_crack(1:n_Cr,4)-d_crack(1:n_Cr:end,8);
opening = sqrt(opening_x.^2+opening_y.^2);



% figure
hold on
plot(d_crack(1:n_Cr,1),opening_y,'b.')


name = ['xc_' num2str(-alpha) '.mat'];
eval(['load out_FEM/' name ' xc']);
plot(xc(:,1),xc(:,3),'rsq')


