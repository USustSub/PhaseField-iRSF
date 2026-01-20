%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Continuum model for faults with RSF friction - phase field
% By: Mohsen Goudarzi (2023, Utrecht), Goudarzi.mohsen@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear all ; 
    close all ;
    clc 

    addpath '/home/mohsen/Desktop/MatlabCodes/Interface Cohesive Crack/Test2'
    addpath ('PFEM_functions')

% ! rm -r stress
% ! mkdir stress
! rm -r out
! mkdir out
! mkdir out/paraview
! mkdir out/slip
! mkdir out/slip_n
! mkdir out/slip_rate
! mkdir out/theta
! mkdir out/mu
! mkdir out/deformed

    ifprofile = 0 ;
    if ifprofile
        profile on
    end

    global volume_ LINE_
    addpath ( genpath ( '/home/mohsen/Desktop/MatlabCodes/MASONRY_jointsasinterfaceelements' ) ) ;
    warning('off', 'MATLAB:nargchk:deprecated')
    
     dynamic = 1 ; 
    msh_type = 'Q4' ; 
    
    Time = 0 ;
    print_step1 = 50 ; 
    print_step2 = 200000000 ;
    
    nn = 540  ;
    dt = 1e6 ;
    h__ = 1 ; 
    disp(['Time step is ' num2str(dt/60/60/24/365) ' years'])
%     50/(dt/60/60/24/365)
    
    if strcmp(msh_type,'Q4')
        Input_ ; 
    elseif strcmp(msh_type,'T3')
        Input_T3 ; 
    end
    tol = 1e-8 ; 
    maxVp = 0 ; 
    
    h2 = figure ; 
    hold on
    plot(node(:,1),node(:,2),'bsq')
    axis equal
    cr = plot(xCr(:,1),xCr(:,2),'k-','LineWidth' , 2);
for ii = 1 : length(conint) 
    node_conint = node(conint(ii,:),:);
    plot(node_conint(:,1),node_conint(:,2),'rsq')
%     pause
end



%% RSF inputs
mu0 = 0.2 ; 
V0 = 4e-9 ;  % m/s
a_ = 0.01 ;
L_ = 0.01 ;
b1_ = 0.001 ; 
b2_ = 0.017 ; 
theta_ = L_/V0*exp(-1);
theta_0 = L_/V0*exp(-1);
LL1 = 32000 ;
LL2 = 108000 ;
ll = norm(xCr(2,:)-xCr(1,:)) ;
% LL1 = xCr(1,1)+0.2*ll;
% LL2 = xCr(2,1)-0.2*ll;
plot([LL1 LL1 ] , [0 L],'k-','LineWidth',4);
plot([LL2 LL2] , [0 L],'k-','LineWidth',4);

kkk = 1e-11 ; % small value

pause(1)
%%
% h1 = figure ; 
crack_grow = 0 ; 
mul = [  ] ;
landa = 1 ; 
conint_cohesive = [ conint ] ;
% start time iterations
istep = 1 ; 
uy = 2e-9 ;  


findBoundaries ;

numnode = length(node) ; 
dispNodes1 = botNodes;
dispNodes2 = topNodes; 
dispNodes3 = [botNodes;topNodes;rightNodes;leftNodes] ; 

anode = 1:2*(numnode);    % index vector of all DOF's
dispNodes = [dispNodes1;dispNodes2;dispNodes3];
presDofs = [dispNodes1*2-1 ;dispNodes2*2-1 ; dispNodes3*2 ];
    
nnm = length(node) ; 
nem = length(element) ; 
nim = size(conint_cohesive ,1) ; 
con = element ; 
coord = node ; 

%% apply initial pressure
DOF = length(node)*2 ; 
Fext = zeros(DOF,1); 
%{

leftNodes = find(node(:,1)==0); leftNodes([length(leftNodes)-1 length(leftNodes)]) = leftNodes([length(leftNodes) length(leftNodes)-1]);
rightNodes = find(node(:,1)==L); rightNodes([length(rightNodes)-1 length(rightNodes)]) = rightNodes([length(rightNodes) length(rightNodes)-1]);
botNodes = find(node(:,2)==0);
topNodes = find(node(:,2)==D);
[~,iii]=sort(node(leftNodes,2));
leftNodes = leftNodes(iii);
[~,iii]=sort(node(rightNodes,2));
rightNodes = rightNodes(iii);
leftEdge = []; rightEdge = [] ; topEdge = [] ; botEdge = [] ; 
for ii = 1 : size(rightNodes,1)-1
    rightEdge = [rightEdge  ; rightNodes(ii) rightNodes(ii+1) ] ;
end
for ii = 1 : size(leftNodes,1)-1
    leftEdge = [leftEdge  ; leftNodes(ii) leftNodes(ii+1) ] ;
end
for ii = 1 : size(topNodes,1)-1
    topEdge = [topEdge  ; topNodes(ii) topNodes(ii+1) ] ;
end
for ii = 1 : size(botNodes,1)-1
    botEdge = [botEdge  ; botNodes(ii) botNodes(ii+1) ] ;
end

% left and right edges
sigmato = 5e6;
[W__,Q__]=quadrature(2,'GAUSS',1);
L_ = 0 ;
for e_ = 1 : size(topEdge,1)
    sctr = topEdge(e_,:);
    sctry = sctr.*2 ;
    for q=1:size(W__,1)
        pt = Q__(q,:);
        wt = W__(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2),1)-node(sctr(1),1) )/2;
        Fext(sctry)=Fext(sctry) - N*sigmato*det(J0)*wt;
        L_ = L_ + det(J0)*wt ; 
    end   % of quadrature loop
end       % of element loop
% L_
%-------------------------------------
L_ = 0 ;
for e_ = 1 : size(botEdge,1)
    sctr = botEdge(e_,:);
    sctry = sctr.*2 ;
    for q=1:size(W__,1)
        pt = Q__(q,:);
        wt = W__(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2),1)-node(sctr(1),1) )/2;
        Fext(sctry)=Fext(sctry) + N*sigmato*det(J0)*wt;
        L_ = L_ + det(J0)*wt ; 
    end   % of quadrature loop
end       % of element loop
% L_
%-------------------------------------
L_ = 0 ;
for e_ = 1 : size(leftEdge,1)
    sctr = leftEdge(e_,:);
    sctrx = sctr.*2-1 ;
    l = 0 ;
    for q=1:size(W__,1)
        pt = Q__(q,:);
        wt = W__(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2),2)-node(sctr(1),2) )/2;
        Fext(sctrx)=Fext(sctrx) + N*sigmato*det(J0)*wt;
        L_ = L_ + det(J0)*wt ; 
        l = l + det(J0)*wt;
    end   % of quadrature loop
end       % of element loop
% L_
%-------------------------------------
L_ = 0 ;
for e_ = 1 : size(rightEdge,1)
    sctr = rightEdge(e_,:);
    sctrx = sctr.*2-1 ;
    l = 0 ;
    for q=1:size(W__,1)
        pt = Q__(q,:);
        wt = W__(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2),2)-node(sctr(1),2) )/2;
        Fext(sctrx)=Fext(sctrx) - N*sigmato*det(J0)*wt;
        L_ = L_ + det(J0)*wt ; 
        l = l + det(J0)*wt;
    end   % of quadrature loop
end       % of element loop
% L_
%}
Time0 = 0 ; 


%% solve for phase field parameter
Solve_for_phase_field

%% time stepping
% load all.mat 
DATA = zeros(nim,14);
f_int_ = zeros(2*length(node),1);
for itime = 1 :  1500

%% update dt
 update_dt
%%
    Time = Time0 + dt ;
    if rem(itime,1)==0
       disp(['Time is --> ' num2str(Time/60/60/24/365) ' years'])
    end
    
    u_elastic    = sparse(DOF,1);
    vu = 0*u_elastic ;
    au = 0*u_elastic ; 
    if itime == 1
        u_elastic_0 = sparse(DOF,1);        
        au0 = sparse(DOF,1);
        vu0 = sparse(DOF,1);
    end
    
% displacement applied
    ap1 =  1*( uy + a3*vu0(dispNodes2*2-1) + a5*au0(dispNodes2*2-1) ) / a1 ;
    ap2 =  1*(-uy + a3*vu0(dispNodes1*2-1) + a5*au0(dispNodes1*2-1) ) / a1 ;
    
    if itime > 1
%         u_elastic(1:length(u_elastic_0))=u_elastic_0  ;
        u_elastic=u_elastic_0  ;
        vu = vu0 ; 
        au = au0 ; 
    end
    
    fintglob = zeros(2*nnm,1);    % internal force vector

    xx = [] ;     jjj = 1; 
    volume_ = 0  ; StrainEnergy =0 ; 
    normm = 110 ; itr = 0  ; 
%% 
while normm>tol && itr<10
        itr = itr + 1; 
%% Bulk Integration 
    if (itr == 1) && ( itime == 1 ) 
        if strcmp(msh_type,'Q4')
            Quad_Matrix_Integral_mass
        elseif strcmp(msh_type,'T3')
            Tri_Matrix_Integral_mass
        end
        kglob00 = kglob ; 
    else
        kglob = kglob00; 
    end
% kglob = 0*kglob 
    fintglob = kglob*u_elastic + mglob*a0*u_elastic; 
    fintglob = fintglob(1:DOF);
    
%% dummy solution
    freeDOF = 1:DOF ;    globalF= sparse(DOF,1);  sigmato= 1; 

%% Integration along the interface
LINE_ = 0 ;
jjj = 1; 
Points = [] ;   alpha = 0;  xx =[ ];    uu = [] ; 
%
MaxP = 0 ; 
theta_p = zeros(nim*2,2) ;   cc = 1 ; 
slip_rate = zeros(nim*2,2) ; 
slip_ = zeros(nim*2,3) ; 
mu_ = zeros(nim*2,2) ; 
VP = zeros(nim*2,1) ; 
% alpha = -atan2(xCr(2,2)-xCr(1,2),xCr(2,1)-xCr(1,1));

seg = xCr(2,:)-xCr(1,:);
alpha = atan(seg(2)/seg(1));

% kk = kglob ; 
% old version
%{
for iint = 1 : nim/1
        [coordloc,uloc] = locvalues(coord,conint_cohesive,u_elastic,ndf,nni,iint,lin_quad,1); % read values (coordinates and displacements) of the finite element nodes
        [~,vloc] = locvalues(coord,conint_cohesive,vu,ndf,nni,iint,lin_quad,1); % read values (coordinates and displacements) of the finite element nodes
        [nmat,xyIP_i,detJ] = nmatrix(rsIP_i,coordloc,lin_quad);                 % compatibility matrix (jump - displacements) for each integration point

%         Points = [ Points ; xyIP_i ] ;
        
        rotmat = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha) ] ;
%          rotmat = Get3DRotationMatrix ( [xCr(1,:) 0] ,  [xCr(2,:) 0], 0 ) 
%         rotmat = [-sin(alpha) -cos(alpha) ; cos(alpha) -sin(alpha) ]; 

        [jump] = jumpmatrix(nmat,rotmat,uloc); 
        [jump_v] = jumpmatrix(nmat,rotmat,vloc);                                   % jump field (normal displacement jump / shear displacement jump)
        [jump_xx] = jumpmatrix(nmat,rotmat,uloc);
        
        
% get data along the interface  
        xx(jjj:jjj+(size(xyIP_i,1)-1),1) = [xyIP_i(:,1)]; 
        xx(jjj:jjj+(size(xyIP_i,1)-1),2) = [xyIP_i(:,2)]; 

        uu(jjj:jjj+(size(xyIP_i,1)-1),1) = jump_xx(1,1,:);
        uu(jjj:jjj+(size(xyIP_i,1)-1),2) = jump_xx(2,1,:);
        
        jjj = jjj + (size(xyIP_i,1)) ; 

%             tmat_secant
        wt =  -(jump(1,:));
        wt_v =  -(jump_v(1,:));
        wn =  -(jump(2,:));
%         Opening = [Opening ;  xyIP_i wt']; 

        tmat_secant = zeros(2,2,size(rsIP_i,1));
        tmat_tangent = zeros(2,2,size(rsIP_i,1));
        traction = zeros(2,1,size(rsIP_i,1)) ;
for IP = 1 : size(rsIP_i,1)
        opening_t = wt(IP) ;
        opening_n = wn(IP) ;
        opening_v_t = wt_v(IP) ;
%         opening_v_t
        parameter_t = 0 ; 
        parameter_n = 0 ; 
        
%  Linear Cohesive Law
        KK = +1*1e8*1 ; 
%         mu = 0.2*5e6 ; 
% % %         if opening_t > wc  || opening_t <= 0 
% % %             parameter_t  = 0 ; 
% % %             dparameter_t  = 0 ; 
% % %         elseif opening_t > 0 && opening_t < wc 
% % %             parameter_t  = ft * ( 1 - opening_t/wc );
% % %             dparameter_t  = -ft/wc ; 
% % %         end
% % %         parameter_t = 0.0001; ; 
        
        parameter_n = KK*opening_n;
        dparameter_n = KK;

        parameter_t  = KK*opening_t;
        dparameter_t  = KK ; 
        
        Vp = opening_v_t;
        if itr == 1
%             Vp = 1e-17;
% Vp
        end
            xx_ = xyIP_i(IP,1); 
        if (xx_<LL1 | xx_>LL2)
            b_ = b1_;
        else
            b_ = b2_;
        end
        slip_rate(cc,:) = [xx_ Vp]; 
        Vp = abs(Vp) ; 
        slip_(cc,:) = [xx_ opening_t opening_n]; 

        if itime > 1
            S1 = 1+theta_p0(cc,2)/dt;
            S2 = (1/dt + VP_0(cc)/L_);
            theta_ = S1/S2;
        end
        theta_p(cc,:) = [xx_ theta_] ; 
        VP(cc,1) = Vp ; 
        if itime > 1
            theta_n0 = theta_p0(cc,2);
        else
            theta_n0 = theta_0 ; 
        end
        
%         mu = a_*asinh((Vp*exp((mu0 + b_*log((V0*theta_)/L_))/a_))/(2*V0));
%         dmu = (a_*exp((mu0 + b_*log((V0*theta_)/L_))/a_))/(2*V0*((Vp^2*exp((2*mu0 + 2*b_*log((V0*theta_)/L_))/a_))/(4*V0^2) + 1)^(1/2)); 

        mu = a_*asinh((Vp*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*V0));
        dmu = (a_*(exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_)/(2*V0) - (Vp*b_*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*L_*V0*a_*(Vp/L_ + 1/dt))))/((Vp^2*exp((2*mu0 + 2*b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(4*V0^2) + 1)^(1/2); 
        
        mu_(cc,:) = [xx_ mu]; 
        
        parameter_t  = mu*(1*KK*opening_n+5e6);
        dparameter_t  = (1*KK*opening_n+1*5e6)*a1*dmu;
        dparameter_t2  = (1*KK*mu);

%         parameter_t = mu*(5e6+0*parameter_n);
%         dparameter_t  = 0 ; 
% opening_t
        if isnan(parameter_t)
            error ('NaN Value ---> Cohesive Law' ) 
        end
        cc = cc + 1 ; 

%         parameter_t = -parameter_t;
%         dparameter_t = -dparameter_t;
%         parameter_n = -parameter_n;
%         dparameter_n = -dparameter_n;
        MaxP = max(abs(parameter_t),MaxP);

%         tmat_tangent(:,:,IP) = 1*[dparameter_t dparameter_t2 ; 0 dparameter_n ];
        tmat_tangent(:,:,IP) = 1*[dparameter_t 0 ; 0 dparameter_n ];
        traction(:,1,IP) =  [parameter_t ; parameter_n] ;
end
%         traction = multiprod(tmat_secant,jump);
        [kmat] = kmatrixint(nmat,tmat_tangent,rotmat,detJ,b_b);
        sctrr = conint_cohesive(iint,1:end); 
        dof = [sctrr(1)*2-1 sctrr(1)*2 sctrr(2)*2-1 sctrr(2)*2   sctrr(3)*2-1 sctrr(3)*2   sctrr(4)*2-1 sctrr(4)*2  ];

        kglob(dof,dof ) = kglob(dof,dof) + 1*kmat;

        % internal force field (out-of-balance)
        [fint] = fmatrixint(nmat,traction,rotmat,detJ,b_b);
%         [fintglob] = assemblefintglob(fintglob,fint,ndf,conint_cohesive,iint,lin_quad,1);
        fintglob(dof,1) = fintglob(dof,1) - fint(:,1);
end

%}

% new version
%{
xi_top = [-1  1; 1  1] ;
xi_bot = [-1 -1; 1 -1] ;
rotmat = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha) ] ;
LINE_ = 0 ; 
% loop over coheisve elements
for iint = 1 : nim
    sc_iint = conint_cohesive(iint,:);
    uloc = zeros(8,1);    uloc(1:2:end) = u_elastic(2*sc_iint-1);    uloc(2:2:end) = u_elastic(2*sc_iint-0);
    vloc = zeros(8,1);    vloc(1:2:end) = vu(2*sc_iint-1);    vloc(2:2:end) = vu(2*sc_iint-0);
    coordloc = coord(sc_iint,:);
    kgg = 0 ; 
% loop over integration points (here, 2)
    for IP = 1 : size(rsIP_i,1)
        nmatrix_intra ;
        jump = rotmat*nmat*[uloc;uloc];
        jump_v = rotmat*nmat*[vloc;vloc];

        opening_t =  -(jump(1,:));
        opening_v_t =  -(jump_v(1,:));
        opening_n = -(jump(2,:));

        tmat_secant = zeros(2,2);
        tmat_tangent = zeros(2,2);
        traction = zeros(2,1) ;
        
        parameter_t = 0 ; 
        parameter_n = 0 ; 
        
%  Linear Cohesive Law
        KK = +1*1e8*1 ; 
        parameter_n = KK*opening_n;
        dparameter_n = KK;
        
        Vp = opening_v_t;

        xx_ = xyIP_i(IP,1); 
        if (xx_<LL1 | xx_>LL2)
            b_ = b1_;
        else
            b_ = b2_;
        end
        slip_rate(cc,:) = [xx_ Vp]; 
        Vp = abs(Vp) ; 
        slip_(cc,:) = [xx_ opening_t opening_n]; 

        if itime > 1
            S1 = 1+theta_p0(cc,2)/dt;
            S2 = (1/dt + VP_0(cc)/L_);
            theta_ = S1/S2;
        end
        theta_p(cc,:) = [xx_ theta_] ; 
        VP(cc,1) = Vp ; 
        if itime > 1
            theta_n0 = theta_p0(cc,2);
        else
            theta_n0 = theta_0 ; 
        end
        
        mu = a_*asinh((Vp*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*V0));
        dmu = (a_*(exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_)/(2*V0) - (Vp*b_*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*L_*V0*a_*(Vp/L_ + 1/dt))))/((Vp^2*exp((2*mu0 + 2*b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(4*V0^2) + 1)^(1/2); 
        
        mu_(cc,:) = [xx_ mu]; 
        
        parameter_t  = mu*(1*KK*opening_n+5e6);
        dparameter_t  = (1*KK*opening_n+1*5e6)*a1*dmu;
        dparameter_t2  = (1*KK*mu);

% opening_t
        if isnan(parameter_t)
            error ('NaN Value ---> Cohesive Law' ) 
        end
        cc = cc + 1 ; 

        MaxP = max(abs(parameter_t),MaxP);

        tmat_tangent = 1*[dparameter_t 0 ; 0 dparameter_n ];
        traction =  [parameter_t ; parameter_n] ;
        
        
 %% add stiffness
    dof = zeros(1,8) ;     
    dof(1:2:8) = (2*sc_iint-1); 
    dof(2:2:8) = (2*sc_iint-0); 
    dof = [dof dof] ; 
    
    kmat = 1 * nmat'*(rotmat'*tmat_tangent*rotmat)*nmat*detJ(IP,1)*detJ(IP,2)*b_b;
    fint = 1 * nmat'*(rotmat'*traction)*detJ(IP,1)*detJ(IP,2)*b_b;
    
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat(1:8,1:8);
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat(9:16,9:16);
%     kgg = kgg + kmat(1:8,1:8) ; 
%     kgg = kgg + kmat(9:16,9:16) ; 
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat;
    for ii_ = 1 : 16
        for jj_ = 1:16
            kglob(dof(ii_),dof(jj_)) = kglob(dof(ii_),dof(jj_)) + 1*kmat(ii_,jj_);
        end
    end
    
    kgg = kgg + kmat ; 

    fintglob(dof(1:8),1) = fintglob(dof(1:8),1) - fint(1:8,1);
    fintglob(dof(9:16),1) = fintglob(dof(9:16),1) - fint(9:16,1);

    LINE_ = LINE_ + detJ(IP,1)*detJ(IP,2)*b_b ;

    end
%     dof
%     sparse(kglob)
%     pause
end
% LINE_
% pause
%}

% new version 2
% %{
xi_top = [-1  1; 1  1] ;
xi_bot = [-1 -1; 1 -1] ;
rotmat = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha) ] ;
LINE_ = 0 ; 
% loop over coheisve elements
for iint = 1 : nim
    sc_iint = conint_cohesive(iint,:);
    eps_ =  150; 
    vloc = zeros(8,1);    vloc(1:2:end) = vu(2*sc_iint-1);    vloc(2:2:end) = vu(2*sc_iint-0);
    coordloc = coord(sc_iint,:);
    x_m = mean(coordloc(:,1));
%     pause
    eps_ =  0.1*D/51/2-0.0001 ;
    eps_ =  100; 
    x_top = [ x_m D/2 + eps_ ]; 
    x_bot = [ x_m D/2 - eps_ ]; 
    
    if itime == 1
        [iel_t,sc_iint_t,local_t] = find_which_elem ( x_top , element , node );
        [iel_b,sc_iint_b,local_b] = find_which_elem ( x_bot , element , node );
        DATA(iint,:) = [iel_t sc_iint_t local_t iel_b sc_iint_b local_b];
    else
        iel_t = DATA(iint,1); 
        sc_iint_t = DATA(iint,2:5); 
        local_t = DATA(iint,6:7); 
        iel_b = DATA(iint,8); 
        sc_iint_b = DATA(iint,9:12); 
        local_b = DATA(iint,13:14); 
    end
    xi_top(:,2) = local_t(2);
    xi_bot(:,2) = local_b(2);
    
    
    uloc_t = zeros(8,1);    uloc_t(1:2:end) = u_elastic(2*sc_iint_t-1);    uloc_t(2:2:end) = u_elastic(2*sc_iint_t-0);
    vloc_t = zeros(8,1);    vloc_t(1:2:end) = vu(2*sc_iint_t-1);    vloc_t(2:2:end) = vu(2*sc_iint_t-0);
    
    uloc_b = zeros(8,1);    uloc_b(1:2:end) = u_elastic(2*sc_iint_b-1);    uloc_b(2:2:end) = u_elastic(2*sc_iint_b-0);
    vloc_b = zeros(8,1);    vloc_b(1:2:end) = vu(2*sc_iint_b-1);    vloc_b(2:2:end) = vu(2*sc_iint_b-0);

    kgg = 0 ; 
% loop over integration points (here, 2)
    for IP = 1 : size(rsIP_i,1)
        nmatrix_intra ;
        jump = rotmat*nmat*[uloc_b;uloc_t];
        jump_v = rotmat*nmat*[vloc_b;vloc_t];

        opening_t =  -(jump(1,:));
        opening_v_t =  -(jump_v(1,:));
        opening_n = -(jump(2,:));

        tmat_secant = zeros(2,2);
        tmat_tangent = zeros(2,2);
        traction = zeros(2,1) ;
        
        parameter_t = 0 ; 
        parameter_n = 0 ; 
        
%  Linear Cohesive Law
        KK = +1*1e8*1 ; 
        parameter_n = KK*opening_n;
        dparameter_n = KK;
        
        Vp = opening_v_t;

        xx_ = xyIP_i(IP,1); 
        if (xx_<LL1 | xx_>LL2)
            b_ = b1_;
        else
            b_ = b2_;
        end
        slip_rate(cc,:) = [xx_ Vp]; 
        Vp = abs(Vp) ; 
        slip_(cc,:) = [xx_ opening_t opening_n]; 

        if itime > 1
            S1 = 1+theta_p0(cc,2)/dt;
            S2 = (1/dt + VP_0(cc)/L_);
            theta_ = S1/S2;
        end
        theta_p(cc,:) = [xx_ theta_] ; 
        VP(cc,1) = Vp ; 
        if itime > 1
            theta_n0 = theta_p0(cc,2);
        else
            theta_n0 = theta_0 ; 
        end
        
        mu = a_*asinh((Vp*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*V0));
        dmu = (a_*(exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_)/(2*V0) - (Vp*b_*exp((mu0 + b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(2*L_*V0*a_*(Vp/L_ + 1/dt))))/((Vp^2*exp((2*mu0 + 2*b_*log((V0*(theta_n0/dt + 1))/(L_*(Vp/L_ + 1/dt))))/a_))/(4*V0^2) + 1)^(1/2); 
        
        mu_(cc,:) = [xx_ mu]; 
        
        parameter_t  = mu*(1*KK*opening_n+5e6);
        dparameter_t  = (1*KK*opening_n+1*5e6)*a1*dmu;
        dparameter_t2  = (1*KK*mu);

% opening_t
        if isnan(parameter_t)
            error ('NaN Value ---> Cohesive Law' ) 
        end
        cc = cc + 1 ; 

        MaxP = max(abs(parameter_t),MaxP);

        tmat_tangent = 1*[dparameter_t 0 ; 0 dparameter_n ];
        traction =  [parameter_t ; parameter_n] ;
        
        
 %% add stiffness
    dof_b = zeros(1,8) ;     
    dof_b(1:2:8) = (2*sc_iint_b-1); 
    dof_b(2:2:8) = (2*sc_iint_b-0); 
    dof_t = zeros(1,8) ;     
    dof_t(1:2:8) = (2*sc_iint_t-1); 
    dof_t(2:2:8) = (2*sc_iint_t-0); 
    dof = [dof_b dof_t] ; 
    
    kmat = 1 * nmat'*(rotmat'*tmat_tangent*rotmat)*nmat*detJ(IP,1)*detJ(IP,2)*b_b;
    fint = 1 * nmat'*(rotmat'*traction)*detJ(IP,1)*detJ(IP,2)*b_b;
    
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat(1:8,1:8);
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat(9:16,9:16);
%     kgg = kgg + kmat(1:8,1:8) ; 
%     kgg = kgg + kmat(9:16,9:16) ; 
%     kglob(dof,dof ) = kglob(dof,dof) + 1*kmat;
    for ii_ = 1 : 16
        for jj_ = 1:16
            kglob(dof(ii_),dof(jj_)) = kglob(dof(ii_),dof(jj_)) + 1*kmat(ii_,jj_);
        end
    end
    
    kgg = kgg + kmat ; 

    fintglob(dof(1:8),1) = fintglob(dof(1:8),1) - fint(1:8,1);
    fintglob(dof(9:16),1) = fintglob(dof(9:16),1) - fint(9:16,1);

    LINE_ = LINE_ + detJ(IP,1)*detJ(IP,2)*b_b ;

    end
%     dof
%     sparse(kglob)
%     pause
end
% LINE_
% pause
%}

% new version with volume integral
%{
imat = 0 ; 
rotmat = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha) ] ;
LINE_ = 0 ;
kgg = 0 ; QQ = [] ; QQ_ = [] ; 
for iel = 1 : nem 
    if rem(iel,1000)==0
%         iel/nem
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
    [W,Q] = quadrature(2,'GAUSS',2);    
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        
        % Shape functions and its derivatives
        [N,dNdxi] = Lag_basis('Q4',pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
        
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        xgp = N'*node(sctr,:) ; 

%         get phase field parameter at this gp
        phi = N'*PHI(sctr) ;
        dphi = dNdx'*PHI(sctr);
                
        if abs(phi)>1e-11
            gamma_ = (1/2/l_)*phi*phi + l_/2*(dphi'*dphi) ;

            iint = find_which_interface ( xgp , node , conint_cohesive );

            sc_iint = conint_cohesive(iint,:);
            uloc = zeros(8,1);    uloc(1:2:end) = u_elastic(2*sc_iint-1);    uloc(2:2:end) = u_elastic(2*sc_iint-0);
            vloc = zeros(8,1);    vloc(1:2:end) = vu(2*sc_iint-1);    vloc(2:2:end) = vu(2*sc_iint-0);
            coordloc = coord(sc_iint,:);

            add_interface_terms_phasefield ;
            QQ = [QQ ; xgp ] ; 
%             QQ_ = [ QQ_ ; seg_t ] ;
%             QQ_ = [ QQ_ ; seg_b ] ; 
        end
        
    end
end
imat/L
%}
% LINE_

% figure
% plot(QQ(:,1),QQ(:,2),'rsq')
% hold on
% plot(QQ_(:,1),QQ_(:,2),'bsq')
% 
% pause

rhs = fintglob; % no external forces, only displacements


if itr == 1 
   fintglob = f_int_ ;  
end

%% Solve system   

%   total prescribed disp vector 
    up = 0*anode'; 
    up (dispNodes1*2-1) = ap2 ; 
    up (dispNodes2*2-1) = ap1 ; 
%  prescribed values 
    ap = zeros(length(presDofs),1) ; 
    ap(1:length(dispNodes1))=ap2 ; 
    ap(length(dispNodes1)+1:length(dispNodes1)+length(dispNodes2)) = ap1 ;
    A = kglob + a0*mglob ;
    
    Converged_last_step = mglob*a0*u_elastic_0 + mglob*a2*vu0 + mglob*a4*au0 ;

    r = [rhs-dynamic*Converged_last_step+1*Fext] ;
    [A,r] = Essential_FEM (A,r,dispNodes,presDofs,anode,itr,ap);
    dd = A\r;
    
%     if itime == 1 && itr == 1
%         [LLLL,UUUU] = ilu(A,struct('type','ilutp','droptol',1e-11));
%     end
%     tol = 1e-13;  
%     maxit = 1000; 

%     [X,FLAG,RELRES,ITER,RESVEC] = gmres(A,r,[],tol,maxit,LLLL,UUUU);
%     norm(X-dd)
%     dd = X ; 
%     figure('visible','on')
%     semilogy([1:length(RESVEC)],RESVEC/norm(r),'rsq-')
%     pause
%     if rem(itime,10)==0
%         length(RESVEC)
%     end
%     du = dd(1:DOF) ;
    du = 0*up ;
    du(setdiff(anode,presDofs))=dd ; 

% update solution
    u_elastic = u_elastic - du ;
 
% replace actual values for essential boundaries
    u_elastic (dispNodes1*2-1) = u_elastic_0 (dispNodes1*2-1) + ap2 ; 
    u_elastic (dispNodes2*2-1) = u_elastic_0 (dispNodes2*2-1) + ap1 ; 

    normm = norm(du);
    
%     f_int_ = kglob * u_elastic ; %+ a0* mglob * u_elastic ;         
%     Reaction = sum(f_int_(dispNodes2*2-1)) ; 
    disp ([ 'normm (' num2str(itr) ') = ' num2str(normm) ' at  itime = ' num2str(itime) ' Reaction: ' num2str(0)  ]); 

    %     disp ([ 'syy_tip = ' num2str(stress_tip(2)) ]); 
    DISP = u_elastic;
%}

%         MaxP

%% calculate acceleration, velocity, ...
        au = a0*(u_elastic - u_elastic_0) - (a2*vu0) - (a4*au0) ;
        vu = a1*(u_elastic - u_elastic_0) - (a3*vu0) - (a5*au0) ;
        
        
end
% figure
% plot(slip_(:,1),slip_(:,2),'sq')
% pause

        u_elastic_0 = u_elastic ; 
        Time0 = Time ; 
        au0 = au ;
        vu0 = vu ;
        theta_p0 = theta_p ;
        maxVp = max(slip_rate(:,2)) ; 
        VP_0 = VP ; 
%% post process
        f_int_ = kglob * u_elastic ; %+ a0* mglob * u_elastic ; 
        
%         F_INTERNAL;
        
        Reaction = sum(f_int_(dispNodes2*2-1)) ; 
%         mul=[mul; Reaction  u_elastic(2*uln-1)-u_elastic(2*lln-1)]; 
        mul=[mul; Reaction  u_elastic(2*uln-1)-u_elastic(2*lln-1) Time/60/60/24/365 dt maxVp]; 
%         figure(h1)
%         clf
%         plot (mul(:,2),mul(:,1),'r.-')
%         hold on
%         load LD_ref.mat 'LD'
%         plot (LD(:,2),LD(:,1),'bsq')
%         pause    (0.01)  
if rem(itime , print_step1 ) == 0    

    figure('visible','off')
    plot(slip_(:,1),slip_(:,2),'rsq')
    ylabel('slip')
    eval(['print -djpeg out/slip/jj' num2str(itime) , '.jpeg'])
    save_mat ( 'slip' , [slip_(:,1),slip_(:,2)] , itime )
    
    figure('visible','off')
    plot(slip_(:,1),slip_(:,3),'rsq')
    ylabel('slip')
    eval(['print -djpeg out/slip_n/jj' num2str(itime) , '.jpeg'])

    figure('visible','off')
    plot(slip_rate(:,1),slip_rate(:,2),'rsq')
    ylabel('slip_rate')
    eval(['print -djpeg out/slip_rate/jj' num2str(itime) , '.jpeg'])
    save_mat ( 'slip_rate' , [slip_rate(:,1),slip_rate(:,2)] , itime )

    figure('visible','off')
    plot(mu_(:,1),mu_(:,2),'bsq')
    ylabel('mu')
    eval(['print -djpeg out/mu/jj' num2str(itime) , '.jpeg'])

    figure('visible','off')
    plot(theta_p(:,1),theta_p(:,2),'rsq')
    ylabel('mu')
    eval(['print -djpeg out/theta/jj' num2str(itime) , '.jpeg'])

    figure('visible','off')
    plot(mul(:,3),mul(:,1)/1e6/L,'r-sq','MarkerSize',1)
    ylabel('LD')
    hold on
%     load /home/mohsen/Desktop/PhD/MyMatlabToolkit/Herr1.mat
%     plot(Herr1(:,1),Herr1(:,2),'b-')
    ylabel('Stress[MPa]')
    eval(['print -djpeg out/LD.jpeg'])

    fac = 10000 ; 
    DISP = u_elastic;
    figure('visible','off')
    hold on
    DISP_ =  [DISP(1:2:end) DISP(2:2:end)];
    node_deformed = node + fac * [DISP(1:2:end) DISP(2:2:end)] ; 
    plot(node_deformed(:,1),node_deformed(:,2),'r.')
    axis equal
    eval(['print -djpeg out/deformed.jpeg'])
    
    h1 = figure('visible','off');
    semilogy(mul(:,3),mul(:,end),'bsq-','MarkerSize',1)
    eval(['print -djpeg out/maxVp' , '.jpeg'])
    close(h1)    
    
close all 

end

    dlmwrite(['out/LD' '.out'],...
            [mul(:,3) mul(:,end) mul(:,4)] ,'delimiter' , '\t','precision',11  )

%     pause

if rem(itime , print_step2 ) == 0    
%%
    q = [] ;
    Disp = [] ; STRESS = [ ];  vol_ = 0 ; 
    for iel = 1 : size(element,1) 
    %     iel/size(element,1) 
        sctr = element(iel,:); % element connectivity
        nn   = length(sctr);   % number of nodes per element
        ke   = 0 ;             % elementary stiffness matrix
        order = 2 ;
        [W,Q] = quadrature(order,'GAUSS',2);

    % -----------------------------------------------
        for igp = 1 : size(W,1)
            gpnt = Q(igp,:);
            [N,dNdxi]=lagrange_basis('Q4',gpnt);
            Gpnt = N' * node(sctr,:); % global GP
            q = [q;Gpnt];
        end
        % -----------------------------------------------

            N1  = element(iel,1);                                                  % Node 1 for current element
            N2  = element(iel,2);                                                  % Node 2 for current element
            N3  = element(iel,3);                                                  % Node 3 for current element
            N4  = element(iel,4);                                                  % Node 4 for current element
    %       
        sctrBstd  = [N1*2-1 N1*2 N2*2-1 N2*2 N3*2-1 N3*2 N4*2-1 N4*2];             % Traditional index locations
        sctrBstd = 0*sctrBstd ;
        sctrBstd(1:2:end) = 2*element(iel,:)-1;
        sctrBstd(2:2:end) = 2*element(iel,:);
        U     = u_elastic(sctrBstd);
%         U     = 1e9*du(sctrBstd);


        for kk = 1 : size(W,1)
    %         pt = Q(kk,:);                             % quadrature point
            % B matrix
            gpnt = Q(kk,:);
            [N,dNdxi]=lagrange_basis('Q4',gpnt);
            J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
            invJ0 = inv(J0);
            dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
            Nx = dNdx(:,1) ; 
            Ny = dNdx(:,2) ; 

            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0      Nx(4)   0 ;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)  0   Ny(4) ;...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];

            dispx  = N' * U(1:2:end) ; 
            dispy  = N' * U(2:2:end) ; 
            Disp = [ Disp ; dispx dispy ]; 
            strain = Bu*U;
            stress = Cm*strain;
            STRESS =  [ STRESS ; stress' ] ; 
            vol_ = vol_  + W(kk,1) * det(J0) ;
        end                  % end of looping on GPs
    end                      % end of looping on elements

        QQQ = q; 
        tri = delaunay(QQQ(:,1),QQQ(:,2));
    %     vtfile =['out/stress' num2str(itime)] ;
        vtfile =['out/stress' ] ;
        VTKPostProcess(QQQ+0*Disp,tri,1,'Tri3',vtfile,STRESS,Disp)
    ! paraview out/stress.vtu &

% 
pause(1)
end
    
    
%%
fac = 1000 ; 
DISP = u_elastic;
% figure(h2)
% clf
% hold on
% DISP_ =  [DISP(1:2:end) DISP(2:2:end)];
% node_deformed = node + fac * [DISP(1:2:end) DISP(2:2:end)] ; 
% plot(node_deformed(:,1),node_deformed(:,2),'rsq')
% axis equal

if rem(itime,100)==0
    save all.mat
end

% disp "converged!" 
end

if ifprofile
    profile viewer
    profsave
end
return
%%
figure
load 'out_l6000_51/slip/jj500.mat'
% plot(ff(:,1),ff(:,2),'ro')
hold on
load 'out_l6000_251/slip/jj500.mat'
plot(ff(:,1),ff(:,2),'b.')
load 'out_l6000_551/slip/jj500.mat'
plot(ff(:,1),ff(:,2),'k.')
load 'out_l6000_751/slip/jj500.mat'
plot(ff(:,1),ff(:,2),'g.')
% load 'out_l3000_751/slip/jj500.mat'
% plot(ff(:,1),ff(:,2),'c.')
load 'out_l12000_1151/slip/jj500.mat'
plot(ff(:,1),ff(:,2),'c.')
load '../out_ref/slip/jj500.mat'
plot(ff(:,1),ff(:,2),'r.')


%%

h1 = figure('visible','on');
ff = dlmread('out_l6000_251/LD.out');
semilogy(ff(:,1),ff(:,2),'bsq-','MarkerSize',1)
hold on
ff = dlmread('out_l6000_551/LD.out');
semilogy(ff(:,1),ff(:,2),'ksq-','MarkerSize',1)
ff = dlmread('out_l6000_751/LD.out');
semilogy(ff(:,1),ff(:,2),'gsq-','MarkerSize',1)
ff = dlmread('../out_ref/LD.out');
semilogy(ff(:,1),ff(:,2),'rsq-','MarkerSize',1)

%%
figure
plot(node(:,1),node(:,2),'k.')
hold on
plot([0 L],[D/2 D/2],'b-')
plot([0 L],[D/2+eps_ D/2+eps_],'k-','LineWidth',2)
plot([0 L],[D/2-eps_ D/2-eps_],'k-','LineWidth',2)
plot([0 L],[D/2+l_ D/2+l_],'c-','LineWidth',2)
plot([0 L],[D/2-l_ D/2-l_],'c-','LineWidth',2)
%%
figure
% load '/home/mohsen/Desktop/Interface_phasefield/out_80/slip/jj450.mat'
% plot(ff(:,1),ff(:,2),'b-','MarkerSize',1)
hold on
load '/home/mohsen/Desktop/Interface_phasefield/out_180/slip/jj450.mat'
plot(ff(:,1),ff(:,2),'r-','MarkerSize',1)
% load '/home/mohsen/Desktop/Interface_phasefield/out_280/slip/jj450.mat'
% plot(ff(:,1),ff(:,2),'k-','MarkerSize',1)
load '/home/mohsen/Desktop/Interface_phasefield/out_180_151/slip/jj450.mat'
plot(ff(:,1),ff(:,2),'c-','MarkerSize',1)
%%

figure
% load '/home/mohsen/Desktop/Interface_phasefield/out_80/slip/jj450.mat'
% plot(ff(:,1),ff(:,2),'g-','MarkerSize',1)
hold on
% load '/home/mohsen/Desktop/Interface_phasefield/out_180/slip/jj410.mat'
% plot(ff(:,1),ff(:,2),'r-','MarkerSize',1)
% load '/home/mohsen/Desktop/Interface_phasefield/out_180_eps2/slip/jj410.mat'
% plot(ff(:,1),ff(:,2),'b-','MarkerSize',1)
% load '/home/mohsen/Desktop/Interface_phasefield/out_180_151/slip/jj410.mat'
% plot(ff(:,1),ff(:,2),'c-','MarkerSize',1)
load /home/mohsen/Desktop/out_ref/slip/jj450.mat
plot(ff(:,1),ff(:,2),'k-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_5/out/slip/jj450.mat
plot(ff(:,1),ff(:,2),'r-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_6/out/slip/jj450.mat
plot(ff(:,1),ff(:,2),'g-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_7/out/slip/jj450.mat
plot(ff(:,1),ff(:,2),'y-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_8/out/slip/jj450.mat
plot(ff(:,1),ff(:,2),'m-','LineWidth',1)
% load '/home/mohsen/Desktop/Interface_phasefield/out/slip/jj410.mat'
% plot(ff(:,1),ff(:,2),'-','LineWidth',2)
%%

figure
hold on
load /home/mohsen/Desktop/Interface_phasefield_2/out/slip/jj810.mat
plot(ff(:,1),ff(:,2),'r-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_2/out/slip/jj810.mat
plot(ff(:,1),ff(:,2),'g-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_3/out/slip/jj810.mat
plot(ff(:,1),ff(:,2),'y-','LineWidth',1)
load /home/mohsen/Desktop/Interface_phasefield_2_4/out/slip/jj810.mat
plot(ff(:,1),ff(:,2),'m-','LineWidth',1)

