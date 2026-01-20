%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Phase-field continuum model for faults with RSF  
% By: Mohsen Goudarzi (2023, Utrecht), Goudarzi.mohsen@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear all ; 
    close all ;
    clc 

    addpath 'PFEM_functions'
    
    ! rm -r out
    ! mkdir out
    ! mkdir out/paraview
    ! mkdir out/slip
    ! mkdir out/slip_n
    ! mkdir out/slip_rate
    ! mkdir out/theta
    ! mkdir out/mu
    ! mkdir out/deformed

    global volume_ LINE_
    warning('off', 'MATLAB:nargchk:deprecated')
    
    dynamic = 1 ; 
    msh_type = 'Q4' ; 
    
    Time = 0 ;
    print_step1 = 50 ; 
    print_step2 = 200000000 ;
    
    nn = 540  ; % mesh density
    dt = 1e6 ;
    h__ = 1 ; 
    disp(['Time step is ' num2str(dt/60/60/24/365) ' years'])
    
    if strcmp(msh_type,'Q4')
        Input_ ; 
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
    plot([LL1 LL1 ] , [0 L],'k-','LineWidth',4);
    plot([LL2 LL2] , [0 L],'k-','LineWidth',4);

    kkk = 1e-11 ; % small value

pause(1)
%%
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

    Time0 = 0 ; 


%% solve for phase field parameter
    Solve_for_phase_field

%% time stepping
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

    

    xi_top = [-1  1; 1  1] ;
    xi_bot = [-1 -1; 1 -1] ;
    rotmat = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha) ] ;
    LINE_ = 0 ; 
    
    
% loop over coheisve elements
    for iint = 1 : nim
        sc_iint = conint_cohesive(iint,:);
        uloc = zeros(8,1);    uloc(1:2:end) = u_elastic(2*sc_iint-1);    uloc(2:2:end) = u_elastic(2*sc_iint-0);
        vloc = zeros(8,1);    vloc(1:2:end) = vu(2*sc_iint-1);    vloc(2:2:end) = vu(2*sc_iint-0);        coordloc = coord(sc_iint,:);
        x_m = mean(coordloc(:,1));

        
    eps_ =  150; 
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

    end


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

    du = 0*up ;
    du(setdiff(anode,presDofs))=dd ; 

% update solution
    u_elastic = u_elastic - du ;
 
% replace actual values for essential boundaries
    u_elastic (dispNodes1*2-1) = u_elastic_0 (dispNodes1*2-1) + ap2 ; 
    u_elastic (dispNodes2*2-1) = u_elastic_0 (dispNodes2*2-1) + ap1 ; 

    normm = norm(du);
    
    disp ([ 'normm (' num2str(itr) ') = ' num2str(normm) ' at  itime = ' num2str(itime) ' Reaction: ' num2str(0)  ]); 

    DISP = u_elastic;

%% calculate acceleration, velocity, ...
        au = a0*(u_elastic - u_elastic_0) - (a2*vu0) - (a4*au0) ;
        vu = a1*(u_elastic - u_elastic_0) - (a3*vu0) - (a5*au0) ;
        
end

        u_elastic_0 = u_elastic ; 
        Time0 = Time ; 
        au0 = au ;
        vu0 = vu ;
        theta_p0 = theta_p ;
        maxVp = max(slip_rate(:,2)) ; 
        VP_0 = VP ; 
%% post process
        f_int_ = kglob * u_elastic ; %+ a0* mglob * u_elastic ; 
        
        Reaction = sum(f_int_(dispNodes2*2-1)) ; 
        mul=[mul; Reaction  u_elastic(2*uln-1)-u_elastic(2*lln-1) Time/60/60/24/365 dt maxVp]; 

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
%% post process 2

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
        vtfile =['out/stress' ] ;
        VTKPostProcess(QQQ+0*Disp,tri,1,'Tri3',vtfile,STRESS,Disp)
        ! paraview out/stress.vtu &
        pause(1)
end
    
    
%%
fac = 1000 ; 
DISP = u_elastic;

if rem(itime,100)==0
    save all.mat
end

end



