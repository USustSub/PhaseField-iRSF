
        xi_top = [pt(1)  1 ] ;
        xi_bot = [pt(1) -1 ] ;

        x1 = coordloc(1,1);
        y1 = coordloc(1,2);
        x2 = coordloc(2,1);
        y2 = coordloc(2,2);

        nmat = zeros(2,16);
% bot
        local_ = xi_bot(:,:);
        [N,~]=lagrange_basis('Q4',local_);
        nmat(1,1:2:8) = N ;
        nmat(2,2:2:8) = N ;
        seg_b = N'*coordloc; 

% top
        local_ = xi_top(:,:);
        [N,~]=lagrange_basis('Q4',local_);
        nmat(1,8+[1:2:8]) = -N ;
        nmat(2,8+[2:2:8]) = -N ;
        seg_t = N'*coordloc; 

% center
        local_(:,2) = 0 ;
        [N,~]=lagrange_basis('Q4',local_);

        seg_ = N'*coordloc; 
%         xyIP_i(IP,:) = seg_;
        
        length_i = sqrt((y2-y1)^2 + (x2-x1)^2);
        
%         detJ(IP,1) = (1/2)*length_i;
%         detJ(IP,2) = rsIP_i(IP,2);


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

        xx_ = xgp(1); 
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
    
    kmat = gamma_ * nmat'*(rotmat'*tmat_tangent*rotmat)*nmat*W(kk)*det(J0);
    fint = gamma_ * nmat'*(rotmat'*traction)*W(kk)*det(J0);
    
    imat = imat + gamma_ * W(kk)*det(J0) ;
    
    for ii_ = 1 : 16
        for jj_ = 1:16
            kglob(dof(ii_),dof(jj_)) = kglob(dof(ii_),dof(jj_)) + 1*kmat(ii_,jj_);
        end
    end
    
    kgg = kgg + kmat ; 

    fintglob(dof(1:8),1) = fintglob(dof(1:8),1) - fint(1:8,1);
    fintglob(dof(9:16),1) = fintglob(dof(9:16),1) - fint(9:16,1);
    
    LINE_ = LINE_ + W(kk)*det(J0);