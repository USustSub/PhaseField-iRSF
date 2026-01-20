        
%         r = rsIP_i(IP,1);
        x1 = coordloc(1,1);
        y1 = coordloc(1,2);
        x2 = coordloc(2,1);
        y2 = coordloc(2,2);
%                 
%         N1 = 0.5*(1-r);
%         N2 = 0.5*(1+r);        
        
        
        nmat = zeros(2,16);
% bot
        local_ = xi_bot(IP,:);
        [N,~]=lagrange_basis('Q4',local_);
        nmat(1,1:2:8) = N ;
        nmat(2,2:2:8) = N ;
% top
        local_ = xi_top(IP,:);
        [N,~]=lagrange_basis('Q4',local_);
        nmat(1,8+[1:2:8]) = -N ;
        nmat(2,8+[2:2:8]) = -N ;

% center
        local_(:,2) = 0 ;
        [N,~]=lagrange_basis('Q4',local_);

        seg_ = N'*coordloc; 
        xyIP_i(IP,:) = seg_;
        
        length_i = sqrt((y2-y1)^2 + (x2-x1)^2);
        
        detJ(IP,1) = (1/2)*length_i;
        detJ(IP,2) = rsIP_i(IP,2);