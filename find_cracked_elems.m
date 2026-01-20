function [split_elem , tip_elem , split_elem_all] = find_cracked_elems ( node , element,  xCr) 
%% first for real crack
xTip = xCr(2,:) ; 
seg   = xCr(2,:) - xCr(1,:);   % tip segment

% find elements that needs to be cracked:
% the tangent LS. The i-th row is for node I
x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t   = 1/norm(seg)*seg;
for i = 1 : size(node,1)
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip)*t';  % tangent LS
end

% enrich_node = zeros(size(node,1),1);
split_elem = [] ; 
tip_elem = [] ; 
count1 = 0;
count2 = 0;
for iel = 1 : size(element,1)
    sctr = element(iel,:);
    phi  = ls(sctr,1);
    psi  = ls(sctr,2);
    if ( max(phi)*min(phi) < 0 )
        if max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1) = iel;
%             enrich_node(sctr)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem(count2) = iel;
%             enrich_node(sctr)   = 2;
        end
    end
end
% enrich_node(enrich_node==2) = 1 ; 

%% now for all cracked elements
xCr_all = xCr ; 
xCr_all(2,1) = max(node(:,1))+1 ; 
xTip_all = xCr_all(2,:) ; 
seg_all   = xCr_all(2,:) - xCr_all(1,:);   % tip segment


% find elements that needs to be cracked:
% the tangent LS. The i-th row is for node I
x0  = xCr_all(1,1); y0 = xCr_all(1,2);
x1  = xCr_all(2,1); y1 = xCr_all(2,2);
t   = 1/norm(seg_all)*seg_all;
for i = 1 : size(node,1)
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip_all)*t';  % tangent LS
end

% enrich_node = zeros(size(node,1),1);
split_elem_all = [] ; 
tip_elem_all = [] ; 
count1 = 0;
count2 = 0;
for iel = 1 : size(element,1)
    sctr = element(iel,:);
    phi  = ls(sctr,1);
    psi  = ls(sctr,2);
    if ( max(phi)*min(phi) < 0 )
        if max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem_all(count1) = iel;
%             enrich_node(sctr)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem_all(count2) = iel;
%             enrich_node(sctr)   = 2;
        end
    end
end
% enrich_node(enrich_node==2) = 1 ; 
split_elem_all=  [split_elem_all tip_elem_all]; 
end