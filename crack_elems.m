% elem_candidates = 
dmax  = 1.1 ;              % radius = dmax * nodal spacing
numnode = size(node,1);

deltaX = L/(nnx-1);
di     = ones(1,numnode)*dmax*deltaX ;


seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
xTip  = xCr(2,:);

% Only normal level sets are computed and stored in nodes
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
    dt(i)   = norm([x y]-xTip); % distances from node I to tip
end

set1 = find(abs(ls(:,1)) - di' < 0);
set2 = find(ls(:,2) < 0);
set3 = find(dt - di < 0);        % tip nodes
set4 = intersect(set1,set2);     % split nodes

% % some nodes belong to both sets, remove tip
% split_nodes = setdiff(set4,set3);
% split_nodes   = [split_nodes ;set3']; 
% 
% split_nodes((node(split_nodes,1)<0.2))=[];
% 
% [conn_nodes ] = find_conn_nodes (node , element ) ;
% 
% desir_elems =conn_nodes(split_nodes,2:end);
% desir_elems = unique(desir_elems(:));
% desir_elems(desir_elems==0) = [] ; 





enrich_node = zeros(numnode,1);

count1 = 0;
count2 = 0;
for iel = 1 : size(element,1)
    sctr = element(iel,:);
    phi  = ls(sctr,1);
    psi  = ls(sctr,2);
    if all(node(sctr,1)>0.2)
    if ( max(phi)*min(phi) < 0 )
        if max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1) = iel;
            enrich_node(sctr)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem(count2) = iel;
            enrich_node(sctr)   = 2;
        end
    end
    end
end
split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);
split_nodes = [split_nodes ; tip_nodes ] ; 

% desir_elems =[split_elem tip_elem];
% desir_elems = unique(desir_elems(:));
% desir_elems(desir_elems==0) = [] ; 

desir_elems =conn_nodes(split_nodes,2:end);
desir_elems = unique(desir_elems(:));
desir_elems(desir_elems==0) = [] ; 


%% find crack nodes
%     crack_nodes = [0:0.01:0.5];
%     crack_nodes = [crack_nodes' crack_nodes'*0+0.5];
    crack_nodes = []  ;
    crack_nodes2 = []  ;
    crack_nodes3 = []  ;
%     eps_ = 0.00012 ; 
    eps_ = deltaX/2*2. ; 
    for i_ = 1 : n_Cr
        crack_nodes(i_,:) = xCr(1,:)+seg*i_/n_Cr;
%         if P_order > 1
            crack_nodes2(i_,:) = xCr(1,:)+null(seg)'*eps_+seg*i_/n_Cr;
            crack_nodes3(i_,:) = xCr(1,:)-null(seg)'*eps_+seg*i_/n_Cr;
%         end
    end
    
%     figure
    plot(crack_nodes(:,1),crack_nodes(:,2),'k.')
    hold on
    plot(crack_nodes2(:,1),crack_nodes2(:,2),'r.')
    plot(crack_nodes3(:,1),crack_nodes3(:,2),'b.')
    if P_order > 1
        crack_nodes = [crack_nodes ; crack_nodes2 ; crack_nodes3 ] ;
    end
    
    crack_nodes2 = [crack_nodes2 ; crack_nodes3 ] ;

%% find local coordinate of each node
crack_elems_ = crack_nodes(:,1)*0 ; 
crack_local = crack_nodes*0 ; 
for ip = 1 : length(crack_nodes)
    xx = crack_nodes(ip,:);
    
    for ii = 1 : length(desir_elems)
            sc = node(element(desir_elems(ii),:),:);
            min_x = min(sc(:,1));
            max_x = max(sc(:,1));
            min_y = min(sc(:,2));
            max_y = max(sc(:,2));
            C.type = 'Q4';
            if xx(1)>=min_x && xx(1)<=max_x && xx(2)>=min_y && xx(2)<=max_y
                [local] = element_coordinate (xx , element(desir_elems(ii),:) , node , C);
                crack_elems_(ip) = desir_elems(ii);
                crack_local(ip,:) = local ; 
            end
    end
    
end



crack_elems_2 = crack_nodes2(:,1)*0 ; 
crack_local2 = crack_nodes2*0 ; 
for ip = 1 : length(crack_nodes2)
    xx = crack_nodes2(ip,:);
    
    for ii = 1 : length(desir_elems)
            sc = node(element(desir_elems(ii),:),:);
            min_x = min(sc(:,1));
            max_x = max(sc(:,1));
            min_y = min(sc(:,2));
            max_y = max(sc(:,2));
            C.type = 'Q4';
            if xx(1)>=min_x && xx(1)<=max_x && xx(2)>=min_y && xx(2)<=max_y
                [local] = element_coordinate (xx , element(desir_elems(ii),:) , node , C);
                crack_elems_2(ip) = desir_elems(ii);
                crack_local2(ip,:) = local ; 
            end
    end
    
end
