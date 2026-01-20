function  [node_p2 , element_p2]  =  P1toP2_func_v3 ( node , element , method ,  p_step , element_des)
% This function transforms a linear Quad element to higher order Quads 
% It works with P-FEM. 
% Mohsen Goudarzi September, 2016, Delft 

%% Find all the elements connected to each element
% conn_nodes =  find_connectivity_nodes ( node , element  );
[conn_nodes ] = find_conn_nodes (node , element ) ;
[conn_elems ] = find_conn_elements (node , element ,conn_nodes) ;


%% increase order of each element 
if strcmp(method,'Lagrange')
        edges = [ 1 2 ; 2 3 ; 3 4 ; 1 4 ] ;
elseif strcmp(method,'Lobatto')
        edges = [ 1 4 ; 2 3 ; 1 2 ; 3 4 ] ;
end

        node_p2 = node ; 
        nd = 4+4*(p_step-1)+(p_step-1).^2 ;
        element_p2 = [ element  zeros(size(element,1),nd-4)]   ; 

for pp = 2 : p_step 
    pp_ = pp - 1 ; 
    init = 4+4*(pp_-1)+(pp_-1).^2;
for nn = 1 : length( element_des) 
        kk = element_des(nn) ;
        sctr = element(kk,:) ;
        c_elems = conn_elems{kk} ;  c_elems(1) = [] ; 
    for ll = 1 : 4  % loop over edges
        edge = edges(ll,:);
        sctr_edge = sctr(edge);
        node_new = mean(node(sctr_edge,:)) ;
        
        if element_p2(kk,init+ll) == 0 
            node_p2 = [ node_p2 ; node_new ] ; 
            sc_new = size(node_p2,1);
            element_p2(kk,init+ll) = sc_new ; 
        
        
            for ii = 1 : length(c_elems)
                if c_elems(ii) == kk 
                    continue 
                end
                [yy,uu] = ismember(sctr_edge,element(c_elems(ii),:));
                if all( yy )
                    edge_ = unique([uu]) ;
                    [~,oo] = ismember(edge_,edges,'rows') ; 
                    if element_p2(c_elems(ii),init+oo) == 0 
                        element_p2(c_elems(ii),init+oo) = sc_new ; 
                    end
%                     pause
                end
            end
        end
    end

    
%     add extra node in the center
    node_new = [(node(sctr(1),1)+node(sctr(2),1))/2 (node(sctr(1),2)+node(sctr(4),2))/2];

        B_num = (pp-1)*2-1 ; % num new bubble funcs
        for kj = 1 : B_num 
            node_p2 = [ node_p2 ; node_new ] ; 
            sc_new = size(node_p2,1);
            element_p2(kk,init+4+kj) = sc_new ; 
        end


end

end