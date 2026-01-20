function [conn_elems ] = find_conn_elements (node , element ,conn_nodes) 

    conn_elems = []; 
        cc = zeros(size(element,1),1) ; 
    for jj = 1 : size(element,1)
    %     jj
        sctr = element(jj,:) ; 
    %     if ( sctr(1)==sctr(2)) 
    %         continue
    %     end
        conn_elems_tmp = [ ]; 
        for kk = 1 : length(sctr) 
            for nn = 1 : conn_nodes(sctr(kk),1)
            cc(jj) = cc(jj)+1 ;             
    %         conn_elems(jj,cc(jj)+1) = conn_nodes2(sctr(kk),nn+1); 
    %         conn_elems(jj,1) = cc(jj); 
            conn_elems_tmp = [conn_elems_tmp ; conn_nodes(sctr(kk),nn+1)]; 
            end
        end
        conn_elems_tmp = unique(conn_elems_tmp);
    %     conn_elems_tmp = setdiff(conn_elems_tmp,jj); 
    %     conn_elems_tmp
        conn_elems{jj} =[length(conn_elems_tmp) conn_elems_tmp']  ; 
    %  plot_mesh(node,element( conn_elems{jj}(2:end),:),'Q4',[rand(1,3)])
    %   plot_mesh(node,element( jj,:),'Q4','b')

    %     pause
    end