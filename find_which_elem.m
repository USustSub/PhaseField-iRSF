function [iel_,sc_iint,local] = find_which_elem ( xx_ , element , node )
iel_ = -1 ; 
for ij_ = 1 : size(element,1)
    
   n_sc = node(element(ij_,:),:) ;
   min_x_sc = min(n_sc(:,1));
   max_x_sc = max(n_sc(:,1));
   min_y_sc = min(n_sc(:,2));
   max_y_sc = max(n_sc(:,2));
   
   if xx_(1)>=min_x_sc && xx_(1)<=max_x_sc && xx_(2)>=min_y_sc && xx_(2)<=max_y_sc
       iel_ = ij_ ;
       break;
   end
    
end

    sc_iint = element(iel_,:);

     CC.type = 'Q4';
    [local] = element_coordinate (xx_ , sc_iint , node , CC);

