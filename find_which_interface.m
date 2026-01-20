function iint = find_which_interface ( xgp , node , conint_cohesive )

for ii = 1 : size(conint_cohesive,1)
    
    node_sc = node(conint_cohesive(ii,:),:);
    min_ = min(node_sc(:,1));
    max_ = max(node_sc(:,1));
    if xgp(1)>=min_ && xgp(1)<=max_
        iint = ii ; 
        return ;
    end
end