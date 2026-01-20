
    topNodes = find(abs(node(:,2)-max(node(:,2)))<0.001) ; 
    botNodes = find(abs(node(:,2)-min(node(:,2)))<0.001) ; 
    rightNodes = find(abs(node(:,1)-max(node(:,1)))<0.001) ; 
    leftNodes = find(abs(node(:,1)-min(node(:,1)))<0.001) ; 
    
    [c,I] = sortrows(node(topNodes,1),1) ;
    topNodes = topNodes(I ) ; 
    [c,I] = sortrows(node(botNodes,1),1) ;
    botNodes = botNodes(I ) ; 
    [c,I] = sortrows(node(rightNodes,2),1) ;
    rightNodes = rightNodes(I ) ; 
    [c,I] = sortrows(node(leftNodes,2),1) ;
    leftNodes = leftNodes(I ) ; 
    
    uln = find(node(:,1) == min(node(:,1)) & (node(:,2) == max(node(:,2))));    
    urn = find(node(:,1) == max(node(:,1)) & (node(:,2) == max(node(:,2))));    
    lln = find(node(:,1) == min(node(:,1)) & (node(:,2) == min(node(:,2))));
    lrn = find(node(:,1) == max(node(:,1)) & (node(:,2) == min(node(:,2)))); 
