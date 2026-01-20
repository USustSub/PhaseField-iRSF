function [kmat] = kmatrixint(nmat,tmat,rotmat,detJ_IP,b_b)
global LINE_
kmat = zeros(size(nmat,2),size(nmat,2));

for IP = 1:size(detJ_IP,1)
           
    kmat = kmat + nmat(:,:,IP)'*(rotmat'*tmat(:,:,IP)*rotmat)*nmat(:,:,IP)*detJ_IP(IP,1)*detJ_IP(IP,2)*b_b;
    LINE_ = LINE_ + detJ_IP(IP,1)*detJ_IP(IP,2)*b_b ;
end 
