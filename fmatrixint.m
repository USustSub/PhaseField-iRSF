function [fint] = fmatrixint(nmat,traction,rotmat,detJ_IP,b_b)

fint = zeros(size(nmat,2),1);

for IP = 1:size(detJ_IP,1)
   
    fint = fint + nmat(:,:,IP)'*(rotmat'*traction(:,:,IP))*detJ_IP(IP,1)*detJ_IP(IP,2)*b_b;
    
end