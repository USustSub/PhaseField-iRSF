function [jump] = jumpmatrix(nmat,rotmat,uloc)

jump = zeros(2,1,size(nmat,3));
    
for IP = 1 : size(nmat,3)
    
    jump(:,:,IP) = rotmat*nmat(:,:,IP)*uloc; % horizontal jump / vertical jump    
    
end