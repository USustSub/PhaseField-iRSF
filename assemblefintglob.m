function [fintglob] = assemblefintglob(fintglob,fint,ndf,con,ielm,lin_quad,interface)

if (lin_quad == 1)  % four-node elements (and four-node interfaces)
    
    fintglob(ndf*con(ielm,1)-1,1) =  fintglob(ndf*con(ielm,1)-1,1) - fint(1,1);
    fintglob(ndf*con(ielm,1)  ,1) =  fintglob(ndf*con(ielm,1)  ,1) - fint(2,1);
    
    fintglob(ndf*con(ielm,2)-1,1) =  fintglob(ndf*con(ielm,2)-1,1) - fint(3,1);
    fintglob(ndf*con(ielm,2)  ,1) =  fintglob(ndf*con(ielm,2)  ,1) - fint(4,1);
    
    fintglob(ndf*con(ielm,3)-1,1) =  fintglob(ndf*con(ielm,3)-1,1) - fint(5,1);
    fintglob(ndf*con(ielm,3)  ,1) =  fintglob(ndf*con(ielm,3)  ,1) - fint(6,1);
    
    fintglob(ndf*con(ielm,4)-1,1) =  fintglob(ndf*con(ielm,4)-1,1) - fint(7,1);
    fintglob(ndf*con(ielm,4)  ,1) =  fintglob(ndf*con(ielm,4)  ,1) - fint(8,1);
        
else
    
    if (interface == 0)  % eight-node elements
        
    fintglob(ndf*con(ielm,1)-1,1) =  fintglob(ndf*con(ielm,1)-1,1) - fint(1,1);
    fintglob(ndf*con(ielm,1)  ,1) =  fintglob(ndf*con(ielm,1)  ,1) - fint(2,1);
    
    fintglob(ndf*con(ielm,2)-1,1) =  fintglob(ndf*con(ielm,2)-1,1) - fint(3,1);
    fintglob(ndf*con(ielm,2)  ,1) =  fintglob(ndf*con(ielm,2)  ,1) - fint(4,1);
    
    fintglob(ndf*con(ielm,3)-1,1) =  fintglob(ndf*con(ielm,3)-1,1) - fint(5,1);
    fintglob(ndf*con(ielm,3)  ,1) =  fintglob(ndf*con(ielm,3)  ,1) - fint(6,1);
    
    fintglob(ndf*con(ielm,4)-1,1) =  fintglob(ndf*con(ielm,4)-1,1) - fint(7,1);
    fintglob(ndf*con(ielm,4)  ,1) =  fintglob(ndf*con(ielm,4)  ,1) - fint(8,1);
    
    fintglob(ndf*con(ielm,5)-1,1) =  fintglob(ndf*con(ielm,5)-1,1) - fint(9,1);
    fintglob(ndf*con(ielm,5)  ,1) =  fintglob(ndf*con(ielm,5)  ,1) - fint(10,1);
    
    fintglob(ndf*con(ielm,6)-1,1) =  fintglob(ndf*con(ielm,6)-1,1) - fint(11,1);
    fintglob(ndf*con(ielm,6)  ,1) =  fintglob(ndf*con(ielm,6)  ,1) - fint(12,1);
        
    fintglob(ndf*con(ielm,7)-1,1) =  fintglob(ndf*con(ielm,7)-1,1) - fint(13,1);
    fintglob(ndf*con(ielm,7)  ,1) =  fintglob(ndf*con(ielm,7)  ,1) - fint(14,1);
    
    fintglob(ndf*con(ielm,8)-1,1) =  fintglob(ndf*con(ielm,8)-1,1) - fint(15,1);
    fintglob(ndf*con(ielm,8)  ,1) =  fintglob(ndf*con(ielm,8)  ,1) - fint(16,1);
     
    else  % six-node interfaces
     
    fintglob(ndf*con(ielm,1)-1,1) =  fintglob(ndf*con(ielm,1)-1,1) - fint(1,1);
    fintglob(ndf*con(ielm,1)  ,1) =  fintglob(ndf*con(ielm,1)  ,1) - fint(2,1);
    
    fintglob(ndf*con(ielm,2)-1,1) =  fintglob(ndf*con(ielm,2)-1,1) - fint(3,1);
    fintglob(ndf*con(ielm,2)  ,1) =  fintglob(ndf*con(ielm,2)  ,1) - fint(4,1);
    
    fintglob(ndf*con(ielm,3)-1,1) =  fintglob(ndf*con(ielm,3)-1,1) - fint(5,1);
    fintglob(ndf*con(ielm,3)  ,1) =  fintglob(ndf*con(ielm,3)  ,1) - fint(6,1);
    
    fintglob(ndf*con(ielm,4)-1,1) =  fintglob(ndf*con(ielm,4)-1,1) - fint(7,1);
    fintglob(ndf*con(ielm,4)  ,1) =  fintglob(ndf*con(ielm,4)  ,1) - fint(8,1);
    
    fintglob(ndf*con(ielm,5)-1,1) =  fintglob(ndf*con(ielm,5)-1,1) - fint(9,1);
    fintglob(ndf*con(ielm,5)  ,1) =  fintglob(ndf*con(ielm,5)  ,1) - fint(10,1);
    
    fintglob(ndf*con(ielm,6)-1,1) =  fintglob(ndf*con(ielm,6)-1,1) - fint(11,1);
    fintglob(ndf*con(ielm,6)  ,1) =  fintglob(ndf*con(ielm,6)  ,1) - fint(12,1);       
        
    end
    
end