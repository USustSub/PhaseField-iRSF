
function phi = lobatto ( p , xi ) 
phi = 0 ; 
        if p == 0 
           phi = l0(xi) ;  
        elseif p == 1 
           phi = l1(xi) ; 
        elseif p == 2 
           phi = l2(xi) ; 
        elseif p == 3 
           phi = l3(xi) ; 
        elseif p == 4 
           phi = l4(xi) ; 
        elseif p == 5 
           phi = l5(xi) ; 
        elseif p == 6 
           phi = l6(xi) ; 
        elseif p == 7 
           phi = l7(xi) ; 
       elseif p == 8 
           phi = l8(xi) ;
        elseif p == 9 
           phi = l9(xi) ; 
       elseif p == 10 
           phi = l10(xi) ;
        end




function phi = l0( xi ) 
    phi = (1-xi)/2 ; 

function phi = l1( xi ) 
    phi = (1+xi)/2 ; 

function phi = l2( xi ) 
    phi = (1/2)*sqrt(3/2)*(xi^2-1) ; 

function phi = l3( xi ) 
    phi = (1/2)*sqrt(5/2)*(xi^2-1)*xi ; 

function phi = l4( xi ) 
    phi = (1/8)*sqrt(7/2)*(xi^2-1)*(5*xi^2-1) ; 

function phi = l5( xi ) 
    phi = (1/8)*sqrt(9/2)*(xi^2-1)*(7*xi^2-3)*xi ; 
 
function phi = l6( xi ) 
    phi = (1/16)*sqrt(11/2)*(xi^2-1)*(21*xi^4-14*xi^2+1) ; 

function phi = l7( xi ) 
    phi = (1/16)*sqrt(13/2)*(xi^2-1)*(33*xi^4-30*xi^2+5)*xi ; 
        
function phi = l8( xi ) 
    phi = (1/128)*sqrt(15/2)*(xi^2-1)*(429*xi^6-495*xi^4+135*xi^2-5) ; 

function phi = l9( xi ) 
    phi = (1/128)*sqrt(17/2)*(xi^2-1)*(715*xi^6-1001*xi^4+385*xi^2-35)*xi ; 

function phi = l10( xi ) 
    phi = (1/256)*sqrt(19/2)*(xi^2-1)*(2431*xi^8-4004*xi^6+2002*xi^4-308*xi^2+7) ; 

    
    
