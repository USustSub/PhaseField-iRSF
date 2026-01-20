
function dphi = dlobatto ( p , xi ) 
        dphi = 0 ; 
        
        if p == 0 
           dphi = l0(xi) ;  
        elseif p == 1 
           dphi = l1(xi) ; 
        elseif p == 2 
           dphi = l2(xi) ; 
        elseif p == 3 
           dphi = l3(xi) ; 
        elseif p == 4 
           dphi = l4(xi) ; 
        elseif p == 5 
           dphi = l5(xi) ; 
        elseif p == 6 
           dphi = l6(xi) ; 
        elseif p == 7 
           dphi = l7(xi) ; 
       elseif p == 8 
           dphi = l8(xi) ;
        elseif p == 9 
           dphi = l9(xi) ; 
       elseif p == 10 
           dphi = l10(xi) ;
        end




function dphi = l0( xi ) 
    dphi = (-1)/2 ; 

function dphi = l1( xi ) 
    dphi = (+1)/2 ; 

function dphi = l2( xi ) 
    dphi = 0.5*sqrt(3/2)*(2*xi) ; 

function dphi = l3( xi ) 
    dphi = 0.5*sqrt(5/2)*(3*xi^2-1); 

function dphi = l4( xi ) 
    dphi = (1/8)*sqrt(7/2)*(10*xi*(xi^2 - 1) + 2*xi*(5*xi^2 - 1) ); 

function dphi = l5( xi ) 
%     dphi = (1/8)*sqrt(9/2)*(xi^2 - 1)*(7*xi^2 - 3) + 14*xi^2*(xi^2 - 1) + 2*xi^2*(7*xi^2 - 3) ; 
    dphi = (1/8)*sqrt(9/2)*(35*xi^4 - 30*xi^2 + 3); 

    function dphi = l6( xi ) 
%     dphi = (1/16)*sqrt(11/2)*2*xi*(21*xi^4 - 14*xi^2 + 1) - (xi^2 - 1)*(- 84*xi^3 + 28*xi) ; 
    dphi = (1/16)*sqrt(11/2)* ( 126*xi^5 - 140*xi^3 + 30*xi); 

function dphi = l7( xi ) 
%     dphi = (1/16)*sqrt(13/2)*2*xi^2*(33*xi^4 - 30*xi^2 + 5) + (xi^2 - 1)*(33*xi^4 - 30*xi^2 + 5) - xi*(xi^2 - 1)*(- 132*xi^3 + 60*xi) ; 
        dphi = (1/16)*sqrt(13/2)* ( 231*xi^6 - 150*xi^4 + 15*xi^2 - 165*xi^4 + 90*xi^2 - 5 );

function dphi = l8( xi ) 
%     dphi = (1/128)*sqrt(15/2)*(xi^2 - 1)*(2574*xi^5 - 1980*xi^3 + 270*xi) + 2*xi*(429*xi^6 - 495*xi^4 + 135*xi^2 - 5) ; 
    dphi = (1/128)*sqrt(15/2)*((2*xi)*(429*xi^6-495*xi^4+135*xi^2-5) + (xi^2-1)*(2574*xi^5-1980*xi^3+270*xi^1) ); 

function dphi = l9( xi ) 
%     dphi = (1/128)*sqrt(17/2)*2*xi^2*(715*xi^6 - 1001*xi^4 + 385*xi^2 - 35) + (xi^2 - 1)*(715*xi^6 - 1001*xi^4 + 385*xi^2 - 35) + xi*(xi^2 - 1)*(4290*xi^5 - 4004*xi^3 + 770*xi) ; 
    dphi = (1/128)*sqrt(17/2)*(2*xi)*(715*xi^6-1001*xi^4+385*xi^2-35)*xi+ ... 
            (1/128)*sqrt(17/2)*(xi^2-1)*(4290*xi^5-4004*xi^3+770*xi^1)*xi + ... 
            (1/128)*sqrt(17/2)*(xi^2-1)*(715*xi^6-1001*xi^4+385*xi^2-35)*1 ;

function dphi = l10( xi ) 
    dphi = (1/256)*sqrt(19/2)*2*xi*(2431*xi^8 - 4004*xi^6 + 2002*xi^4 - 308*xi^2 + 7) - (xi^2 - 1)*(- 19448*xi^7 + 24024*xi^5 - 8008*xi^3 + 616*xi); 

    