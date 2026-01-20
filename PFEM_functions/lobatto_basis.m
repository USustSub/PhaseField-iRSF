% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lobatto basis functions
% Mohsen Goudarzi September 2016 - Delft 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nv,dN]=lobatto_basis(type,coord,dim,order)

if ( nargin == 2 )
    dim=1;
end
       xi = coord (1) ; 
       eta = coord(2) ;
%     phi = l ( p , xi )  

    if length(coord) == 3 
        zeta = coord(3) ; 
    end

switch type
   case 'Q4'
       
       
    if (order == 1) 
        N = zeros(4,1) ; 
        dN = zeros(4,2) ; 
    elseif (order == 2) 
        N = zeros(9,1) ; 
        dN = zeros(9,2) ; 
    elseif (order == 3) 
        N = zeros(16,1) ; 
        dN = zeros(16,2) ; 
    elseif (order == 4) 
        N = zeros(25,1) ; 
        dN = zeros(25,2) ; 
    elseif (order == 5) 
        N = zeros(36,1) ; 
        dN = zeros(36,2) ;
    elseif (order == 6 )
        N = zeros(49,1) ;
        dN = zeros(49,2) ;
    elseif (order == 7)
        N = zeros(64,1) ; 
        dN = zeros(64,2) ; 
    end       
       
       
       
       if order >= 1 
%         N = [] ; 
        N(1,1) = lobatto(0,xi)*lobatto(0,eta) ;
        N(2,1) = lobatto(1,xi)*lobatto(0,eta) ;
        N(3,1) = lobatto(1,xi)*lobatto(1,eta) ;
        N(4,1) = lobatto(0,xi)*lobatto(1,eta) ;

        dN(1,1) =  dlobatto(0,xi)*lobatto(0,eta) ;
        dN(1,2) =  lobatto(0,xi)*dlobatto(0,eta) ;
        dN(2,1) =  dlobatto(1,xi)*lobatto(0,eta) ;
        dN(2,2) =  lobatto(1,xi)*dlobatto(0,eta) ;
        dN(3,1) =  dlobatto(1,xi)*lobatto(1,eta) ;
        dN(3,2) =  lobatto(1,xi)*dlobatto(1,eta) ;
        dN(4,1) =  dlobatto(0,xi)*lobatto(1,eta) ;
        dN(4,2) =  lobatto(0,xi)*dlobatto(1,eta) ;
       end
       if order >= 2 
%     case 2
%         N = [] ; 
%         N(1,1) = lobatto(0,xi)*lobatto(0,eta) ;
%         N(2,1) = lobatto(1,xi)*lobatto(0,eta) ;
%         N(3,1) = lobatto(1,xi)*lobatto(1,eta) ;
%         N(4,1) = lobatto(0,xi)*lobatto(1,eta) ;
% 
%         dN(1,1) =  dlobatto(0,xi)*lobatto(0,eta) ;
%         dN(1,2) =  lobatto(0,xi)*dlobatto(0,eta) ;
%         dN(2,1) =  dlobatto(1,xi)*lobatto(0,eta) ;
%         dN(2,2) =  lobatto(1,xi)*dlobatto(0,eta) ;
%         dN(3,1) =  dlobatto(1,xi)*lobatto(1,eta) ;
%         dN(3,2) =  lobatto(1,xi)*dlobatto(1,eta) ;
%         dN(4,1) =  dlobatto(0,xi)*lobatto(1,eta) ;
%         dN(4,2) =  lobatto(0,xi)*dlobatto(1,eta) ;

        k = 2 ; 
        % edge functions
        N(5,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(6,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(7,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(8,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(5,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(5,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(6,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(6,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(7,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(7,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(8,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(8,2) =  lobatto (k,xi)*dlobatto(1,eta) ;

% bubble function
        N(9,1) = lobatto(2,xi)*lobatto(2,eta) ;
        dN(9,1) = dlobatto(2,xi)*lobatto(2,eta) ;
        dN(9,2) = lobatto(2,xi)*dlobatto(2,eta) ;
%     case 3
        end
        if order >= 3 

        k = 3 ; 
% edge functions
        N(10,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(11,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(12,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(13,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(10,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(10,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(11,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(11,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(12,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(12,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(13,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(13,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        N(14,1) = lobatto(2,xi)*lobatto(3,eta) ;
        dN(14,1) = dlobatto(2,xi)*lobatto(3,eta) ;
        dN(14,2) = lobatto(2,xi)*dlobatto(3,eta) ;
        
        N(15,1) = lobatto(3,xi)*lobatto(2,eta) ;
        dN(15,1) = dlobatto(3,xi)*lobatto(2,eta) ;
        dN(15,2) = lobatto(3,xi)*dlobatto(2,eta) ;
                
        N(16,1) = lobatto(3,xi)*lobatto(3,eta) ;
        dN(16,1) = dlobatto(3,xi)*lobatto(3,eta) ;
        dN(16,2) = lobatto(3,xi)*dlobatto(3,eta) ;
        end
        if order >= 4 
        k = 4 ; 
% edge functions
        N(17,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(18,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(19,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(20,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(17,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(17,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(18,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(18,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(19,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(19,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(20,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(20,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        N(21,1) = lobatto(2,xi)*lobatto(4,eta) ;
        dN(21,1) = dlobatto(2,xi)*lobatto(4,eta) ;
        dN(21,2) = lobatto(2,xi)*dlobatto(4,eta) ;
        
        N(22,1) = lobatto(4,xi)*lobatto(2,eta) ;
        dN(22,1) = dlobatto(4,xi)*lobatto(2,eta) ;
        dN(22,2) = lobatto(4,xi)*dlobatto(2,eta) ;
                
        N(23,1) = lobatto(3,xi)*lobatto(4,eta) ;
        dN(23,1) = dlobatto(3,xi)*lobatto(4,eta) ;
        dN(23,2) = lobatto(3,xi)*dlobatto(4,eta) ;
        
        N(24,1) = lobatto(4,xi)*lobatto(3,eta) ;
        dN(24,1) = dlobatto(4,xi)*lobatto(3,eta) ;
        dN(24,2) = lobatto(4,xi)*dlobatto(3,eta) ;
        
        N(25,1) = lobatto(4,xi)*lobatto(4,eta) ;
        dN(25,1) = dlobatto(4,xi)*lobatto(4,eta) ;
        dN(25,2) = lobatto(4,xi)*dlobatto(4,eta) ;
        
        end
        
    if order >= 5 
        k = 5 ; 
% edge functions
        N(26,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(27,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(28,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(29,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(26,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(26,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(27,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(27,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(28,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(28,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(29,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(29,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        N(30,1) = lobatto(2,xi)*lobatto(5,eta) ;
        dN(30,1) = dlobatto(2,xi)*lobatto(5,eta) ;
        dN(30,2) = lobatto(2,xi)*dlobatto(5,eta) ;
        
        N(31,1) = lobatto(5,xi)*lobatto(2,eta) ;
        dN(31,1) = dlobatto(5,xi)*lobatto(2,eta) ;
        dN(31,2) = lobatto(5,xi)*dlobatto(2,eta) ;
                
        N(32,1) = lobatto(3,xi)*lobatto(5,eta) ;
        dN(32,1) = dlobatto(3,xi)*lobatto(5,eta) ;
        dN(32,2) = lobatto(3,xi)*dlobatto(5,eta) ;
        
        N(33,1) = lobatto(5,xi)*lobatto(3,eta) ;
        dN(33,1) = dlobatto(5,xi)*lobatto(3,eta) ;
        dN(33,2) = lobatto(5,xi)*dlobatto(3,eta) ;
        
        N(34,1) = lobatto(4,xi)*lobatto(5,eta) ;
        dN(34,1) = dlobatto(4,xi)*lobatto(5,eta) ;
        dN(34,2) = lobatto(4,xi)*dlobatto(5,eta) ;
        
        N(35,1) = lobatto(5,xi)*lobatto(4,eta) ;
        dN(35,1) = dlobatto(5,xi)*lobatto(4,eta) ;
        dN(35,2) = lobatto(5,xi)*dlobatto(4,eta) ;
        
        N(36,1) = lobatto(5,xi)*lobatto(5,eta) ;
        dN(36,1) = dlobatto(5,xi)*lobatto(5,eta) ;
        dN(36,2) = lobatto(5,xi)*dlobatto(5,eta) ;
        
    end
        
    if order >= 6 
        k = 6 ; 
% edge functions
        N(37,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(38,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(39,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(40,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(37,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(37,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(38,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(38,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(39,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(39,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(40,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(40,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        N(41,1) = lobatto(2,xi)*lobatto(6,eta) ;
        dN(41,1) = dlobatto(2,xi)*lobatto(6,eta) ;
        dN(41,2) = lobatto(2,xi)*dlobatto(6,eta) ;
        
        N(42,1) = lobatto(6,xi)*lobatto(2,eta) ;
        dN(42,1) = dlobatto(6,xi)*lobatto(2,eta) ;
        dN(42,2) = lobatto(6,xi)*dlobatto(2,eta) ;
                
        N(43,1) = lobatto(3,xi)*lobatto(6,eta) ;
        dN(43,1) = dlobatto(3,xi)*lobatto(6,eta) ;
        dN(43,2) = lobatto(3,xi)*dlobatto(6,eta) ;
        
        N(44,1) = lobatto(6,xi)*lobatto(3,eta) ;
        dN(44,1) = dlobatto(6,xi)*lobatto(3,eta) ;
        dN(44,2) = lobatto(6,xi)*dlobatto(3,eta) ;
        
        N(45,1) = lobatto(4,xi)*lobatto(6,eta) ;
        dN(45,1) = dlobatto(4,xi)*lobatto(6,eta) ;
        dN(45,2) = lobatto(4,xi)*dlobatto(6,eta) ;
        
        N(46,1) = lobatto(6,xi)*lobatto(4,eta) ;
        dN(46,1) = dlobatto(6,xi)*lobatto(4,eta) ;
        dN(46,2) = lobatto(6,xi)*dlobatto(4,eta) ;
        
        N(47,1) = lobatto(5,xi)*lobatto(6,eta) ;
        dN(47,1) = dlobatto(5,xi)*lobatto(6,eta) ;
        dN(47,2) = lobatto(5,xi)*dlobatto(6,eta) ;
        
        N(48,1) = lobatto(6,xi)*lobatto(5,eta) ;
        dN(48,1) = dlobatto(6,xi)*lobatto(5,eta) ;
        dN(48,2) = lobatto(6,xi)*dlobatto(5,eta) ;
        
        N(49,1) = lobatto(6,xi)*lobatto(6,eta) ;
        dN(49,1) = dlobatto(6,xi)*lobatto(6,eta) ;
        dN(49,2) = lobatto(6,xi)*dlobatto(6,eta) ;
        
    end
    
    if order >= 7 
        k = 7 ; 
% edge functions
        N(50,1) = lobatto(0,xi)*lobatto(k,eta) ;
        N(51,1) = lobatto(1,xi)*lobatto(k,eta) ;
        N(52,1) = lobatto(k,xi)*lobatto(0,eta) ;
        N(53,1) = lobatto(k,xi)*lobatto(1,eta) ;

        dN(50,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        dN(50,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        dN(51,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        dN(51,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        dN(52,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        dN(52,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        dN(53,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        dN(53,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        N(54,1) = lobatto(2,xi)*lobatto(7,eta) ;
        dN(54,1) = dlobatto(2,xi)*lobatto(7,eta) ;
        dN(54,2) = lobatto(2,xi)*dlobatto(7,eta) ;

        N(55,1) = lobatto(7,xi)*lobatto(2,eta) ;
        dN(55,1) = dlobatto(7,xi)*lobatto(2,eta) ;
        dN(55,2) = lobatto(7,xi)*dlobatto(2,eta) ;

        N(56,1) = lobatto(3,xi)*lobatto(7,eta) ;
        dN(56,1) = dlobatto(3,xi)*lobatto(7,eta) ;
        dN(56,2) = lobatto(3,xi)*dlobatto(7,eta) ;

        N(57,1) = lobatto(7,xi)*lobatto(3,eta) ;
        dN(57,1) = dlobatto(7,xi)*lobatto(3,eta) ;
        dN(57,2) = lobatto(7,xi)*dlobatto(3,eta) ;

        N(58,1) = lobatto(4,xi)*lobatto(7,eta) ;
        dN(58,1) = dlobatto(4,xi)*lobatto(7,eta) ;
        dN(58,2) = lobatto(4,xi)*dlobatto(7,eta) ;
        
        N(59,1) = lobatto(7,xi)*lobatto(4,eta) ;
        dN(59,1) = dlobatto(7,xi)*lobatto(4,eta) ;
        dN(59,2) = lobatto(7,xi)*dlobatto(4,eta) ;
        
        N(60,1) = lobatto(5,xi)*lobatto(7,eta) ;
        dN(60,1) = dlobatto(5,xi)*lobatto(7,eta) ;
        dN(60,2) = lobatto(5,xi)*dlobatto(7,eta) ;
        
        N(61,1) = lobatto(7,xi)*lobatto(5,eta) ;
        dN(61,1) = dlobatto(7,xi)*lobatto(5,eta) ;
        dN(61,2) = lobatto(7,xi)*dlobatto(5,eta) ;
        
        N(62,1) = lobatto(6,xi)*lobatto(7,eta) ;
        dN(62,1) = dlobatto(6,xi)*lobatto(7,eta) ;
        dN(62,2) = lobatto(6,xi)*dlobatto(7,eta) ;
        
        N(63,1) = lobatto(7,xi)*lobatto(6,eta) ;
        dN(63,1) = dlobatto(7,xi)*lobatto(6,eta) ;
        dN(63,2) = lobatto(7,xi)*dlobatto(6,eta) ;
        
        N(64,1) = lobatto(7,xi)*lobatto(7,eta) ;
        dN(64,1) = dlobatto(7,xi)*lobatto(7,eta) ;
        dN(64,2) = lobatto(7,xi)*dlobatto(7,eta) ;
        
    end
    
    case 'H8'
       if order >= 1 
        N = [] ; 
        N(1,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(2,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(3,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(4,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(5,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(6,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(7,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(8,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        
        
        dN(1,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(1,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(1,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(2,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(2,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(2,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;
        
        dN(3,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(3,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(3,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;
        
        dN(4,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(4,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(4,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;
       
        dN(5,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(5,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(5,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(6,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(6,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(6,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(7,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(7,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(7,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
       
        dN(8,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(8,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(8,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
       
       end
        if order >= 2 

        k = 2 ; 
        % edge functions
        N(9,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(10,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(11,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(12,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(13,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(14,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(15,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(16,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(17,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(18,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(19,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(20,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(9,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(9,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(9,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(10,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(10,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(10,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(11,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(11,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(11,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(12,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(12,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(12,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(13,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(13,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(13,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(14,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(14,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(14,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(15,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(15,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(15,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(16,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(16,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(16,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(17,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(17,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(17,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(18,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(18,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(18,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(19,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(19,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(19,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(20,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(20,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(20,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 2 ; 
        N(21,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(22,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(23,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(24,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(25,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(26,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(21,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(21,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(21,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(22,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(22,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(22,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(23,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(23,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(23,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(24,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(24,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(24,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(25,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(25,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(25,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(26,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(26,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(26,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

% bubble function
        N(27,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dN(27,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(27,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(27,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

%     case 3
        end
        
        if order >= 3 

        k = 3 ; 
        % edge functions
        N(28,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(29,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(30,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(31,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(32,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(33,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(34,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(35,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(36,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(37,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(38,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(39,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(28,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(28,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(28,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(29,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(29,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(29,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(30,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(30,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(30,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(31,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(31,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(31,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(32,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(32,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(32,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(33,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(33,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(33,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(34,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(34,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(34,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(35,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(35,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(35,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(36,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(36,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(36,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(37,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(37,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(37,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(38,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(38,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(38,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(39,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(39,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(39,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 3 ; 
        N(40,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(41,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(42,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(43,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(44,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(45,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(40,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(40,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(40,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(41,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(41,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(41,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(42,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(42,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(42,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(43,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(43,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(43,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(44,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(44,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(44,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(45,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(45,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(45,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 3 ; n2 = 2 ; 
        N(46,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(47,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(48,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(49,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(50,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(51,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(46,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(46,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(46,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(47,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(47,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(47,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(48,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(48,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(48,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(49,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(49,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(49,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(50,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(50,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(50,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(51,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(51,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(51,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 3 ; 
        N(52,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(53,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(54,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(55,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(56,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(57,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(52,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(52,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(52,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(53,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(53,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(53,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(54,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(54,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(54,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(55,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(55,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(55,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(56,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(56,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(56,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(57,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(57,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(57,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function
        N(58,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        N(59,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        N(60,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        N(61,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        N(62,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        N(63,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        N(64,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        
        dN(58,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(58,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(58,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;

        dN(59,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(59,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(59,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;

        dN(60,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(60,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(60,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

        dN(61,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(61,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(61,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;

        dN(62,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(62,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(62,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;

        dN(63,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(63,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(63,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

        dN(64,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(64,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(64,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

        end
        
        if order >= 4
        k = 4 ; 
        % edge functions
        N(65,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(66,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(67,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(68,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(69,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(70,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(71,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(72,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(73,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(74,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(75,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(76,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(65,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(65,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(65,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(66,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(66,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(66,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(67,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(67,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(67,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(68,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(68,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(68,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(69,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(69,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(69,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(70,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(70,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(70,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(71,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(71,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(71,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(72,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(72,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(72,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(73,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(73,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(73,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(74,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(74,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(74,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(75,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(75,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(75,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(76,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(76,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(76,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 4 ; 
        N(77,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(78,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(79,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(80,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(81,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(82,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(77,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(77,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(77,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(78,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(78,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(78,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(79,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(79,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(79,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(80,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(80,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(80,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(81,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(81,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(81,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(82,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(82,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(82,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 4 ; n2 = 2 ; 
        N(83,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(84,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(85,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(86,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(87,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(88,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(83,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(83,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(83,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(84,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(84,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(84,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(85,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(85,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(85,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(86,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(86,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(86,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(87,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(87,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(87,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(88,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(88,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(88,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 3 ; n2 = 4 ; 
        N(89,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(90,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(91,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(92,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(93,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(94,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(89,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(89,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(89,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(90,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(90,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(90,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(91,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(91,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(91,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(92,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(92,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(92,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(93,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(93,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(93,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(94,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(94,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(94,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 4 ; n2 = 3 ; 
        N(95,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(96,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(97,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(98,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(99,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(100,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(95,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(95,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(95,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(96,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(96,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(96,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(97,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(97,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(97,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(98,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(98,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(98,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(99,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(99,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(99,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(100,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(100,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(100,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           

        n1 = 4 ; n2 = 4 ; 
        N(101,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(102,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(103,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(104,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(105,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(106,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(101,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(101,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(101,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(102,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(102,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(102,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(103,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(103,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(103,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(104,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(104,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(104,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(105,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(105,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(105,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(106,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(106,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(106,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function
% 2,3,4
% 2 2 
        N(107,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(107,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(107,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(107,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;


        N(108,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dN(108,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(108,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(108,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;


        N(109,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dN(109,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(109,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(109,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

        N(110,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dN(110,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(110,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(110,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;

        N(111,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(111,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(111,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(111,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;

        N(112,1)  =  lobatto(2,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(112,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(112,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(112,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;

        N(113,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dN(113,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(113,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(113,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;


        N(114,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(114,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(114,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(114,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;

        N(115,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dN(115,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(115,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(115,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

        N(116,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(116,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(116,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(116,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;

        N(117,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dN(117,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(117,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(117,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;

        N(118,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(118,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(118,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(118,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;

        N(119,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(119,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(119,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(119,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
        
        N(120,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(4,zeta) ;        
        dN(120,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(120,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(120,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
        
        N(121,1)  = lobatto(3,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(121,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(121,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(121,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;


        N(122,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(2,zeta) ;        
        dN(122,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(122,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(122,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
        
        N(123,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dN(123,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(123,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(123,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
        

        N(124,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(124,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(124,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(124,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
        
        N(125,1)  = lobatto(3,xi)*lobatto(4,eta)*lobatto(2,zeta) ;        
        dN(125,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(125,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(125,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
        end
        
        
        if order >= 5
        k = 5 ; 
        % edge functions
        N(126,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(127,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(128,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(129,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(130,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(131,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(132,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(133,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(134,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(135,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(136,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(137,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(126,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(126,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(126,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(127,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(127,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(127,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(128,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(128,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(128,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(129,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(129,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(129,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(130,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(130,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(130,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(131,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(131,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(131,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(132,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(132,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(132,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(133,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(133,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(133,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(134,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(134,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(134,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(135,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(135,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(135,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(136,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(136,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(136,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(137,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(137,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(137,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 5 ; 
        N(138,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(139,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(140,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(141,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(142,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(143,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(138,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(138,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(138,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(139,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(139,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(139,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(140,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(140,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(140,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(141,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(141,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(141,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(142,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(142,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(142,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(143,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(143,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(143,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 5 ; n2 = 2 ; 
        N(144,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(145,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(146,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(147,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(148,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(149,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(144,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(144,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(144,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(145,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(145,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(145,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(146,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(146,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(146,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(147,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(147,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(147,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(148,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(148,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(148,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(149,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(149,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(149,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 5 ; n2 = 3 ; 
        N(150,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(151,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(152,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(153,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(154,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(155,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(150,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(150,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(150,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(151,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(151,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(151,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(152,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(152,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(152,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(153,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(153,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(153,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(154,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(154,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(154,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(155,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(155,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(155,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 5 ; 
        N(156,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(157,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(158,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(159,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(160,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(161,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(156,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(156,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(156,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(157,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(157,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(157,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(158,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(158,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(158,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(159,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(159,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(159,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(160,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(160,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(160,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(161,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(161,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(161,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
       
        n1 = 5 ; n2 = 4 ; 
        N(162,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(163,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(164,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(165,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(166,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(167,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(162,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(162,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(162,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(163,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(163,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(163,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(164,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(164,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(164,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(165,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(165,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(165,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(166,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(166,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(166,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(167,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(167,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(167,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         
        n1 = 4 ; n2 = 5 ; 
        N(168,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(169,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(170,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(171,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(172,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(173,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(168,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(168,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(168,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(169,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(169,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(169,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(170,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(170,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(170,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(171,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(171,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(171,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(172,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(172,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(172,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(173,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(173,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(173,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         



        n1 = 5 ; n2 = 5 ; 
        N(174,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(175,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(176,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(177,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(178,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(179,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(174,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(174,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(174,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(175,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(175,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(175,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(176,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(176,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(176,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(177,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(177,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(177,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(178,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(178,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(178,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(179,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(179,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(179,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function (needs to be changed)
% 2,3,4,5
% 225
        N(180,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(180,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(180,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(180,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;

% 335
        N(181,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(181,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(181,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(181,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 445
        N(182,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(182,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(182,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(182,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 235
        N(183,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(183,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(183,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(183,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 325
        N(184,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(184,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(184,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(184,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 245
        N(185,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(185,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(185,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(185,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 425
        N(186,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(186,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(186,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(186,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 345
        N(187,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(187,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(187,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(187,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 435
        N(188,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(188,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(188,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(188,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 522
        N(189,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dN(189,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(189,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(189,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;
% 533
        N(190,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dN(190,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(190,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(190,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;
% 544
        N(191,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(191,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(191,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(191,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
% 523
        N(192,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dN(192,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(192,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(192,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
% 532
        N(193,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        dN(193,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(193,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(193,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
% 524
        N(194,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(194,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(194,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(194,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;
% 542
        N(195,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dN(195,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(195,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(195,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
% 534
        N(196,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dN(196,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(196,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(196,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
% 543
        N(197,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(197,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(197,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(197,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
% 252
        N(198,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(198,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(198,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(198,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 353
        N(199,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(199,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(199,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(199,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 454
        N(200,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(200,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(200,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(200,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 253
        N(201,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(201,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(201,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(201,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 352
        N(202,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(202,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(202,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(202,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 254
        N(203,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(203,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(203,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(203,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 452
        N(204,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(204,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(204,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(204,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 354
        N(205,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(205,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(205,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(205,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 453
        N(206,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(206,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(206,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(206,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 255
        N(207,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(207,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(207,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(207,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 355
        N(208,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(208,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(208,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(208,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 455
        N(209,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(209,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(209,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(209,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 525
        N(210,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(210,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(210,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(210,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 535
        N(211,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(211,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(211,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(211,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 545
        N(212,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(212,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(212,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(212,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 552
        N(213,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(213,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(213,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(213,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 553
        N(214,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(214,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(214,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(214,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 554
        N(215,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(215,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(215,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(215,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 555
        N(216,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(216,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(216,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(216,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;

        end
        
        
        if order >= 6
        k = 6 ; 
        % edge functions
        N(217,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(218,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(219,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(220,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(221,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(222,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(223,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(224,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(225,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(226,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(227,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(228,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(217,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(217,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(217,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(218,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(218,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(218,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(219,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(219,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(219,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(220,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(220,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(220,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(221,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(221,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(221,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(222,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(222,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(222,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(223,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(223,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(223,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(224,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(224,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(224,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(225,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(225,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(225,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(226,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(226,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(226,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(227,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(227,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(227,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(228,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(228,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(228,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 6 ; 
        N(229,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(230,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(231,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(232,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(233,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(234,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(229,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(229,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(229,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(230,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(230,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(230,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(231,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(231,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(231,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(232,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(232,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(232,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(233,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(233,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(233,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(234,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(234,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(234,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 6 ; n2 = 2 ; 
        N(235,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(236,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(237,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(238,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(239,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(240,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(235,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(235,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(235,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(236,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(236,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(236,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(237,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(237,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(237,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(238,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(238,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(238,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(239,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(239,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(239,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(240,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(240,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(240,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 6 ; n2 = 3 ; 
        N(241,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(242,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(243,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(244,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(245,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(246,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(241,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(241,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(241,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(242,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(242,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(242,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(243,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(243,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(243,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(244,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(244,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(244,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(245,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(245,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(245,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(246,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(246,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(246,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 6 ; 
        N(247,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(248,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(249,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(250,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(251,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(252,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(247,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(247,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(247,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(248,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(248,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(248,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(249,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(249,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(249,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(250,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(250,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(250,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(251,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(251,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(251,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(252,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(252,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(252,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
       
        n1 = 6 ; n2 = 4 ; 
        N(253,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(254,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(255,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(256,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(257,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(258,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(253,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(253,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(253,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(254,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(254,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(254,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(255,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(255,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(255,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(256,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(256,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(256,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(257,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(257,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(257,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(258,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(258,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(258,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         
        n1 = 4 ; n2 = 6 ; 
        N(259,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(260,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(261,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(262,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(263,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(264,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(259,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(259,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(259,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(260,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(260,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(260,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(261,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(261,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(261,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(262,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(262,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(262,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(263,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(263,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(263,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(264,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(264,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(264,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         

        n1 = 6 ; n2 = 5 ; 
        N(265,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(266,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(267,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(268,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(269,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(270,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(265,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(265,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(265,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(266,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(266,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(266,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(267,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(267,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(267,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(268,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(268,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(268,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(269,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(269,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(269,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(270,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(270,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(270,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

        n1 = 5 ; n2 = 6 ; 
        N(271,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(272,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(273,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(274,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(275,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(276,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(271,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(271,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(271,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(272,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(272,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(272,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(273,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(273,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(273,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(274,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(274,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(274,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(275,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(275,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(275,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(276,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(276,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(276,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

        n1 = 6 ; n2 = 6 ; 
        N(277,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(278,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(279,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(280,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(281,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(282,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(277,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(277,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(277,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(278,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(278,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(278,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(279,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(279,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(279,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(280,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(280,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(280,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(281,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(281,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(281,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(282,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(282,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(282,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

          
% bubble function (needs to be changed)
% 2,3,4,5,6
% 226
        N(283,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(6,zeta) ;
        dN(283,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(6,zeta)  ;
        dN(283,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(6,zeta)  ;
        dN(283,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(6,zeta)  ;

% 336 
        N(284,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(6,zeta) ;
        dN(284,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(6,zeta)  ;
        dN(284,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(6,zeta)  ;
        dN(284,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(6,zeta)  ;
% 446
        N(285,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(6,zeta) ;
        dN(285,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(6,zeta)  ;
        dN(285,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(6,zeta)  ;
        dN(285,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(6,zeta)  ;
% 236     
        N(286,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(6,zeta) ;
        dN(286,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(6,zeta)  ;
        dN(286,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(6,zeta)  ;
        dN(286,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(6,zeta)  ;
% 326     
        N(287,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(6,zeta) ;
        dN(287,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(6,zeta)  ;
        dN(287,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(6,zeta)  ;
        dN(287,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(6,zeta)  ;
% 246
        N(288,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(6,zeta) ;
        dN(288,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(6,zeta)  ;
        dN(288,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(6,zeta)  ;
        dN(288,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(6,zeta)  ;
% 426   
        N(289,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(6,zeta) ;
        dN(289,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(6,zeta)  ;
        dN(289,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(6,zeta)  ;
        dN(289,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(6,zeta)  ;
% 346   
        N(290,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(6,zeta) ;
        dN(290,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(6,zeta)  ;
        dN(290,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(6,zeta)  ;
        dN(290,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(6,zeta)  ;
% 436   
        N(291,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(6,zeta) ;
        dN(291,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(6,zeta)  ;
        dN(291,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(6,zeta)  ;
        dN(291,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(6,zeta)  ;
% 622   
        N(292,1) = lobatto(6,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dN(292,1) =  dlobatto(6,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(292,2) =  lobatto(6,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(292,3) =  lobatto(6,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;
% 633   
        N(293,1) = lobatto(6,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dN(293,1) =  dlobatto(6,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(293,2) =  lobatto(6,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(293,3) =  lobatto(6,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;
% 644
        N(294,1) = lobatto(6,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(294,1) =  dlobatto(6,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(294,2) =  lobatto(6,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(294,3) =  lobatto(6,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
% 623     
        N(295,1) = lobatto(6,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dN(295,1) =  dlobatto(6,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(295,2) =  lobatto(6,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(295,3) =  lobatto(6,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
% 632     
        N(296,1) = lobatto(6,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        dN(296,1) =  dlobatto(6,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(296,2) =  lobatto(6,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(296,3) =  lobatto(6,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
% 624
        N(297,1) = lobatto(6,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(297,1) =  dlobatto(6,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(297,2) =  lobatto(6,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(297,3) =  lobatto(6,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;
% 642
        N(298,1) = lobatto(6,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dN(298,1) =  dlobatto(6,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(298,2) =  lobatto(6,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(298,3) =  lobatto(6,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
% 634
        N(299,1) = lobatto(6,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dN(299,1) =  dlobatto(6,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(299,2) =  lobatto(6,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(299,3) =  lobatto(6,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
% 643
        N(300,1) = lobatto(6,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(300,1) =  dlobatto(6,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(300,2) =  lobatto(6,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(300,3) =  lobatto(6,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
% 262
        N(301,1) = lobatto(2,xi)*lobatto(6,eta)*lobatto(2,zeta) ;
        dN(301,1) =  dlobatto(2,xi)*lobatto(6,eta)*lobatto(2,zeta)  ;
        dN(301,2) =  lobatto(2,xi)*dlobatto(6,eta)*lobatto(2,zeta)  ;
        dN(301,3) =  lobatto(2,xi)*lobatto(6,eta)*dlobatto(2,zeta)  ;
% 363
        N(302,1) = lobatto(3,xi)*lobatto(6,eta)*lobatto(3,zeta) ;
        dN(302,1) =  dlobatto(3,xi)*lobatto(6,eta)*lobatto(3,zeta)  ;
        dN(302,2) =  lobatto(3,xi)*dlobatto(6,eta)*lobatto(3,zeta)  ;
        dN(302,3) =  lobatto(3,xi)*lobatto(6,eta)*dlobatto(3,zeta)  ;
% 464
        N(303,1) = lobatto(4,xi)*lobatto(6,eta)*lobatto(4,zeta) ;
        dN(303,1) =  dlobatto(4,xi)*lobatto(6,eta)*lobatto(4,zeta)  ;
        dN(303,2) =  lobatto(4,xi)*dlobatto(6,eta)*lobatto(4,zeta)  ;
        dN(303,3) =  lobatto(4,xi)*lobatto(6,eta)*dlobatto(4,zeta)  ;
% 263
        N(304,1) = lobatto(2,xi)*lobatto(6,eta)*lobatto(3,zeta) ;
        dN(304,1) =  dlobatto(2,xi)*lobatto(6,eta)*lobatto(3,zeta)  ;
        dN(304,2) =  lobatto(2,xi)*dlobatto(6,eta)*lobatto(3,zeta)  ;
        dN(304,3) =  lobatto(2,xi)*lobatto(6,eta)*dlobatto(3,zeta)  ;
% 362
        N(305,1) = lobatto(3,xi)*lobatto(6,eta)*lobatto(2,zeta) ;
        dN(305,1) =  dlobatto(3,xi)*lobatto(6,eta)*lobatto(2,zeta)  ;
        dN(305,2) =  lobatto(3,xi)*dlobatto(6,eta)*lobatto(2,zeta)  ;
        dN(305,3) =  lobatto(3,xi)*lobatto(6,eta)*dlobatto(2,zeta)  ;
% 264
        N(306,1) = lobatto(2,xi)*lobatto(6,eta)*lobatto(4,zeta) ;
        dN(306,1) =  dlobatto(2,xi)*lobatto(6,eta)*lobatto(4,zeta)  ;
        dN(306,2) =  lobatto(2,xi)*dlobatto(6,eta)*lobatto(4,zeta)  ;
        dN(306,3) =  lobatto(2,xi)*lobatto(6,eta)*dlobatto(4,zeta)  ;
% 462
        N(307,1) = lobatto(4,xi)*lobatto(6,eta)*lobatto(2,zeta) ;
        dN(307,1) =  dlobatto(4,xi)*lobatto(6,eta)*lobatto(2,zeta)  ;
        dN(307,2) =  lobatto(4,xi)*dlobatto(6,eta)*lobatto(2,zeta)  ;
        dN(307,3) =  lobatto(4,xi)*lobatto(6,eta)*dlobatto(2,zeta)  ;
% 364
        N(308,1) = lobatto(3,xi)*lobatto(6,eta)*lobatto(4,zeta) ;
        dN(308,1) =  dlobatto(3,xi)*lobatto(6,eta)*lobatto(4,zeta)  ;
        dN(308,2) =  lobatto(3,xi)*dlobatto(6,eta)*lobatto(4,zeta)  ;
        dN(308,3) =  lobatto(3,xi)*lobatto(6,eta)*dlobatto(4,zeta)  ;
% 463
        N(309,1) = lobatto(4,xi)*lobatto(6,eta)*lobatto(3,zeta) ;
        dN(309,1) =  dlobatto(4,xi)*lobatto(6,eta)*lobatto(3,zeta)  ;
        dN(309,2) =  lobatto(4,xi)*dlobatto(6,eta)*lobatto(3,zeta)  ;
        dN(309,3) =  lobatto(4,xi)*lobatto(6,eta)*dlobatto(3,zeta)  ;
% 266
        N(310,1) = lobatto(2,xi)*lobatto(6,eta)*lobatto(6,zeta) ;
        dN(310,1) =  dlobatto(2,xi)*lobatto(6,eta)*lobatto(6,zeta)  ;
        dN(310,2) =  lobatto(2,xi)*dlobatto(6,eta)*lobatto(6,zeta)  ;
        dN(310,3) =  lobatto(2,xi)*lobatto(6,eta)*dlobatto(6,zeta)  ;
% 366
        N(311,1) = lobatto(3,xi)*lobatto(6,eta)*lobatto(6,zeta) ;
        dN(311,1) =  dlobatto(3,xi)*lobatto(6,eta)*lobatto(6,zeta)  ;
        dN(311,2) =  lobatto(3,xi)*dlobatto(6,eta)*lobatto(6,zeta)  ;
        dN(311,3) =  lobatto(3,xi)*lobatto(6,eta)*dlobatto(6,zeta)  ;
% 466
        N(312,1) = lobatto(4,xi)*lobatto(6,eta)*lobatto(6,zeta) ;
        dN(312,1) =  dlobatto(4,xi)*lobatto(6,eta)*lobatto(6,zeta)  ;
        dN(312,2) =  lobatto(4,xi)*dlobatto(6,eta)*lobatto(6,zeta)  ;
        dN(312,3) =  lobatto(4,xi)*lobatto(6,eta)*dlobatto(6,zeta)  ;
% 626
        N(313,1) = lobatto(6,xi)*lobatto(2,eta)*lobatto(6,zeta) ;
        dN(313,1) =  dlobatto(6,xi)*lobatto(2,eta)*lobatto(6,zeta)  ;
        dN(313,2) =  lobatto(6,xi)*dlobatto(2,eta)*lobatto(6,zeta)  ;
        dN(313,3) =  lobatto(6,xi)*lobatto(2,eta)*dlobatto(6,zeta)  ;
% 636
        N(314,1) = lobatto(6,xi)*lobatto(3,eta)*lobatto(6,zeta) ;
        dN(314,1) =  dlobatto(6,xi)*lobatto(3,eta)*lobatto(6,zeta)  ;
        dN(314,2) =  lobatto(6,xi)*dlobatto(3,eta)*lobatto(6,zeta)  ;
        dN(314,3) =  lobatto(6,xi)*lobatto(3,eta)*dlobatto(6,zeta)  ;
% 646
        N(315,1) = lobatto(6,xi)*lobatto(4,eta)*lobatto(6,zeta) ;
        dN(315,1) =  dlobatto(6,xi)*lobatto(4,eta)*lobatto(6,zeta)  ;
        dN(315,2) =  lobatto(6,xi)*dlobatto(4,eta)*lobatto(6,zeta)  ;
        dN(315,3) =  lobatto(6,xi)*lobatto(4,eta)*dlobatto(6,zeta)  ;
% 662
        N(316,1) = lobatto(6,xi)*lobatto(6,eta)*lobatto(2,zeta) ;
        dN(316,1) =  dlobatto(6,xi)*lobatto(6,eta)*lobatto(2,zeta)  ;
        dN(316,2) =  lobatto(6,xi)*dlobatto(6,eta)*lobatto(2,zeta)  ;
        dN(316,3) =  lobatto(6,xi)*lobatto(6,eta)*dlobatto(2,zeta)  ;
% 663
        N(317,1) = lobatto(6,xi)*lobatto(6,eta)*lobatto(3,zeta) ;
        dN(317,1) =  dlobatto(6,xi)*lobatto(6,eta)*lobatto(3,zeta)  ;
        dN(317,2) =  lobatto(6,xi)*dlobatto(6,eta)*lobatto(3,zeta)  ;
        dN(317,3) =  lobatto(6,xi)*lobatto(6,eta)*dlobatto(3,zeta)  ;
% 664
        N(318,1) = lobatto(6,xi)*lobatto(6,eta)*lobatto(4,zeta) ;
        dN(318,1) =  dlobatto(6,xi)*lobatto(6,eta)*lobatto(4,zeta)  ;
        dN(318,2) =  lobatto(6,xi)*dlobatto(6,eta)*lobatto(4,zeta)  ;
        dN(318,3) =  lobatto(6,xi)*lobatto(6,eta)*dlobatto(4,zeta)  ;

% 556 
        N(319,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(6,zeta) ;
        dN(319,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(6,zeta)  ;
        dN(319,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(6,zeta)  ;
        dN(319,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(6,zeta)  ;

% 256
        N(320,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(6,zeta) ;
        dN(320,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(6,zeta)  ;
        dN(320,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(6,zeta)  ;
        dN(320,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(6,zeta)  ;
% 356
        N(321,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(6,zeta) ;
        dN(321,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(6,zeta)  ;
        dN(321,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(6,zeta)  ;
        dN(321,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(6,zeta)  ;
% 536
        N(322,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(6,zeta) ;
        dN(322,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(6,zeta)  ;
        dN(322,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(6,zeta)  ;
        dN(322,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(6,zeta)  ;
% 526
        N(323,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(6,zeta) ;
        dN(323,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(6,zeta)  ;
        dN(323,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(6,zeta)  ;
        dN(323,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(6,zeta)  ;
% 546
        N(324,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(6,zeta) ;
        dN(324,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(6,zeta)  ;
        dN(324,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(6,zeta)  ;
        dN(324,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(6,zeta)  ;
% 456
        N(325,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(6,zeta) ;
        dN(325,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(6,zeta)  ;
        dN(325,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(6,zeta)  ;
        dN(325,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(6,zeta)  ;
% 655
        N(326,1) = lobatto(6,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(326,1) =  dlobatto(6,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(326,2) =  lobatto(6,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(326,3) =  lobatto(6,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 625
        N(327,1) = lobatto(6,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(327,1) =  dlobatto(6,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(327,2) =  lobatto(6,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(327,3) =  lobatto(6,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 635
        N(328,1) = lobatto(6,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(328,1) =  dlobatto(6,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(328,2) =  lobatto(6,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(328,3) =  lobatto(6,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 652
        N(329,1) = lobatto(6,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(329,1) =  dlobatto(6,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(329,2) =  lobatto(6,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(329,3) =  lobatto(6,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 653
        N(330,1) = lobatto(6,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(330,1) =  dlobatto(6,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(330,2) =  lobatto(6,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(330,3) =  lobatto(6,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 654
        N(331,1) = lobatto(6,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(331,1) =  dlobatto(6,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(331,2) =  lobatto(6,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(331,3) =  lobatto(6,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 645
        N(332,1) = lobatto(6,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(332,1) =  dlobatto(6,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(332,2) =  lobatto(6,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(332,3) =  lobatto(6,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 565
        N(333,1) = lobatto(5,xi)*lobatto(6,eta)*lobatto(5,zeta) ;
        dN(333,1) =  dlobatto(5,xi)*lobatto(6,eta)*lobatto(5,zeta)  ;
        dN(333,2) =  lobatto(5,xi)*dlobatto(6,eta)*lobatto(5,zeta)  ;
        dN(333,3) =  lobatto(5,xi)*lobatto(6,eta)*dlobatto(5,zeta)  ;
% 265
        N(334,1) = lobatto(2,xi)*lobatto(6,eta)*lobatto(5,zeta) ;
        dN(334,1) =  dlobatto(2,xi)*lobatto(6,eta)*lobatto(5,zeta)  ;
        dN(334,2) =  lobatto(2,xi)*dlobatto(6,eta)*lobatto(5,zeta)  ;
        dN(334,3) =  lobatto(2,xi)*lobatto(6,eta)*dlobatto(5,zeta)  ;
% 562
        N(335,1) = lobatto(5,xi)*lobatto(6,eta)*lobatto(2,zeta) ;
        dN(335,1) =  dlobatto(5,xi)*lobatto(6,eta)*lobatto(2,zeta)  ;
        dN(335,2) =  lobatto(5,xi)*dlobatto(6,eta)*lobatto(2,zeta)  ;
        dN(335,3) =  lobatto(5,xi)*lobatto(6,eta)*dlobatto(2,zeta)  ;
% 365
        N(336,1) = lobatto(3,xi)*lobatto(6,eta)*lobatto(5,zeta) ;
        dN(336,1) =  dlobatto(3,xi)*lobatto(6,eta)*lobatto(5,zeta)  ;
        dN(336,2) =  lobatto(3,xi)*dlobatto(6,eta)*lobatto(5,zeta)  ;
        dN(336,3) =  lobatto(3,xi)*lobatto(6,eta)*dlobatto(5,zeta)  ;
% 563
        N(337,1) = lobatto(5,xi)*lobatto(6,eta)*lobatto(3,zeta) ;
        dN(337,1) =  dlobatto(5,xi)*lobatto(6,eta)*lobatto(3,zeta)  ;
        dN(337,2) =  lobatto(5,xi)*dlobatto(6,eta)*lobatto(3,zeta)  ;
        dN(337,3) =  lobatto(5,xi)*lobatto(6,eta)*dlobatto(3,zeta)  ;
% 465
        N(338,1) = lobatto(4,xi)*lobatto(6,eta)*lobatto(5,zeta) ;
        dN(338,1) =  dlobatto(4,xi)*lobatto(6,eta)*lobatto(5,zeta)  ;
        dN(338,2) =  lobatto(4,xi)*dlobatto(6,eta)*lobatto(5,zeta)  ;
        dN(338,3) =  lobatto(4,xi)*lobatto(6,eta)*dlobatto(5,zeta)  ;
% 564
        N(339,1) = lobatto(5,xi)*lobatto(6,eta)*lobatto(4,zeta) ;
        dN(339,1) =  dlobatto(5,xi)*lobatto(6,eta)*lobatto(4,zeta)  ;
        dN(339,2) =  lobatto(5,xi)*dlobatto(6,eta)*lobatto(4,zeta)  ;
        dN(339,3) =  lobatto(5,xi)*lobatto(6,eta)*dlobatto(4,zeta)  ;
% 665
        N(340,1) = lobatto(6,xi)*lobatto(6,eta)*lobatto(5,zeta) ;
        dN(340,1) =  dlobatto(6,xi)*lobatto(6,eta)*lobatto(5,zeta)  ;
        dN(340,2) =  lobatto(6,xi)*dlobatto(6,eta)*lobatto(5,zeta)  ;
        dN(340,3) =  lobatto(6,xi)*lobatto(6,eta)*dlobatto(5,zeta)  ;
% 656
        N(341,1) = lobatto(6,xi)*lobatto(5,eta)*lobatto(6,zeta) ;
        dN(341,1) =  dlobatto(6,xi)*lobatto(5,eta)*lobatto(6,zeta)  ;
        dN(341,2) =  lobatto(6,xi)*dlobatto(5,eta)*lobatto(6,zeta)  ;
        dN(341,3) =  lobatto(6,xi)*lobatto(5,eta)*dlobatto(6,zeta)  ;
% 566
        N(342,1) = lobatto(5,xi)*lobatto(6,eta)*lobatto(6,zeta) ;
        dN(342,1) =  dlobatto(5,xi)*lobatto(6,eta)*lobatto(6,zeta)  ;
        dN(342,2) =  lobatto(5,xi)*dlobatto(6,eta)*lobatto(6,zeta)  ;
        dN(342,3) =  lobatto(5,xi)*lobatto(6,eta)*dlobatto(6,zeta)  ;

% 666
        N(343,1) = lobatto(6,xi)*lobatto(6,eta)*lobatto(6,zeta) ;
        dN(343,1) =  dlobatto(6,xi)*lobatto(6,eta)*lobatto(6,zeta)  ;
        dN(343,2) =  lobatto(6,xi)*dlobatto(6,eta)*lobatto(6,zeta)  ;
        dN(343,3) =  lobatto(6,xi)*lobatto(6,eta)*dlobatto(6,zeta)  ;

        end
          
        if order >= 7
        k = 7 ; 
        % edge functions
        N(344,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        N(345,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(346,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        N(347,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        N(348,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(349,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        N(350,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(351,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        N(352,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        N(353,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        N(354,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        N(355,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dN(344,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dN(344,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dN(344,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dN(345,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(345,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(345,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dN(346,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dN(346,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dN(346,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dN(347,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dN(347,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dN(347,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dN(348,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(348,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(348,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(349,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dN(349,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dN(349,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dN(350,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(350,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(350,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(351,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dN(351,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dN(351,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dN(352,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dN(352,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dN(352,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dN(353,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(353,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(353,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dN(354,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dN(354,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dN(354,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dN(355,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dN(355,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dN(355,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 7 ; 
        N(356,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(357,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(358,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(359,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(360,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(361,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(356,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(356,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(356,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(357,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(357,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(357,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(358,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(358,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(358,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(359,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(359,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(359,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(360,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(360,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(360,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(361,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(361,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(361,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 7 ; n2 = 2 ; 
        N(362,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(363,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(364,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(365,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(366,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(367,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(362,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(362,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(362,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(363,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(363,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(363,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(364,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(364,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(364,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(365,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(365,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(365,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(366,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(366,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(366,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(367,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(367,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(367,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 7 ; n2 = 3 ; 
        N(368,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(369,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(370,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(371,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(372,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(373,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(368,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(368,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(368,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(369,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(369,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(369,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(370,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(370,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(370,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(371,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(371,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(371,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(372,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(372,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(372,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(373,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(373,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(373,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 7 ; 
        N(374,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(375,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(376,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(377,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(378,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(379,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(374,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(374,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(374,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(375,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(375,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(375,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(376,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(376,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(376,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(377,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(377,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(377,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(378,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(378,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(378,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(379,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(379,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(379,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
       
        n1 = 7 ; n2 = 4 ; 
        N(380,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(381,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(382,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(383,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(384,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(385,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(380,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(380,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(380,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(381,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(381,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(381,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(382,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(382,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(382,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(383,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(383,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(383,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(384,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(384,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(384,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(385,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(385,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(385,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         
        n1 = 4 ; n2 = 7 ; 
        N(386,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(387,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(388,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(389,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(390,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(391,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(386,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(386,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(386,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(387,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(387,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(387,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(388,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(388,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(388,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(389,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(389,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(389,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(390,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(390,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(390,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(391,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(391,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(391,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         

        n1 = 7 ; n2 = 5 ; 
        N(392,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(393,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(394,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(395,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(396,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(397,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(392,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(392,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(392,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(393,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(393,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(393,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(394,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(394,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(394,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(395,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(395,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(395,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(396,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(396,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(396,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(397,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(397,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(397,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

        n1 = 5 ; n2 = 7 ; 
        N(398,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(399,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(400,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(401,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(402,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(403,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(398,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(398,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(398,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(399,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(399,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(399,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(400,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(400,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(400,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(401,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(401,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(401,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(402,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(402,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(402,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(403,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(403,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(403,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

        n1 = 6 ; n2 = 7 ; 
        N(404,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(405,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(406,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(407,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(408,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(409,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(404,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(404,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(404,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(405,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(405,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(405,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(406,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(406,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(406,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(407,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(407,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(407,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(408,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(408,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(408,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(409,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(409,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(409,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;


        n1 = 7 ; n2 = 6 ; 
        N(410,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(411,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(412,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(413,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(414,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(415,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(410,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(410,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(410,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(411,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(411,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(411,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(412,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(412,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(412,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(413,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(413,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(413,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(414,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(414,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(414,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(415,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(415,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(415,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          


        n1 = 7 ; n2 = 7 ; 
        N(416,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(417,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        N(418,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        N(419,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        N(420,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        N(421,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dN(416,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(416,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(416,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(417,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(417,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dN(417,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dN(418,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(418,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dN(418,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dN(419,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(419,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dN(419,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dN(420,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(420,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dN(420,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dN(421,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(421,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dN(421,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          

% bubble function (needs to be changed)
% 2,3,4,5,6,7
% 227
        N(422,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(7,zeta) ;
        dN(422,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(7,zeta)  ;
        dN(422,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(7,zeta)  ;
        dN(422,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(7,zeta)  ;

% 337 
        N(423,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(7,zeta) ;
        dN(423,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(7,zeta)  ;
        dN(423,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(7,zeta)  ;
        dN(423,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(7,zeta)  ;
% 447
        N(424,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(7,zeta) ;
        dN(424,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(7,zeta)  ;
        dN(424,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(7,zeta)  ;
        dN(424,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(7,zeta)  ;
% 237     
        N(425,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(7,zeta) ;
        dN(425,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(7,zeta)  ;
        dN(425,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(7,zeta)  ;
        dN(425,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(7,zeta)  ;
% 327     
        N(426,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(7,zeta) ;
        dN(426,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(7,zeta)  ;
        dN(426,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(7,zeta)  ;
        dN(426,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(7,zeta)  ;
% 247
        N(427,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(7,zeta) ;
        dN(427,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(7,zeta)  ;
        dN(427,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(7,zeta)  ;
        dN(427,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(7,zeta)  ;
% 427   
        N(428,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(7,zeta) ;
        dN(428,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(7,zeta)  ;
        dN(428,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(7,zeta)  ;
        dN(428,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(7,zeta)  ;
% 347   
        N(429,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(7,zeta) ;
        dN(429,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(7,zeta)  ;
        dN(429,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(7,zeta)  ;
        dN(429,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(7,zeta)  ;
% 437   
        N(430,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(7,zeta) ;
        dN(430,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(7,zeta)  ;
        dN(430,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(7,zeta)  ;
        dN(430,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(7,zeta)  ;
% 722   
        N(431,1) = lobatto(7,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dN(431,1) =  dlobatto(7,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dN(431,2) =  lobatto(7,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dN(431,3) =  lobatto(7,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;
% 733   
        N(432,1) = lobatto(7,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dN(432,1) =  dlobatto(7,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dN(432,2) =  lobatto(7,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dN(432,3) =  lobatto(7,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;
% 744
        N(433,1) = lobatto(7,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dN(433,1) =  dlobatto(7,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dN(433,2) =  lobatto(7,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dN(433,3) =  lobatto(7,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
% 723     
        N(434,1) = lobatto(7,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dN(434,1) =  dlobatto(7,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dN(434,2) =  lobatto(7,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dN(434,3) =  lobatto(7,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
% 732     
        N(435,1) = lobatto(7,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        dN(435,1) =  dlobatto(7,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dN(435,2) =  lobatto(7,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dN(435,3) =  lobatto(7,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
% 724
        N(436,1) = lobatto(7,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dN(436,1) =  dlobatto(7,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dN(436,2) =  lobatto(7,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dN(436,3) =  lobatto(7,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;
% 742
        N(437,1) = lobatto(7,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dN(437,1) =  dlobatto(7,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dN(437,2) =  lobatto(7,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dN(437,3) =  lobatto(7,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
% 734
        N(438,1) = lobatto(7,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dN(438,1) =  dlobatto(7,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dN(438,2) =  lobatto(7,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dN(438,3) =  lobatto(7,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
% 743
        N(439,1) = lobatto(7,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dN(439,1) =  dlobatto(7,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dN(439,2) =  lobatto(7,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dN(439,3) =  lobatto(7,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
% 272
        N(440,1) = lobatto(2,xi)*lobatto(7,eta)*lobatto(2,zeta) ;
        dN(440,1) =  dlobatto(2,xi)*lobatto(7,eta)*lobatto(2,zeta)  ;
        dN(440,2) =  lobatto(2,xi)*dlobatto(7,eta)*lobatto(2,zeta)  ;
        dN(440,3) =  lobatto(2,xi)*lobatto(7,eta)*dlobatto(2,zeta)  ;
% 373
        N(441,1) = lobatto(3,xi)*lobatto(7,eta)*lobatto(3,zeta) ;
        dN(441,1) =  dlobatto(3,xi)*lobatto(7,eta)*lobatto(3,zeta)  ;
        dN(441,2) =  lobatto(3,xi)*dlobatto(7,eta)*lobatto(3,zeta)  ;
        dN(441,3) =  lobatto(3,xi)*lobatto(7,eta)*dlobatto(3,zeta)  ;
% 474
        N(442,1) = lobatto(4,xi)*lobatto(7,eta)*lobatto(4,zeta) ;
        dN(442,1) =  dlobatto(4,xi)*lobatto(7,eta)*lobatto(4,zeta)  ;
        dN(442,2) =  lobatto(4,xi)*dlobatto(7,eta)*lobatto(4,zeta)  ;
        dN(442,3) =  lobatto(4,xi)*lobatto(7,eta)*dlobatto(4,zeta)  ;
% 273
        N(443,1) = lobatto(2,xi)*lobatto(7,eta)*lobatto(3,zeta) ;
        dN(443,1) =  dlobatto(2,xi)*lobatto(7,eta)*lobatto(3,zeta)  ;
        dN(443,2) =  lobatto(2,xi)*dlobatto(7,eta)*lobatto(3,zeta)  ;
        dN(443,3) =  lobatto(2,xi)*lobatto(7,eta)*dlobatto(3,zeta)  ;
% 372
        N(444,1) = lobatto(3,xi)*lobatto(7,eta)*lobatto(2,zeta) ;
        dN(444,1) =  dlobatto(3,xi)*lobatto(7,eta)*lobatto(2,zeta)  ;
        dN(444,2) =  lobatto(3,xi)*dlobatto(7,eta)*lobatto(2,zeta)  ;
        dN(444,3) =  lobatto(3,xi)*lobatto(7,eta)*dlobatto(2,zeta)  ;
% 274
        N(445,1) = lobatto(2,xi)*lobatto(7,eta)*lobatto(4,zeta) ;
        dN(445,1) =  dlobatto(2,xi)*lobatto(7,eta)*lobatto(4,zeta)  ;
        dN(445,2) =  lobatto(2,xi)*dlobatto(7,eta)*lobatto(4,zeta)  ;
        dN(445,3) =  lobatto(2,xi)*lobatto(7,eta)*dlobatto(4,zeta)  ;
% 374
        N(446,1) = lobatto(3,xi)*lobatto(7,eta)*lobatto(4,zeta) ;
        dN(446,1) =  dlobatto(3,xi)*lobatto(7,eta)*lobatto(4,zeta)  ;
        dN(446,2) =  lobatto(3,xi)*dlobatto(7,eta)*lobatto(4,zeta)  ;
        dN(446,3) =  lobatto(3,xi)*lobatto(7,eta)*dlobatto(4,zeta)  ;
% 473
        N(447,1) = lobatto(4,xi)*lobatto(7,eta)*lobatto(3,zeta) ;
        dN(447,1) =  dlobatto(4,xi)*lobatto(7,eta)*lobatto(3,zeta)  ;
        dN(447,2) =  lobatto(4,xi)*dlobatto(7,eta)*lobatto(3,zeta)  ;
        dN(447,3) =  lobatto(4,xi)*lobatto(7,eta)*dlobatto(3,zeta)  ;
% 277
        N(448,1) = lobatto(2,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(448,1) =  dlobatto(2,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(448,2) =  lobatto(2,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(448,3) =  lobatto(2,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;
% 377
        N(449,1) = lobatto(3,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(449,1) =  dlobatto(3,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(449,2) =  lobatto(3,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(449,3) =  lobatto(3,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;
% 477
        N(450,1) = lobatto(4,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(450,1) =  dlobatto(4,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(450,2) =  lobatto(4,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(450,3) =  lobatto(4,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;
% 727
        N(451,1) = lobatto(7,xi)*lobatto(2,eta)*lobatto(7,zeta) ;
        dN(451,1) =  dlobatto(7,xi)*lobatto(2,eta)*lobatto(7,zeta)  ;
        dN(451,2) =  lobatto(7,xi)*dlobatto(2,eta)*lobatto(7,zeta)  ;
        dN(451,3) =  lobatto(7,xi)*lobatto(2,eta)*dlobatto(7,zeta)  ;
% 737
        N(452,1) = lobatto(7,xi)*lobatto(3,eta)*lobatto(7,zeta) ;
        dN(452,1) =  dlobatto(7,xi)*lobatto(3,eta)*lobatto(7,zeta)  ;
        dN(452,2) =  lobatto(7,xi)*dlobatto(3,eta)*lobatto(7,zeta)  ;
        dN(452,3) =  lobatto(7,xi)*lobatto(3,eta)*dlobatto(7,zeta)  ;
% 747
        N(453,1) = lobatto(7,xi)*lobatto(4,eta)*lobatto(7,zeta) ;
        dN(453,1) =  dlobatto(7,xi)*lobatto(4,eta)*lobatto(7,zeta)  ;
        dN(453,2) =  lobatto(7,xi)*dlobatto(4,eta)*lobatto(7,zeta)  ;
        dN(453,3) =  lobatto(7,xi)*lobatto(4,eta)*dlobatto(7,zeta)  ;
% 772
        N(454,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(2,zeta) ;
        dN(454,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(2,zeta)  ;
        dN(454,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(2,zeta)  ;
        dN(454,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(2,zeta)  ;
% 773
        N(455,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(3,zeta) ;
        dN(455,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(3,zeta)  ;
        dN(455,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(3,zeta)  ;
        dN(455,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(3,zeta)  ;
% 774
        N(456,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(4,zeta) ;
        dN(456,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(4,zeta)  ;
        dN(456,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(4,zeta)  ;
        dN(456,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(4,zeta)  ;

% 557 
        N(457,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(7,zeta) ;
        dN(457,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(7,zeta)  ;
        dN(457,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(7,zeta)  ;
        dN(457,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(7,zeta)  ;

% 257
        N(458,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(7,zeta) ;
        dN(458,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(7,zeta)  ;
        dN(458,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(7,zeta)  ;
        dN(458,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(7,zeta)  ;
% 357
        N(459,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(7,zeta) ;
        dN(459,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(7,zeta)  ;
        dN(459,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(7,zeta)  ;
        dN(459,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(7,zeta)  ;
% 537
        N(460,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(7,zeta) ;
        dN(460,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(7,zeta)  ;
        dN(460,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(7,zeta)  ;
        dN(460,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(7,zeta)  ;
% 527
        N(461,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(7,zeta) ;
        dN(461,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(7,zeta)  ;
        dN(461,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(7,zeta)  ;
        dN(461,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(7,zeta)  ;
% 547
        N(462,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(7,zeta) ;
        dN(462,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(7,zeta)  ;
        dN(462,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(7,zeta)  ;
        dN(462,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(7,zeta)  ;
% 457
        N(463,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(7,zeta) ;
        dN(463,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(7,zeta)  ;
        dN(463,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(7,zeta)  ;
        dN(463,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(7,zeta)  ;
% 755
        N(464,1) = lobatto(7,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dN(464,1) =  dlobatto(7,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dN(464,2) =  lobatto(7,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dN(464,3) =  lobatto(7,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 725
        N(465,1) = lobatto(7,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dN(465,1) =  dlobatto(7,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dN(465,2) =  lobatto(7,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dN(465,3) =  lobatto(7,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 735
        N(466,1) = lobatto(7,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dN(466,1) =  dlobatto(7,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dN(466,2) =  lobatto(7,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dN(466,3) =  lobatto(7,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 752
        N(467,1) = lobatto(7,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dN(467,1) =  dlobatto(7,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dN(467,2) =  lobatto(7,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dN(467,3) =  lobatto(7,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 753
        N(468,1) = lobatto(7,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dN(468,1) =  dlobatto(7,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dN(468,2) =  lobatto(7,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dN(468,3) =  lobatto(7,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 754
        N(469,1) = lobatto(7,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dN(469,1) =  dlobatto(7,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dN(469,2) =  lobatto(7,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dN(469,3) =  lobatto(7,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 745
        N(470,1) = lobatto(7,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dN(470,1) =  dlobatto(7,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dN(470,2) =  lobatto(7,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dN(470,3) =  lobatto(7,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 575
        N(471,1) = lobatto(5,xi)*lobatto(7,eta)*lobatto(5,zeta) ;
        dN(471,1) =  dlobatto(5,xi)*lobatto(7,eta)*lobatto(5,zeta)  ;
        dN(471,2) =  lobatto(5,xi)*dlobatto(7,eta)*lobatto(5,zeta)  ;
        dN(471,3) =  lobatto(5,xi)*lobatto(7,eta)*dlobatto(5,zeta)  ;
% 275
        N(472,1) = lobatto(2,xi)*lobatto(7,eta)*lobatto(5,zeta) ;
        dN(472,1) =  dlobatto(2,xi)*lobatto(7,eta)*lobatto(5,zeta)  ;
        dN(472,2) =  lobatto(2,xi)*dlobatto(7,eta)*lobatto(5,zeta)  ;
        dN(472,3) =  lobatto(2,xi)*lobatto(7,eta)*dlobatto(5,zeta)  ;
% 572
        N(473,1) = lobatto(5,xi)*lobatto(7,eta)*lobatto(2,zeta) ;
        dN(473,1) =  dlobatto(5,xi)*lobatto(7,eta)*lobatto(2,zeta)  ;
        dN(473,2) =  lobatto(5,xi)*dlobatto(7,eta)*lobatto(2,zeta)  ;
        dN(473,3) =  lobatto(5,xi)*lobatto(7,eta)*dlobatto(2,zeta)  ;
% 375
        N(474,1) = lobatto(3,xi)*lobatto(7,eta)*lobatto(5,zeta) ;
        dN(474,1) =  dlobatto(3,xi)*lobatto(7,eta)*lobatto(5,zeta)  ;
        dN(474,2) =  lobatto(3,xi)*dlobatto(7,eta)*lobatto(5,zeta)  ;
        dN(474,3) =  lobatto(3,xi)*lobatto(7,eta)*dlobatto(5,zeta)  ;
% 573
        N(475,1) = lobatto(5,xi)*lobatto(7,eta)*lobatto(3,zeta) ;
        dN(475,1) =  dlobatto(5,xi)*lobatto(7,eta)*lobatto(3,zeta)  ;
        dN(475,2) =  lobatto(5,xi)*dlobatto(7,eta)*lobatto(3,zeta)  ;
        dN(475,3) =  lobatto(5,xi)*lobatto(7,eta)*dlobatto(3,zeta)  ;
% 475
        N(476,1) = lobatto(4,xi)*lobatto(7,eta)*lobatto(5,zeta) ;
        dN(476,1) =  dlobatto(4,xi)*lobatto(7,eta)*lobatto(5,zeta)  ;
        dN(476,2) =  lobatto(4,xi)*dlobatto(7,eta)*lobatto(5,zeta)  ;
        dN(476,3) =  lobatto(4,xi)*lobatto(7,eta)*dlobatto(5,zeta)  ;
% 574
        N(477,1) = lobatto(5,xi)*lobatto(7,eta)*lobatto(4,zeta) ;
        dN(477,1) =  dlobatto(5,xi)*lobatto(7,eta)*lobatto(4,zeta)  ;
        dN(477,2) =  lobatto(5,xi)*dlobatto(7,eta)*lobatto(4,zeta)  ;
        dN(477,3) =  lobatto(5,xi)*lobatto(7,eta)*dlobatto(4,zeta)  ;
% 775
        N(478,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(5,zeta) ;
        dN(478,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(5,zeta)  ;
        dN(478,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(5,zeta)  ;
        dN(478,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(5,zeta)  ;
% 757
        N(479,1) = lobatto(7,xi)*lobatto(5,eta)*lobatto(7,zeta) ;
        dN(479,1) =  dlobatto(7,xi)*lobatto(5,eta)*lobatto(7,zeta)  ;
        dN(479,2) =  lobatto(7,xi)*dlobatto(5,eta)*lobatto(7,zeta)  ;
        dN(479,3) =  lobatto(7,xi)*lobatto(5,eta)*dlobatto(7,zeta)  ;
% 577
        N(480,1) = lobatto(5,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(480,1) =  dlobatto(5,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(480,2) =  lobatto(5,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(480,3) =  lobatto(5,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;
% 472
        N(481,1) = lobatto(4,xi)*lobatto(7,eta)*lobatto(2,zeta) ;
        dN(481,1) =  dlobatto(4,xi)*lobatto(7,eta)*lobatto(2,zeta)  ;
        dN(481,2) =  lobatto(4,xi)*dlobatto(7,eta)*lobatto(2,zeta)  ;
        dN(481,3) =  lobatto(4,xi)*lobatto(7,eta)*dlobatto(2,zeta)  ;

% 777
        N(482,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(482,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(482,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(482,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;
% 2,3,4,5,6,7
       
% 627
        N(483,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(483,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(483,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(483,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;


% 672
        N(484,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(484,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(484,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(484,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 637
        N(485,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(485,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(485,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(485,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 673
        N(486,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(486,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(486,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(486,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 647
        N(487,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(487,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(487,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(487,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 674
        N(488,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(488,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(488,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(488,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 657
        N(489,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(489,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(489,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(489,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 675
        N(490,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(490,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(490,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(490,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 267
        N(491,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(491,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(491,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(491,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 762
        N(492,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(492,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(492,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(492,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 367
        N(493,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(493,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(493,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(493,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 763
        N(494,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(494,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(494,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(494,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 467
        N(495,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(495,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(495,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(495,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 764
        N(496,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(496,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(496,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(496,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 567
        N(497,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(497,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(497,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(497,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 765
        N(498,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(498,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(498,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(498,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 276
        N(499,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(499,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(499,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(499,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 726
        N(500,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(500,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(500,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(500,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 376
        N(501,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(501,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(501,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(501,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 736
        N(502,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(502,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(502,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(502,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 476
        N(503,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(503,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(503,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(503,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 746
        N(504,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(504,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(504,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(504,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 576
        N(505,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(505,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(505,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(505,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 756
        N(506,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(506,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(506,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(506,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 677
        N(507,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(507,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(507,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(507,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 766
        N(508,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(508,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(508,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(508,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 676
        N(509,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(509,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(509,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(509,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 767
        N(510,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(510,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(510,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(510,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 776
        N(511,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(511,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(511,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(511,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

% 762
        N(512,1) = lobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta) ;
        dN(512,1) =  dlobatto(7,xi)*lobatto(7,eta)*lobatto(7,zeta)  ;
        dN(512,2) =  lobatto(7,xi)*dlobatto(7,eta)*lobatto(7,zeta)  ;
        dN(512,3) =  lobatto(7,xi)*lobatto(7,eta)*dlobatto(7,zeta)  ;

        end

       %{
        if order >= 7
        k = 6 ; 
% edge functions
        [N(344),dN(344,1),dN(344,2),dN(344,3)] = eval_lobatto_edge (k,0,0,xi,eta,zeta) ;
        [N(345),dN(345,1),dN(345,2),dN(345,3)] = eval_lobatto_edge (1,k,0,xi,eta,zeta) ;
        [N(346),dN(346,1),dN(346,2),dN(346,3)] = eval_lobatto_edge (k,1,0,xi,eta,zeta) ;
        [N(347),dN(347,1),dN(347,2),dN(347,3)] = eval_lobatto_edge (0,k,0,xi,eta,zeta) ;
        [N(348),dN(348,1),dN(348,2),dN(348,3)] = eval_lobatto_edge (0,0,k,xi,eta,zeta) ;
        [N(349),dN(349,1),dN(349,2),dN(349,3)] = eval_lobatto_edge (1,0,k,xi,eta,zeta) ;
        [N(350),dN(350,1),dN(350,2),dN(350,3)] = eval_lobatto_edge (1,1,k,xi,eta,zeta) ;
        [N(351),dN(351,1),dN(351,2),dN(351,3)] = eval_lobatto_edge (0,1,k,xi,eta,zeta) ;
        [N(352),dN(352,1),dN(352,2),dN(352,3)] = eval_lobatto_edge (k,0,1,xi,eta,zeta) ;
        [N(353),dN(353,1),dN(353,2),dN(353,3)] = eval_lobatto_edge (1,k,1,xi,eta,zeta) ;
        [N(354),dN(354,1),dN(354,2),dN(354,3)] = eval_lobatto_edge (k,1,1,xi,eta,zeta) ;
        [N(355),dN(355,1),dN(355,2),dN(355,3)] = eval_lobatto_face (0,k,1,xi,eta,zeta) ;

% face functions
        n1 = 2 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(356:356+5,1) = NN ; dN(336:341,:) = dNN ; 
         
        n1 = 7 ; n2 = 2 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(362:362+5,1) = NN ; dN(336:341,:) = dNN ; 
        
        n1 = 7 ; n2 = 3 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(368:368+5,1) = NN ; dN(336:341,:) = dNN ; 
       

        n1 = 3 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(374:374+5,1) = NN ; dN(336:341,:) = dNN ; 

        n1 = 7 ; n2 = 4 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(380:380+5,1) = NN ; dN(336:341,:) = dNN ; 
               

        n1 = 4 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(386:386+5,1) = NN ; dN(336:341,:) = dNN ; 
                 

        n1 = 7 ; n2 = 5 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(392:392+5,1) = NN ; dN(336:341,:) = dNN ; 
        

        n1 = 5 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(398:398+5,1) = NN ; dN(336:341,:) = dNN ; 
       
        n1 = 7 ; n2 = 6 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(404:404+5,1) = NN ; dN(336:341,:) = dNN ; 
        

        n1 = 6 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(410:410+5,1) = NN ; dN(336:341,:) = dNN ; 
       
        n1 = 7 ; n2 = 7 ; 
        [NN,dNN] = eval_lobatto_face (n1,n2,xi,eta,zeta) ;
        N(416:416+5,1) = NN ; dN(336:341,:) = dNN ; 
          
% bubble function (needs to be changed)
% 2,3,4,5,6
% 226
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 336 
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 446
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 236     
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 326     
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 246
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 426   
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 346   
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 436   
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 622   
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 633   
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 644
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 623     
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 632     
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 624
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 642
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 634
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 643
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 262
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 363
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 464
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 263
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 362
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 264
        [N(422),dN(422,:)] = eval_lobatto_bubble (2,2,6,xi,eta,zeta) ;
% 462

% 364

% 463

% 266

% 366

% 466

% 626

% 636

% 646

% 662

% 663

% 664

% 556 

% 256

% 356

% 536

% 526

% 546

% 456

% 655

% 625

% 635

% 652

% 653

% 654

% 645

% 565

% 265

% 562

% 365

% 563

% 465

% 564

% 665

% 656

% 566


% 666


        end
   %}
    case 'T3'   
      N = [] ;
 if order >= 1 
    lambda1 =  ( eta +  1  )/2;     
    lambda2 = -( xi  + eta )/2 ;
    lambda3 =  ( xi  +  1  )/2 ;
    
    dlambda1dxi = 0 ;
    dlambda1deta = 1/2 ;
    dlambda2dxi = -1/2 ;
    dlambda2deta = -1/2 ;
    dlambda3dxi = 1/2 ;
    dlambda3deta = 0 ;

    f1 = (2*xi + eta + 1)/2 ; 
    f2 = (eta - xi)/2 ; 
    f3 = (-2*eta - xi - 1)/2 ; 
    df1dxi = 2/2 ; 
    df1deta = 1/2 ;
    df2dxi = -1/2; 
    df2deta = 1/2; 
    df3dxi = -1/2; 
    df3deta = -2/2; 
    
    
    N(1,1) = lambda2 ;
    N(2,1) = lambda3 ;
    N(3,1) = lambda1 ;

    dN(1,1) =  -1/2 ;
    dN(1,2) =  -1/2 ;
    dN(2,1) =   1/2 ;
    dN(2,2) =   0   ;
    dN(3,1) =   0   ;
    dN(3,2) =   1/2 ;

 end
      if order >=2 
        k = 2 ; 
% edge functions
        N(4,1) = lambda2*lambda3*kernel(k-2,f1) ;
        N(5,1) = lambda3*lambda1*kernel(k-2,f2) ; 
        N(6,1) = lambda1*lambda2*kernel(k-2,f3) ; 

        dN(4,1) =  dlambda2dxi*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3dxi*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1dxi ;
        
        dN(4,2) = dlambda2deta*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3deta*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1deta ; 
    
        dN(5,1) = dlambda3dxi*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1dxi*kernel(k-2,f2) +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2dxi ; 
    
        dN(5,2) = dlambda3deta*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1deta*kernel(k-2,f2)  +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2deta ; 
    
        dN(6,1) =  dlambda1dxi*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2dxi*kernel(k-2,f3) +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3dxi ;
    
        dN(6,2) =  dlambda1deta*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2deta*kernel(k-2,f3)  +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3deta ; 
        
% bubble function
%         N(9,1) = lobatto(2,xi)*lobatto(2,eta) ;
%         dN(9,1) = dlobatto(2,xi)*lobatto(2,eta) ;
%         dN(9,2) = lobatto(2,xi)*dlobatto(2,eta) ;


      end

      if order >= 3
        k = 3 ; 
% edge functions
        N(7,1)  = lambda2*lambda3*kernel(k-2,f1) ;
        N(8,1)  = lambda3*lambda1*kernel(k-2,f2) ; 
        N(9,1)  = lambda1*lambda2*kernel(k-2,f3) ; 

        dN(7,1) =  dlambda2dxi*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3dxi*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1dxi ;
        
        dN(7,2) = dlambda2deta*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3deta*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1deta ; 
    
        dN(8,1) = dlambda3dxi*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1dxi*kernel(k-2,f2) +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2dxi ; 
    
        dN(8,2) = dlambda3deta*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1deta*kernel(k-2,f2)  +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2deta ; 
    
        dN(9,1) =  dlambda1dxi*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2dxi*kernel(k-2,f3) +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3dxi ;
    
        dN(9,2) =  dlambda1deta*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2deta*kernel(k-2,f3)  +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3deta ;
               
        N(10,1)  = lambda1*lambda2*lambda3 ;
        dN(10,1) = dlambda1dxi*lambda2*lambda3  + lambda1*dlambda2dxi*lambda3 + lambda1*lambda2*dlambda3dxi  ;
        dN(10,2) = dlambda1deta*lambda2*lambda3  + lambda1*dlambda2deta*lambda3 + lambda1*lambda2*dlambda3deta  ;
%        N(10,1)  = lambda1*lambda2*lambda3*kernel(0,f1)*kernel(0,f3) ;
%         dN(10,1) = kernel(0,f1)*kernel(0,f3)*(dlambda1dxi*lambda2*lambda3  + lambda1*dlambda2dxi*lambda3 + lambda1*lambda2*dlambda3dxi)  ;
%         dN(10,2) = kernel(0,f1)*kernel(0,f3)*(dlambda1deta*lambda2*lambda3  + lambda1*dlambda2deta*lambda3 + lambda1*lambda2*dlambda3deta)  ;
      end
      if order >= 4 
        k = 4 ; 
% edge functions
        N(11,1)  = lambda2*lambda3*kernel(k-2,f1) ;
        N(12,1)  = lambda3*lambda1*kernel(k-2,f2) ; 
        N(13,1)  = lambda1*lambda2*kernel(k-2,f3) ; 

        dN(11,1) =  dlambda2dxi*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3dxi*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1dxi ;
        
        dN(11,2) = dlambda2deta*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3deta*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1deta ; 
    
        dN(12,1) = dlambda3dxi*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1dxi*kernel(k-2,f2) +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2dxi ; 
    
        dN(12,2) = dlambda3deta*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1deta*kernel(k-2,f2)  +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2deta ; 
    
        dN(13,1) =  dlambda1dxi*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2dxi*kernel(k-2,f3) +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3dxi ;
    
        dN(13,2) =  dlambda1deta*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2deta*kernel(k-2,f3)  +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3deta ;
               
        N(14,1)  = lambda1*lambda2*lambda3^2 ;
        dN(14,1) = dlambda1dxi*lambda2*lambda3^2  + lambda1*dlambda2dxi*lambda3^2 + lambda1*lambda2*2*lambda3*dlambda3dxi  ;
        dN(14,2) = dlambda1deta*lambda2*lambda3^2  + lambda1*dlambda2deta*lambda3^2 + lambda1*lambda2*2*lambda3*dlambda3deta  ;
        
        N(15,1)  = lambda1*lambda2^2*lambda3 ;
        dN(15,1) = dlambda1dxi*lambda2^2*lambda3  + lambda1*2*lambda2*dlambda2dxi*lambda3 + lambda1*lambda2^2*dlambda3dxi  ;
        dN(15,2) = dlambda1deta*lambda2^2*lambda3  + lambda1*2*lambda2*dlambda2deta*lambda3 + lambda1*lambda2^2*dlambda3deta  ;

      end
    case 'T4'   
      N = [] ;
    if order >= 1 
    l1 =  ( eta +  1  )/2;     
    l2 = -( 1 + xi  + eta + zeta)/2 ;
    l3 =  ( xi  +  1  )/2 ;
    l4 =  ( zeta +  1  )/2 ;
    
    dl1dxi = 0 ;
    dl1deta = 1/2 ;
    dl1dzeta = 0 ;

    dl2dxi = -1/2 ;
    dl2deta = -1/2 ;
    dl2dzeta = -1/2 ;
    
    dl3dxi = 1/2 ;
    dl3deta = 0 ;
    dl3dzeta = 0 ;
    
    dl4dxi = 0 ;
    dl4deta = 0 ;
    dl4dzeta = 1/2 ;
    


%     f1 = (2*xi + eta + 1)/2 ; 
%     f2 = (eta - xi)/2 ; 
%     f3 = (-2*eta - xi - 1)/2 ; 
%     df1dxi = 2/2 ; 
%     df1deta = 1/2 ;
%     df2dxi = -1/2; 
%     df2deta = 1/2; 
%     df3dxi = -1/2; 
%     df3deta = -2/2; 
   
    N(1,1) = l2 ;
    N(2,1) = l3 ;
    N(3,1) = l1 ;
    N(4,1) = l4 ;

    dN(1,1) =  dl2dxi ;
    dN(1,2) =  dl2deta ;
    dN(1,3) =  dl2dzeta ;

    dN(2,1) =  dl3dxi ;
    dN(2,2) =  dl3deta ;
    dN(2,3) =  dl3dzeta ;

    dN(3,1) =  dl1dxi ;
    dN(3,2) =  dl1deta ;
    dN(3,3) =  dl1dzeta ;
    
    dN(4,1) =  dl4dxi ;
    dN(4,2) =  dl4deta ;
    dN(4,3) =  dl4dzeta ;

    end
    
    if order >= 2 

    f1 = - ( 2*xi + eta + zeta + 2 ) / 2;  
    f2 = ( eta - xi ) / 2;  
    f3 = ( xi + 2*eta + zeta + 2 ) / 2;  
    f4 = -(xi + eta + 2*zeta + 2) / 2 ;  
    f5 = (xi + zeta ) / 2 ; 
    f6 = (eta - zeta ) / 2 ; 
    
%     f1 =-f1 ;
%     f2 = -f2 ; 
%     f3 = -f3 ;
%     f4 = -f4 ; 
%     f5 = -f5 ; 
%     f6 = -f6 ; 
    
    df1dxi = -2/2 ; 
    df1deta = -1/2; 
    df1dzeta = -1/2; 
    df2dxi = -1/2 ; 
    df2deta = 1/2; 
    df2dzeta = 0; 
    df3dxi = 1/2 ; 
    df3deta = 2/2; 
    df3dzeta = 1/2; 
    df4dxi = -1/2 ; 
    df4deta = -1/2; 
    df4dzeta = -2/2; 
    df5dxi = 1/2 ; 
    df5deta = 0; 
    df5dzeta = 1/2;
    df6dxi = 0 ; 
    df6deta = 1/2; 
    df6dzeta = -1/2;

    
    
    f7 = -f1 ; 
    df7dxi = -df1dxi ;    
    df7deta = -df1deta ;
    df7dzeta = -df1dzeta ;

    f8 = -f2 ; 
    df8dxi = -df2dxi ;    
    df8deta = -df2deta ;
    df8dzeta = -df2dzeta ;
    
    f9 = -f3 ;
    df9dxi = -df3dxi ;    
    df9deta = -df3deta ;
    df9dzeta = -df3dzeta ;
    
    f10 = -f4 ; 
    df10dxi = -df4dxi ;    
    df10deta = -df4deta ;
    df10dzeta = -df4dzeta ;
    
    f11 = -f5 ; 
    df11dxi = -df5dxi ;    
    df11deta = -df5deta ;
    df11dzeta = -df5dzeta ;
    
    f12 = -f6 ; 
    df12dxi = -df6dxi ;    
    df12deta = -df6deta ;
    df12dzeta = -df6dzeta ;
    

    
%     edge functions
    k = 2 ; 
    N(5,1) = l2*l3*kernel(k-2,f1) ; 
    N(6,1) = l1*l3*kernel(k-2,f2) ; 
    N(7,1) = l1*l2*kernel(k-2,f3) ; 
    N(8,1) = l2*l4*kernel(k-2,f4) ; 
    N(9,1) = l3*l4*kernel(k-2,f5) ; 
    N(10,1)= l1*l4*kernel(k-2,f6) ; 

    dN(5,1) =  dl2dxi*l3*kernel(k-2,f1) +...
                  l2*dl3dxi*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1dxi ;
              
    dN(5,2) =  dl2deta*l3*kernel(k-2,f1) +...
                  l2*dl3deta*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1deta ;
              
    dN(5,3) =  dl2dzeta*l3*kernel(k-2,f1) +...
                  l2*dl3dzeta*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1dzeta ;

    dN(6,1) =  dl1dxi*l3*kernel(k-2,f2) +...
                  l1*dl3dxi*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2dxi ;
              
    dN(6,2) =  dl1deta*l3*kernel(k-2,f2) +...
                  l1*dl3deta*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2deta ;
              
    dN(6,3) =  dl1dzeta*l3*kernel(k-2,f2) +...
                  l1*dl3dzeta*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2dzeta ;


    dN(7,1) =  dl1dxi*l2*kernel(k-2,f3) +...
                  l1*dl2dxi*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3dxi ;
              
    dN(7,2) =  dl1deta*l2*kernel(k-2,f3) +...
                  l1*dl2deta*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3deta ;
              
    dN(7,3) =  dl1dzeta*l2*kernel(k-2,f3) +...
                  l1*dl2dzeta*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3dzeta ;


    dN(8,1) =  dl2dxi*l4*kernel(k-2,f4) +...
                  l2*dl4dxi*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4dxi ;
              
    dN(8,2) =  dl2deta*l4*kernel(k-2,f4) +...
                  l2*dl4deta*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4deta ;
              
    dN(8,3) =  dl2dzeta*l4*kernel(k-2,f4) +...
                  l2*dl4dzeta*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4dzeta ;

    dN(9,1) =  dl3dxi*l4*kernel(k-2,f5) +...
                  l3*dl4dxi*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5dxi ;
              
    dN(9,2) =  dl3deta*l4*kernel(k-2,f5) +...
                  l3*dl4deta*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5deta ;
              
    dN(9,3) =  dl3dzeta*l4*kernel(k-2,f5) +...
                  l3*dl4dzeta*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5dzeta ;

    dN(10,1) =  dl1dxi*l4*kernel(k-2,f6) +...
                  l1*dl4dxi*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6dxi ;
              
    dN(10,2) =  dl1deta*l4*kernel(k-2,f6) +...
                  l1*dl4deta*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6deta ;
              
    dN(10,3) =  dl1dzeta*l4*kernel(k-2,f6) +...
                  l1*dl4dzeta*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6dzeta ;
              
    end
    
    if order >= 3 
%     edge functions
    k = 3 ; 
    N(11,1) = l2*l3*kernel(k-2,f1) ; 
    N(12,1) = l1*l3*kernel(k-2,f2) ; 
    N(13,1) = l1*l2*kernel(k-2,f3) ; 
    N(14,1) = l2*l4*kernel(k-2,f4) ; 
    N(15,1) = l3*l4*kernel(k-2,f5) ; 
    N(16,1)= l1*l4*kernel(k-2,f6) ; 

    dN(11,1) =  dl2dxi*l3*kernel(k-2,f1) +...
                  l2*dl3dxi*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1dxi ;
              
    dN(11,2) =  dl2deta*l3*kernel(k-2,f1) +...
                  l2*dl3deta*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1deta ;
              
    dN(11,3) =  dl2dzeta*l3*kernel(k-2,f1) +...
                  l2*dl3dzeta*kernel(k-2,f1) +...
                  l2*l3*dkernel(k-2,f1)*df1dzeta ;

    dN(12,1) =  dl1dxi*l3*kernel(k-2,f2) +...
                  l1*dl3dxi*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2dxi ;
              
    dN(12,2) =  dl1deta*l3*kernel(k-2,f2) +...
                  l1*dl3deta*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2deta ;
              
    dN(12,3) =  dl1dzeta*l3*kernel(k-2,f2) +...
                  l1*dl3dzeta*kernel(k-2,f2) +...
                  l1*l3*dkernel(k-2,f2)*df2dzeta ;


    dN(13,1) =  dl1dxi*l2*kernel(k-2,f3) +...
                  l1*dl2dxi*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3dxi ;
              
    dN(13,2) =  dl1deta*l2*kernel(k-2,f3) +...
                  l1*dl2deta*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3deta ;
              
    dN(13,3) =  dl1dzeta*l2*kernel(k-2,f3) +...
                  l1*dl2dzeta*kernel(k-2,f3) +...
                  l1*l2*dkernel(k-2,f3)*df3dzeta ;


    dN(14,1) =  dl2dxi*l4*kernel(k-2,f4) +...
                  l2*dl4dxi*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4dxi ;
              
    dN(14,2) =  dl2deta*l4*kernel(k-2,f4) +...
                  l2*dl4deta*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4deta ;
              
    dN(14,3) =  dl2dzeta*l4*kernel(k-2,f4) +...
                  l2*dl4dzeta*kernel(k-2,f4) +...
                  l2*l4*dkernel(k-2,f4)*df4dzeta ;

    dN(15,1) =  dl3dxi*l4*kernel(k-2,f5) +...
                  l3*dl4dxi*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5dxi ;
              
    dN(15,2) =  dl3deta*l4*kernel(k-2,f5) +...
                  l3*dl4deta*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5deta ;
              
    dN(15,3) =  dl3dzeta*l4*kernel(k-2,f5) +...
                  l3*dl4dzeta*kernel(k-2,f5) +...
                  l3*l4*dkernel(k-2,f5)*df5dzeta ;

    dN(16,1) =  dl1dxi*l4*kernel(k-2,f6) +...
                  l1*dl4dxi*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6dxi ;
              
    dN(16,2) =  dl1deta*l4*kernel(k-2,f6) +...
                  l1*dl4deta*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6deta ;
              
    dN(16,3) =  dl1dzeta*l4*kernel(k-2,f6) +...
                  l1*dl4dzeta*kernel(k-2,f6) +...
                  l1*l4*dkernel(k-2,f6)*df6dzeta ;
% face functions 
    n1 = 1 ; n2 = 1 ; 
    N(17,1) = l1*l2*l4*kernel(n1-1,f9)*kernel(n2-1,f6) ; 
    N(18,1) = l2*l3*l4*kernel(n1-1,f7)*kernel(n2-1,f4) ; 
    N(19,1) = l1*l3*l4*kernel(n1-1,f8)*kernel(n2-1,f6) ; 
    N(20,1) = l1*l2*l3*kernel(n1-1,f9)*kernel(n2-1,f2) ; 
    
    dN(17,1) =    dl1dxi*l2*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*dl2dxi*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*l2*dl4dxi*kernel(n1-1,f9)*kernel(n2-1,f6) +... 
                  l1*l2*l4*dkernel(n1-1,f9)*df9dxi*kernel(n2-1,f6) +...
                  l1*l2*l4*kernel(n1-1,f9)*dkernel(n2-1,f6)*df6dxi ;

    dN(17,2) =    dl1deta*l2*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*dl2deta*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*l2*dl4deta*kernel(n1-1,f9)*kernel(n2-1,f6) +... 
                  l1*l2*l4*dkernel(n1-1,f9)*df9deta*kernel(n2-1,f6) +...
                  l1*l2*l4*kernel(n1-1,f9)*dkernel(n2-1,f6)*df6deta ;
              
    dN(17,3) =    dl1dzeta*l2*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*dl2dzeta*l4*kernel(n1-1,f9)*kernel(n2-1,f6) +...
                  l1*l2*dl4dzeta*kernel(n1-1,f9)*kernel(n2-1,f6) +... 
                  l1*l2*l4*dkernel(n1-1,f9)*df9dzeta*kernel(n2-1,f6) +...
                  l1*l2*l4*kernel(n1-1,f9)*dkernel(n2-1,f6)*df6dzeta ;
    
              
    dN(18,1) =    dl2dxi*l3*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*dl3dxi*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*l3*dl4dxi*kernel(n1-1,f7)*kernel(n2-1,f4) +... 
                  l2*l3*l4*dkernel(n1-1,f7)*df7dxi*kernel(n2-1,f4) +...
                  l2*l3*l4*kernel(n1-1,f7)*dkernel(n2-1,f4)*df4dxi ;

    dN(18,2) =    dl2deta*l3*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*dl3deta*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*l3*dl4deta*kernel(n1-1,f7)*kernel(n2-1,f4) +... 
                  l2*l3*l4*dkernel(n1-1,f7)*df7deta*kernel(n2-1,f4) +...
                  l2*l3*l4*kernel(n1-1,f7)*dkernel(n2-1,f4)*df4deta ;
              
    dN(18,3) =    dl2dzeta*l3*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*dl3dzeta*l4*kernel(n1-1,f7)*kernel(n2-1,f4) +...
                  l2*l3*dl4dzeta*kernel(n1-1,f7)*kernel(n2-1,f4) +... 
                  l2*l3*l4*dkernel(n1-1,f7)*df7dzeta*kernel(n2-1,f4) +...
                  l2*l3*l4*kernel(n1-1,f7)*dkernel(n2-1,f4)*df4dzeta ;
    
              
    dN(19,1) =    dl1dxi*l3*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*dl3dxi*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*l3*dl4dxi*kernel(n1-1,f8)*kernel(n2-1,f6) +... 
                  l1*l3*l4*dkernel(n1-1,f8)*df8dxi*kernel(n2-1,f6) +...
                  l1*l3*l4*kernel(n1-1,f8)*dkernel(n2-1,f6)*df6dxi ;

    dN(19,2)  =    dl1deta*l3*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*dl3deta*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*l3*dl4deta*kernel(n1-1,f8)*kernel(n2-1,f6) +... 
                  l1*l3*l4*dkernel(n1-1,f8)*df8deta*kernel(n2-1,f6) +...
                  l1*l3*l4*kernel(n1-1,f8)*dkernel(n2-1,f6)*df6deta ;
              
    dN(19,3)  =    dl1dzeta*l3*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*dl3dzeta*l4*kernel(n1-1,f8)*kernel(n2-1,f6) +...
                  l1*l3*dl4dzeta*kernel(n1-1,f8)*kernel(n2-1,f6) +... 
                  l1*l3*l4*dkernel(n1-1,f8)*df8dzeta*kernel(n2-1,f6) +...
                  l1*l3*l4*kernel(n1-1,f8)*dkernel(n2-1,f6)*df6dzeta ;
    
    dN(20,1) =    dl1dxi*l2*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*dl2dxi*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*l2*dl3dxi*kernel(n1-1,f9)*kernel(n2-1,f2) +... 
                  l1*l2*l3*dkernel(n1-1,f9)*df9dxi*kernel(n2-1,f2) +...
                  l1*l2*l3*kernel(n1-1,f9)*dkernel(n2-1,f2)*df2dxi ;

    dN(20,2)  =    dl1deta*l2*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*dl2deta*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*l2*dl3deta*kernel(n1-1,f9)*kernel(n2-1,f2) +... 
                  l1*l2*l3*dkernel(n1-1,f9)*df9deta*kernel(n2-1,f2) +...
                  l1*l2*l3*kernel(n1-1,f9)*dkernel(n2-1,f2)*df2deta ;
              
    dN(20,3) =    dl1dzeta*l2*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*dl2dzeta*l3*kernel(n1-1,f9)*kernel(n2-1,f2) +...
                  l1*l2*dl3dzeta*kernel(n1-1,f9)*kernel(n2-1,f2) +... 
                  l1*l2*l3*dkernel(n1-1,f9)*df9dzeta*kernel(n2-1,f2) +...
                  l1*l2*l3*kernel(n1-1,f9)*dkernel(n2-1,f2)*df2dzeta ;
        
              
    end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
end
% N
    I=eye(dim);
    Nv=[];
    for i=1:size(N,1)
        Nv=[Nv;I*N(i)];
    end
%{
%   if ( dim == 1 )
%     B=dNdxi;
%   elseif ( dim == 2 )
%     B=zeros(dim*size(N,1),3);
%     
%     B(1:dim:dim*size(N,1)-1,1) = dNdxi(:,1);
%     B(2:dim:dim*size(N,1),2)   = dNdxi(:,2);
%     
%     B(1:dim:dim*size(N,1)-1,3) = dNdxi(:,2);
%     B(2:dim:dim*size(N,1),3)   = dNdxi(:,1);
%   elseif ( dim == 3 )
%     B=zeros(dim*size(N,1),6);
%     
%     disp('Error: need to add 3D N and dNdxi')
%     
%     B(1:dim:dim*size(N,1)-2,1) = dNdxi(:,1);
%     B(2:dim:dim*size(N,1)-1,2) = dNdxi(:,2);
%     B(3:dim:dim*size(N,1),3)   = dNdxi(:,3);
%     
%     B(2:dim:dim*size(N,1)-1,4) = dNdxi(:,3);
%     B(3:dim:dim*size(N,1),4)   = dNdxi(:,2);
%     
%     B(3:dim:dim*size(N,1),5)   = dNdxi(:,1);
%     B(1:dim:dim*size(N,1)-2,5) = dNdxi(:,3);
%     
%     B(1:dim:dim*size(N,1)-2,6) = dNdxi(:,2);
%     B(2:dim:dim*size(N,1)-1,6) = dNdxi(:,1);
%     
%   end  
end
%}


% dlobatto
%% lobatto test
%     phi = [] ;
%     figure
%     hold on
%     for p = 0 : 10
% %     p = 2 ;
%     cc = 1 ;
%     for kk = -1 : 0.01 : 1 
%        phi(cc,:)  = [ kk  lobatto(p,kk) ]; 
%        cc = cc + 1; 
%     end
%     
%     plot(phi(:,1),phi(:,2),'-')
%     
% %     axis equal
% %     axis tight 
%     end
%% 2D lobatto test:
%{
QQ =[ ] ; ss = [] ;
for xi = -1:0.1/5:1 
    for eta = -1:0.1/5:1 
%         for zeta = -1:0.2:1 
            
        % vertex functions
        ff(1,1) = lobatto(0,xi)*lobatto(0,eta) ;
        ff(2,1) = lobatto(1,xi)*lobatto(0,eta) ;
        ff(3,1) = lobatto(1,xi)*lobatto(1,eta) ;
        ff(4,1) = lobatto(0,xi)*lobatto(1,eta) ;


        % edge functions
        k = 2 ; 
        ff(5,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(6,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(7,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(8,1) = lobatto(k,xi)*lobatto(1,eta) ;

% bubble function
        ff(9,1) = lobatto(2,xi)*lobatto(2,eta) ;


k = 3 ; 
% edge functions
        ff(10,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(11,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(12,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(13,1) = lobatto(k,xi)*lobatto(1,eta) ;

        % dN(10,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        % dN(10,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        % dN(11,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        % dN(11,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        % dN(12,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        % dN(12,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        % dN(13,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        % dN(13,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        ff(14,1) = lobatto(2,xi)*lobatto(3,eta) ;
        % dN(14,1) = dlobatto(2,xi)*lobatto(3,eta) ;
        % dN(14,2) = lobatto(2,xi)*dlobatto(3,eta) ;
        
        ff(15,1) = lobatto(3,xi)*lobatto(2,eta) ;
        % dN(15,1) = dlobatto(3,xi)*lobatto(2,eta) ;
        % dN(15,2) = lobatto(3,xi)*dlobatto(2,eta) ;
                
        ff(16,1) = lobatto(3,xi)*lobatto(3,eta) ;
        % dN(16,1) = dlobatto(3,xi)*lobatto(3,eta) ;
        % dN(16,2) = lobatto(3,xi)*dlobatto(3,eta) ;

        k = 4 ; 
% edge functions
        ff(17,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(18,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(19,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(20,1) = lobatto(k,xi)*lobatto(1,eta) ;

        % dN(17,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        % dN(17,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        % dN(18,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        % dN(18,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        % dN(19,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        % dN(19,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        % dN(20,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        % dN(20,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        ff(21,1) = lobatto(2,xi)*lobatto(4,eta) ;
        % dN(21,1) = dlobatto(2,xi)*lobatto(4,eta) ;
        % dN(21,2) = lobatto(2,xi)*dlobatto(4,eta) ;
        
        ff(22,1) = lobatto(4,xi)*lobatto(2,eta) ;
        % dN(22,1) = dlobatto(4,xi)*lobatto(2,eta) ;
        % dN(22,2) = lobatto(4,xi)*dlobatto(2,eta) ;
                
        ff(23,1) = lobatto(3,xi)*lobatto(4,eta) ;
        % dN(23,1) = dlobatto(3,xi)*lobatto(4,eta) ;
        % dN(23,2) = lobatto(3,xi)*dlobatto(4,eta) ;
        
        ff(24,1) = lobatto(4,xi)*lobatto(3,eta) ;
        % dN(24,1) = dlobatto(4,xi)*lobatto(3,eta) ;
        % dN(24,2) = lobatto(4,xi)*dlobatto(3,eta) ;
        
        ff(25,1) = lobatto(4,xi)*lobatto(4,eta) ;
        % dN(25,1) = dlobatto(4,xi)*lobatto(4,eta) ;
        % dN(25,2) = lobatto(4,xi)*dlobatto(4,eta) ;
        

        
        k = 5 ; 
% edge functions
        ff(26,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(27,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(28,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(29,1) = lobatto(k,xi)*lobatto(1,eta) ;

        % dN(26,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        % dN(26,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        % dN(27,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        % dN(27,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        % dN(28,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        % dN(28,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        % dN(29,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        % dN(29,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        ff(30,1) = lobatto(2,xi)*lobatto(5,eta) ;
        % dN(30,1) = dlobatto(2,xi)*lobatto(5,eta) ;
        % dN(30,2) = lobatto(2,xi)*dlobatto(5,eta) ;
        
        ff(31,1) = lobatto(5,xi)*lobatto(2,eta) ;
        % dN(31,1) = dlobatto(5,xi)*lobatto(2,eta) ;
        % dN(31,2) = lobatto(5,xi)*dlobatto(2,eta) ;
                
        ff(32,1) = lobatto(3,xi)*lobatto(5,eta) ;
        % dN(32,1) = dlobatto(3,xi)*lobatto(5,eta) ;
        % dN(32,2) = lobatto(3,xi)*dlobatto(5,eta) ;
        
        ff(33,1) = lobatto(5,xi)*lobatto(3,eta) ;
        % dN(33,1) = dlobatto(5,xi)*lobatto(3,eta) ;
        % dN(33,2) = lobatto(5,xi)*dlobatto(3,eta) ;
        
        ff(34,1) = lobatto(4,xi)*lobatto(5,eta) ;
        % dN(34,1) = dlobatto(4,xi)*lobatto(5,eta) ;
        % dN(34,2) = lobatto(4,xi)*dlobatto(5,eta) ;
        
        ff(35,1) = lobatto(5,xi)*lobatto(4,eta) ;
        % dN(35,1) = dlobatto(5,xi)*lobatto(4,eta) ;
        % dN(35,2) = lobatto(5,xi)*dlobatto(4,eta) ;
        
        ff(36,1) = lobatto(5,xi)*lobatto(5,eta) ;
        % dN(36,1) = dlobatto(5,xi)*lobatto(5,eta) ;
        % dN(36,2) = lobatto(5,xi)*dlobatto(5,eta) ;
        
        
        k = 6 ; 
% edge functions
        ff(37,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(38,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(39,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(40,1) = lobatto(k,xi)*lobatto(1,eta) ;

        % dN(37,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        % dN(37,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        % dN(38,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        % dN(38,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        % dN(39,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        % dN(39,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        % dN(40,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        % dN(40,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        ff(41,1) = lobatto(2,xi)*lobatto(6,eta) ;
        % dN(41,1) = dlobatto(2,xi)*lobatto(6,eta) ;
        % dN(41,2) = lobatto(2,xi)*dlobatto(6,eta) ;
        
        ff(42,1) = lobatto(6,xi)*lobatto(2,eta) ;
        % dN(42,1) = dlobatto(6,xi)*lobatto(2,eta) ;
        % dN(42,2) = lobatto(6,xi)*dlobatto(2,eta) ;
                
        ff(43,1) = lobatto(3,xi)*lobatto(6,eta) ;
        % dN(43,1) = dlobatto(3,xi)*lobatto(6,eta) ;
        % dN(43,2) = lobatto(3,xi)*dlobatto(6,eta) ;
        
        ff(44,1) = lobatto(6,xi)*lobatto(3,eta) ;
        % dN(44,1) = dlobatto(6,xi)*lobatto(3,eta) ;
        % dN(44,2) = lobatto(6,xi)*dlobatto(3,eta) ;
        
        ff(45,1) = lobatto(4,xi)*lobatto(6,eta) ;
        % dN(45,1) = dlobatto(4,xi)*lobatto(6,eta) ;
        % dN(45,2) = lobatto(4,xi)*dlobatto(6,eta) ;
        
        ff(46,1) = lobatto(6,xi)*lobatto(4,eta) ;
        % dN(46,1) = dlobatto(6,xi)*lobatto(4,eta) ;
        % dN(46,2) = lobatto(6,xi)*dlobatto(4,eta) ;
        
        ff(47,1) = lobatto(5,xi)*lobatto(6,eta) ;
        % dN(47,1) = dlobatto(5,xi)*lobatto(6,eta) ;
        % dN(47,2) = lobatto(5,xi)*dlobatto(6,eta) ;
        
        ff(48,1) = lobatto(6,xi)*lobatto(5,eta) ;
        % dN(48,1) = dlobatto(6,xi)*lobatto(5,eta) ;
        % dN(48,2) = lobatto(6,xi)*dlobatto(5,eta) ;
        
        ff(49,1) = lobatto(6,xi)*lobatto(6,eta) ;
        % dN(49,1) = dlobatto(6,xi)*lobatto(6,eta) ;
        % dN(49,2) = lobatto(6,xi)*dlobatto(6,eta) ;
        
    
        k = 7 ; 
% edge functions
        ff(50,1) = lobatto(0,xi)*lobatto(k,eta) ;
        ff(51,1) = lobatto(1,xi)*lobatto(k,eta) ;
        ff(52,1) = lobatto(k,xi)*lobatto(0,eta) ;
        ff(53,1) = lobatto(k,xi)*lobatto(1,eta) ;

        % dN(50,1) =  dlobatto(0,xi)*lobatto(k,eta) ;
        % dN(50,2) =  lobatto (0,xi)*dlobatto(k,eta) ;
        % dN(51,1) =  dlobatto(1,xi)*lobatto(k,eta) ;
        % dN(51,2) =  lobatto (1,xi)*dlobatto(k,eta) ;
        % dN(52,1) =  dlobatto(k,xi)*lobatto(0,eta) ;
        % dN(52,2) =  lobatto (k,xi)*dlobatto(0,eta) ;
        % dN(53,1) =  dlobatto(k,xi)*lobatto(1,eta) ;
        % dN(53,2) =  lobatto (k,xi)*dlobatto(1,eta) ;
        
% bubble function
        ff(54,1) = lobatto(2,xi)*lobatto(7,eta) ;
        % dN(54,1) = dlobatto(2,xi)*lobatto(7,eta) ;
        % dN(54,2) = lobatto(2,xi)*dlobatto(7,eta) ;

        ff(55,1) = lobatto(7,xi)*lobatto(2,eta) ;
        % dN(55,1) = dlobatto(7,xi)*lobatto(2,eta) ;
        % dN(55,2) = lobatto(7,xi)*dlobatto(2,eta) ;

        ff(56,1) = lobatto(3,xi)*lobatto(7,eta) ;
        % dN(56,1) = dlobatto(3,xi)*lobatto(7,eta) ;
        % dN(56,2) = lobatto(3,xi)*dlobatto(7,eta) ;

        ff(57,1) = lobatto(7,xi)*lobatto(3,eta) ;
        % dN(57,1) = dlobatto(7,xi)*lobatto(3,eta) ;
        % dN(57,2) = lobatto(7,xi)*dlobatto(3,eta) ;

        ff(58,1) = lobatto(4,xi)*lobatto(7,eta) ;
        % dN(58,1) = dlobatto(4,xi)*lobatto(7,eta) ;
        % dN(58,2) = lobatto(4,xi)*dlobatto(7,eta) ;
        
        ff(59,1) = lobatto(7,xi)*lobatto(4,eta) ;
        % dN(59,1) = dlobatto(7,xi)*lobatto(4,eta) ;
        % dN(59,2) = lobatto(7,xi)*dlobatto(4,eta) ;
        
        ff(60,1) = lobatto(5,xi)*lobatto(7,eta) ;
        % dN(60,1) = dlobatto(5,xi)*lobatto(7,eta) ;
        % dN(60,2) = lobatto(5,xi)*dlobatto(7,eta) ;
        
        ff(61,1) = lobatto(7,xi)*lobatto(5,eta) ;
        % dN(61,1) = dlobatto(7,xi)*lobatto(5,eta) ;
        % dN(61,2) = lobatto(7,xi)*dlobatto(5,eta) ;
        
        ff(62,1) = lobatto(6,xi)*lobatto(7,eta) ;
        % dN(62,1) = dlobatto(6,xi)*lobatto(7,eta) ;
        % dN(62,2) = lobatto(6,xi)*dlobatto(7,eta) ;
        
        ff(63,1) = lobatto(7,xi)*lobatto(6,eta) ;
        % dN(63,1) = dlobatto(7,xi)*lobatto(6,eta) ;
        % dN(63,2) = lobatto(7,xi)*dlobatto(6,eta) ;
        
        ff(64,1) = lobatto(7,xi)*lobatto(7,eta) ;
        % dN(64,1) = dlobatto(7,xi)*lobatto(7,eta) ;
        % dN(64,2) = lobatto(7,xi)*dlobatto(7,eta) ;


        
%             ff = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta);
            QQ = [QQ ;  xi eta ] ; 
%             ss = [ss ;  ff1 ff2 ff3 ff4 ff5 ff6 ff7 ff8  ff9];% ff10 ff11 ff12 ff13 ff14 ff15 ff16 ff17 ff18 ff19 ] ; 
            ss = [ss ; ff'] ;

%         end
    end
end

tri = delaunay(QQ(:,1),QQ(:,2));

%%
% stress = [ss(:,1) 1*ss(:,2) 1*ss(:,3) 1*ss(:,4) ss(:,5) 1*ss(:,6) 1*ss(:,7) 1*ss(:,8)] ;
stress =ss ; 
disp = [ss(:,4) ss(:,5) ];
VTKPostProcess_new(QQ,tri,1,'Tri3','stress',stress,0*disp)
!paraview stress.vtu&
%}     
%% 2D lobatto test: T3
%{
    p1 = [-1 -1];
    p2 = [1 -1];
    p3 = [-1 1];
    tri = [p1;p2;p3];
    % check if points are inside
    ff = [] ;
QQ =[ ] ; ss = [] ;
for xi = -1:0.1/5:1 
    for eta = -1:0.1/5:1 
%         for zeta = -1:0.2:1 
if (isPointInTriangle([xi eta], tri))            
        % vertex functions
lambda1 =  ( eta +  1  )/2;     
    lambda2 = -( xi  + eta )/2 ;
    lambda3 =  ( xi  +  1  )/2 ;
    
    dlambda1dxi = 0 ;
    dlambda1deta = 1/2 ;
    dlambda2dxi = -1/2 ;
    dlambda2deta = -1/2 ;
    dlambda3dxi = 1/2 ;
    dlambda3deta = 0 ;

    f1 = (2*xi + eta + 1)/2 ; 
    f2 = (eta - xi)/2 ; 
    f3 = (-2*eta - xi - 1)/2 ; 
    df1dxi = 2/2 ; 
    df1deta = 1/2 ;
    df2dxi = -1/2; 
    df2deta = 1/2; 
    df3dxi = -1/2; 
    df3deta = -2/2; 
    
    k = 3 ;
    ff(1) = lambda2 ;
    ff(2) = lambda3 ;
    ff(3) = lambda1 ;
    ff(4) = lambda2*lambda3*kernel(k-2,f1) ;
    ff(5) = lambda3*lambda1*kernel(k-2,f2) ; 
    ff(6) = lambda1*lambda2*kernel(k-2,f3) ; 
    
    ff(7)  = lambda1*lambda2*lambda3 ;
%     dN(10,1) = dlambda1dxi*lambda2*lambda3  + lambda1*dlambda2dxi*lambda3 + lambda1*lambda2*dlambda3dxi  ;
%     dN(10,2) = dlambda1deta*lambda2*lambda3  + lambda1*dlambda2deta*lambda3 + lambda1*lambda2*dlambda3deta  ;
    %      
        ff(8) =  dlambda2dxi*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3dxi*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1dxi ;
        
        ff(9) = dlambda2deta*lambda3*kernel(k-2,f1) +...
                  lambda2*dlambda3deta*kernel(k-2,f1) +...
                  lambda2*lambda3*dkernel(k-2,f1)*df1deta ; 
    
        ff(10) = dlambda3dxi*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1dxi*kernel(k-2,f2) +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2dxi ; 
    
        ff(11) = dlambda3deta*lambda1*kernel(k-2,f2) +...
                  lambda3*dlambda1deta*kernel(k-2,f2)  +...
                  lambda3*lambda1*dkernel(k-2,f2)*df2deta ; 
    
        ff(12) =  dlambda1dxi*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2dxi*kernel(k-2,f3) +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3dxi ;
    
        ff(13) =  dlambda1deta*lambda2*kernel(k-2,f3) + ... 
                   lambda1*dlambda2deta*kernel(k-2,f3)  +...
                   lambda1*lambda2*dkernel(k-2,f3)*df3deta ; 

               
    ss = [ss ; ff] ;
    QQ = [QQ ;  xi eta ] ; 

end
    end
end

tri = delaunay(QQ(:,1),QQ(:,2));

%%
% stress = [ss(:,1) 1*ss(:,2) 1*ss(:,3) 1*ss(:,4) ss(:,5) 1*ss(:,6) 1*ss(:,7) 1*ss(:,8)] ;
stress =ss ; 
disp = [ss(:,1) ss(:,2) ];
VTKPostProcess_new(QQ,tri,1,'Tri3','stress',stress,0*disp)
!paraview stress.vtu&
%}     
%% 3D lobatto test:
%{
QQ =zeros(length(rr)^3,3) ; ss = zeros(length(rr)^3,216) ;
rr = -1:0.1/2:1 ;
cc = 1 ;
for xi = rr
    xi
    for eta = rr
        for zeta = rr
            
% %        k = 2 ; 
% %         % edge functions
% %         ff1  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
% %         ff2 = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
% %         ff3 = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
% %         ff4 = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
% %         ff5 = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
% %         ff6 = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
% %         ff7 = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
% %         ff8 = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
% %         ff9 = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
% %         ff10 = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
% %         ff11 = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
% %         ff12 = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
% % 
% %         
% % % face functions
% %         n1 = 4 ; n2 = 3 ; 
% %         ff13 = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
% %         ff14 = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
% %         ff15 = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
% %         ff16 = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
% %         ff17 = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
% %         ff18 = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
% % 
% % % bubble function
% %         ff19 = lobatto(4,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
            

% if order >= 1 
        ff = [] ;  dff = [] ; 
        ff(1,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(2,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(3,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(4,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(5,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(6,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(7,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        ff(8,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        
        
        dff(1,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(1,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(1,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dff(2,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(2,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(2,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;
        
        dff(3,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(3,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(3,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;
        
        dff(4,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(4,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(4,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;
       
        dff(5,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(5,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(5,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(6,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(6,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(6,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(7,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(7,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(7,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
       
        dff(8,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(8,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(8,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
       
       % effd
        % if order >= 2 

        k = 2 ; 
        % edge fuffctioffs
        ff(9,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(10,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(11,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(12,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(13,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(14,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(15,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(16,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(17,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(18,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        ff(19,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        ff(20,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dff(9,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(9,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(9,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dff(10,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(10,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(10,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dff(11,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(11,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(11,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dff(12,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(12,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(12,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dff(13,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(13,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(13,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(14,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(14,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(14,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(15,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(15,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(15,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(16,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(16,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(16,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(17,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(17,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(17,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(18,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(18,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(18,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dff(19,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(19,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(19,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dff(20,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(20,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(20,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 2 ; 
        ff(21,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(22,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(23,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(24,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(25,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(26,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(21,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(21,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(21,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(22,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(22,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(22,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(23,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(23,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(23,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(24,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(24,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(24,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(25,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(25,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(25,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(26,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(26,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(26,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;

% bubble function
        ff(27,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dff(27,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dff(27,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dff(27,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

%     case 3
        % end
        
        % if order >= 3 

        k = 3 ; 
        % edge functions
        ff(28,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(29,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(30,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(31,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(32,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(33,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(34,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(35,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(36,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(37,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        ff(38,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        ff(39,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dff(28,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(28,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(28,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dff(29,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(29,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(29,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dff(30,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(30,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(30,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dff(31,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(31,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(31,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dff(32,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(32,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(32,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(33,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(33,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(33,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(34,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(34,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(34,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(35,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(35,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(35,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(36,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(36,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(36,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(37,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(37,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(37,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dff(38,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(38,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(38,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dff(39,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(39,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(39,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 3 ; 
        ff(40,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(41,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(42,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(43,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(44,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(45,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(40,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(40,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(40,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(41,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(41,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(41,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(42,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(42,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(42,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(43,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(43,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(43,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(44,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(44,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(44,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(45,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(45,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(45,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 3 ; n2 = 2 ; 
        ff(46,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(47,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(48,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(49,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(50,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(51,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(46,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(46,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(46,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(47,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(47,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(47,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(48,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(48,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(48,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(49,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(49,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(49,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(50,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(50,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(50,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(51,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(51,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(51,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 3 ; 
        ff(52,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(53,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(54,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(55,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(56,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(57,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(52,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(52,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(52,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(53,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(53,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(53,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(54,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(54,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(54,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(55,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(55,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(55,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(56,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(56,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(56,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(57,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(57,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(57,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function
        ff(58,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        ff(59,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        ff(60,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        ff(61,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        ff(62,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        ff(63,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        ff(64,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
%         ff(64,1)
        dff(58,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dff(58,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dff(58,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;

        dff(59,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dff(59,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dff(59,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;

        dff(60,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dff(60,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dff(60,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

        dff(61,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dff(61,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dff(61,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;

        dff(62,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dff(62,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dff(62,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;

        dff(63,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dff(63,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dff(63,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

        dff(64,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dff(64,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dff(64,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

       % end
   % if order >= 4 
        k = 4 ; 
        % edge functions
        ff(65,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(66,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(67,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(68,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(69,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(70,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(71,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(72,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(73,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(74,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        ff(75,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        ff(76,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dff(65,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(65,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(65,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dff(66,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(66,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(66,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dff(67,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(67,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(67,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dff(68,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(68,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(68,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dff(69,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(69,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(69,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(70,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(70,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(70,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(71,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(71,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(71,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(72,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(72,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(72,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(73,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(73,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(73,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(74,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(74,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(74,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dff(75,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(75,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(75,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dff(76,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(76,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(76,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 4 ; 
        ff(77,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(78,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(79,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(80,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(81,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(82,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(77,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(77,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(77,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(78,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(78,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(78,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(79,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(79,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(79,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(80,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(80,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(80,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(81,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(81,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(81,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(82,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(82,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(82,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 4 ; n2 = 2 ; 
        ff(83,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(84,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(85,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(86,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(87,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(88,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(83,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(83,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(83,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(84,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(84,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(84,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(85,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(85,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(85,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(86,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(86,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(86,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(87,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(87,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(87,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(88,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(88,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(88,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 3 ; n2 = 4 ; 
        ff(89,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(90,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(91,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(92,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(93,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(94,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(89,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(89,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(89,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(90,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(90,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(90,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(91,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(91,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(91,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(92,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(92,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(92,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(93,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(93,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(93,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(94,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(94,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(94,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 4 ; n2 = 3 ; 
        ff(95,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(96,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(97,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(98,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(99,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(100,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(95,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(95,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(95,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(96,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(96,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(96,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(97,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(97,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(97,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(98,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(98,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(98,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(99,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(99,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(99,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(100,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(100,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(100,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           

        n1 = 4 ; n2 = 4 ; 
        ff(101,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(102,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(103,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(104,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(105,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(106,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(101,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(101,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(101,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(102,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(102,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(102,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(103,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(103,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(103,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(104,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(104,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(104,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(105,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(105,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(105,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(106,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(106,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(106,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function
% 2,3,4
% 2 2 
        ff(107,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dff(107,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dff(107,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dff(107,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;


        ff(108,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dff(108,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dff(108,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dff(108,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;


        ff(109,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dff(109,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dff(109,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dff(109,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;

        ff(110,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dff(110,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dff(110,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dff(110,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;

        ff(111,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dff(111,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dff(111,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dff(111,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;

        ff(112,1)  =  lobatto(2,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dff(112,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dff(112,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dff(112,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;

        ff(113,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dff(113,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dff(113,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dff(113,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;


        ff(114,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dff(114,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dff(114,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dff(114,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;

        ff(115,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dff(115,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dff(115,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dff(115,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;

        ff(116,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dff(116,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dff(116,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dff(116,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;

        ff(117,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dff(117,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dff(117,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dff(117,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;

        ff(118,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dff(118,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dff(118,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dff(118,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;

        ff(119,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dff(119,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dff(119,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dff(119,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
        
        ff(120,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(4,zeta) ;        
        dff(120,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dff(120,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dff(120,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
        
        ff(121,1)  = lobatto(3,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dff(121,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dff(121,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dff(121,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;


        ff(122,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(2,zeta) ;        
        dff(122,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dff(122,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dff(122,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
        
        ff(123,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dff(123,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dff(123,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dff(123,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
        

        ff(124,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dff(124,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dff(124,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dff(124,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
        
        ff(125,1)  = lobatto(3,xi)*lobatto(4,eta)*lobatto(2,zeta) ;        
        dff(125,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dff(125,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dff(125,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
   % end
        
    % if order >= 5 
        k = 5 ; 
        % edge functions
        ff(126,1)  = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta) ;
        ff(127,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(128,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta) ;
        ff(129,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta) ;
        ff(130,1) = lobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(131,1) = lobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta) ;
        ff(132,1) = lobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(133,1) = lobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta) ;
        ff(134,1) = lobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta) ;
        ff(135,1) = lobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta) ;
        ff(136,1) = lobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta) ;
        ff(137,1) = lobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta) ;

        dff(126,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta)  ;
        dff(126,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(0,zeta)  ;
        dff(126,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(0,zeta)  ;

        dff(127,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(127,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(127,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;

        dff(128,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(0,zeta)  ;
        dff(128,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(0,zeta)  ;
        dff(128,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(0,zeta)  ;

        dff(129,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(0,zeta)  ;
        dff(129,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(0,zeta)  ;
        dff(129,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(0,zeta)  ;
        
        dff(130,1) =  dlobatto(0,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(130,2) =  lobatto(0,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(130,3) =  lobatto(0,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(131,1) =  dlobatto(1,xi)*lobatto(0,eta)*lobatto(k,zeta)  ;
        dff(131,2) =  lobatto(1,xi)*dlobatto(0,eta)*lobatto(k,zeta)  ;
        dff(131,3) =  lobatto(1,xi)*lobatto(0,eta)*dlobatto(k,zeta)  ;
        
        dff(132,1) =  dlobatto(1,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(132,2) =  lobatto(1,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(132,3) =  lobatto(1,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(133,1) =  dlobatto(0,xi)*lobatto(1,eta)*lobatto(k,zeta)  ;
        dff(133,2) =  lobatto(0,xi)*dlobatto(1,eta)*lobatto(k,zeta)  ;
        dff(133,3) =  lobatto(0,xi)*lobatto(1,eta)*dlobatto(k,zeta)  ;
        
        dff(134,1) =  dlobatto(k,xi)*lobatto(0,eta)*lobatto(1,zeta)  ;
        dff(134,2) =  lobatto(k,xi)*dlobatto(0,eta)*lobatto(1,zeta)  ;
        dff(134,3) =  lobatto(k,xi)*lobatto(0,eta)*dlobatto(1,zeta)  ;
        
        dff(135,1) =  dlobatto(1,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(135,2) =  lobatto(1,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(135,3) =  lobatto(1,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        dff(136,1) =  dlobatto(k,xi)*lobatto(1,eta)*lobatto(1,zeta)  ;
        dff(136,2) =  lobatto(k,xi)*dlobatto(1,eta)*lobatto(1,zeta)  ;
        dff(136,3) =  lobatto(k,xi)*lobatto(1,eta)*dlobatto(1,zeta)  ;
          
        dff(137,1) =  dlobatto(0,xi)*lobatto(k,eta)*lobatto(1,zeta)  ;
        dff(137,2) =  lobatto(0,xi)*dlobatto(k,eta)*lobatto(1,zeta)  ;
        dff(137,3) =  lobatto(0,xi)*lobatto(k,eta)*dlobatto(1,zeta)  ;
          
        
        
% face functions
        n1 = 2 ; n2 = 5 ; 
        ff(138,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(139,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(140,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(141,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(142,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(143,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(138,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(138,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(138,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(139,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(139,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(139,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(140,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(140,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(140,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(141,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(141,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(141,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(142,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(142,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(142,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(143,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(143,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(143,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
            
        n1 = 5 ; n2 = 2 ; 
        ff(144,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(145,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(146,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(147,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(148,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(149,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(144,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(144,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(144,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(145,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(145,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(145,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(146,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(146,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(146,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(147,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(147,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(147,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(148,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(148,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(148,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(149,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(149,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(149,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
        n1 = 5 ; n2 = 3 ; 
        ff(150,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(151,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(152,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(153,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(154,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(155,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(150,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(150,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(150,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(151,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(151,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(151,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(152,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(152,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(152,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(153,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(153,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(153,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(154,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(154,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(154,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(155,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(155,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(155,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
           
        n1 = 3 ; n2 = 5 ; 
        ff(156,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(157,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(158,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(159,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(160,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(161,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(156,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(156,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(156,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(157,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(157,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(157,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(158,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(158,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(158,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(159,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(159,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(159,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(160,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(160,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(160,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(161,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(161,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(161,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
       
        n1 = 5 ; n2 = 4 ; 
        ff(162,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(163,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(164,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(165,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(166,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(167,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(162,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(162,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(162,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(163,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(163,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(163,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(164,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(164,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(164,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(165,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(165,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(165,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(166,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(166,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(166,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(167,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(167,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(167,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         
        n1 = 4 ; n2 = 5 ; 
        ff(168,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(169,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(170,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(171,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(172,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(173,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(168,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(168,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(168,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(169,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(169,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(169,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(170,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(170,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(170,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(171,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(171,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(171,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(172,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(172,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(172,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(173,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(173,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(173,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
         



        n1 = 5 ; n2 = 5 ; 
        ff(174,1) = lobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(175,1) = lobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta) ;
        ff(176,1) = lobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta) ;
        ff(177,1) = lobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta) ;
        ff(178,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta) ;
        ff(179,1) = lobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta) ;
        
        dff(174,1) =  dlobatto(0,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(174,2) =  lobatto(0,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(174,3) =  lobatto(0,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(175,1) =  dlobatto(1,xi)*lobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(175,2) =  lobatto(1,xi)*dlobatto(n1,eta)*lobatto(n2,zeta)  ;
        dff(175,3) =  lobatto(1,xi)*lobatto(n1,eta)*dlobatto(n2,zeta)  ;
          
        dff(176,1) =  dlobatto(n1,xi)*lobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(176,2) =  lobatto(n1,xi)*dlobatto(0,eta)*lobatto(n2,zeta)  ;
        dff(176,3) =  lobatto(n1,xi)*lobatto(0,eta)*dlobatto(n2,zeta)  ;
  
        dff(177,1) =  dlobatto(n1,xi)*lobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(177,2) =  lobatto(n1,xi)*dlobatto(1,eta)*lobatto(n2,zeta)  ;
        dff(177,3) =  lobatto(n1,xi)*lobatto(1,eta)*dlobatto(n2,zeta)  ;
          
        dff(178,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(178,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(0,zeta)  ;
        dff(178,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(0,zeta)  ;
          
        dff(179,1) =  dlobatto(n1,xi)*lobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(179,2) =  lobatto(n1,xi)*dlobatto(n2,eta)*lobatto(1,zeta)  ;
        dff(179,3) =  lobatto(n1,xi)*lobatto(n2,eta)*dlobatto(1,zeta)  ;
          
% bubble function 
% 2,3,4,5
% 225
        ff(180,1) = lobatto(2,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dff(180,1) =  dlobatto(2,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dff(180,2) =  lobatto(2,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dff(180,3) =  lobatto(2,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;

% 335
        ff(181,1) = lobatto(3,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dff(181,1) =  dlobatto(3,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dff(181,2) =  lobatto(3,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dff(181,3) =  lobatto(3,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 445
        ff(182,1) = lobatto(4,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dff(182,1) =  dlobatto(4,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dff(182,2) =  lobatto(4,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dff(182,3) =  lobatto(4,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 235
        ff(183,1) = lobatto(2,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dff(183,1) =  dlobatto(2,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dff(183,2) =  lobatto(2,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dff(183,3) =  lobatto(2,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 325
        ff(184,1) = lobatto(3,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dff(184,1) =  dlobatto(3,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dff(184,2) =  lobatto(3,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dff(184,3) =  lobatto(3,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 245
        ff(185,1) = lobatto(2,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dff(185,1) =  dlobatto(2,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dff(185,2) =  lobatto(2,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dff(185,3) =  lobatto(2,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 425
        ff(186,1) = lobatto(4,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dff(186,1) =  dlobatto(4,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dff(186,2) =  lobatto(4,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dff(186,3) =  lobatto(4,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 345
        ff(187,1) = lobatto(3,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dff(187,1) =  dlobatto(3,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dff(187,2) =  lobatto(3,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dff(187,3) =  lobatto(3,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 435
        ff(188,1) = lobatto(4,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dff(188,1) =  dlobatto(4,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dff(188,2) =  lobatto(4,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dff(188,3) =  lobatto(4,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 522
        ff(189,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(2,zeta) ;
        dff(189,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(2,zeta)  ;
        dff(189,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(2,zeta)  ;
        dff(189,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(2,zeta)  ;
% 533
        ff(190,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(3,zeta) ;
        dff(190,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(3,zeta)  ;
        dff(190,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(3,zeta)  ;
        dff(190,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(3,zeta)  ;
% 544
        ff(191,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(4,zeta) ;
        dff(191,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(4,zeta)  ;
        dff(191,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(4,zeta)  ;
        dff(191,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(4,zeta)  ;
% 523
        ff(192,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(3,zeta) ;
        dff(192,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(3,zeta)  ;
        dff(192,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(3,zeta)  ;
        dff(192,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(3,zeta)  ;
% 532
        ff(193,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(2,zeta) ;
        dff(193,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(2,zeta)  ;
        dff(193,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(2,zeta)  ;
        dff(193,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(2,zeta)  ;
% 524
        ff(194,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(4,zeta) ;
        dff(194,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(4,zeta)  ;
        dff(194,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(4,zeta)  ;
        dff(194,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(4,zeta)  ;
% 542
        ff(195,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(2,zeta) ;
        dff(195,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(2,zeta)  ;
        dff(195,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(2,zeta)  ;
        dff(195,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(2,zeta)  ;
% 534
        ff(196,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(4,zeta) ;
        dff(196,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(4,zeta)  ;
        dff(196,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(4,zeta)  ;
        dff(196,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(4,zeta)  ;
% 543
        ff(197,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(3,zeta) ;
        dff(197,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(3,zeta)  ;
        dff(197,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(3,zeta)  ;
        dff(197,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(3,zeta)  ;
% 252
        ff(198,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dff(198,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dff(198,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dff(198,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 353
        ff(199,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dff(199,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dff(199,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dff(199,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 454
        ff(200,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dff(200,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dff(200,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dff(200,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 253
        ff(201,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dff(201,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dff(201,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dff(201,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 352
        ff(202,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dff(202,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dff(202,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dff(202,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 254
        ff(203,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dff(203,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dff(203,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dff(203,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 452
        ff(204,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dff(204,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dff(204,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dff(204,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 354
        ff(205,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dff(205,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dff(205,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dff(205,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 453
        ff(206,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dff(206,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dff(206,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dff(206,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 255
        ff(207,1) = lobatto(2,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dff(207,1) =  dlobatto(2,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dff(207,2) =  lobatto(2,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dff(207,3) =  lobatto(2,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 355
        ff(208,1) = lobatto(3,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dff(208,1) =  dlobatto(3,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dff(208,2) =  lobatto(3,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dff(208,3) =  lobatto(3,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 455
        ff(209,1) = lobatto(4,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dff(209,1) =  dlobatto(4,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dff(209,2) =  lobatto(4,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dff(209,3) =  lobatto(4,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;
% 525
        ff(210,1) = lobatto(5,xi)*lobatto(2,eta)*lobatto(5,zeta) ;
        dff(210,1) =  dlobatto(5,xi)*lobatto(2,eta)*lobatto(5,zeta)  ;
        dff(210,2) =  lobatto(5,xi)*dlobatto(2,eta)*lobatto(5,zeta)  ;
        dff(210,3) =  lobatto(5,xi)*lobatto(2,eta)*dlobatto(5,zeta)  ;
% 535
        ff(211,1) = lobatto(5,xi)*lobatto(3,eta)*lobatto(5,zeta) ;
        dff(211,1) =  dlobatto(5,xi)*lobatto(3,eta)*lobatto(5,zeta)  ;
        dff(211,2) =  lobatto(5,xi)*dlobatto(3,eta)*lobatto(5,zeta)  ;
        dff(211,3) =  lobatto(5,xi)*lobatto(3,eta)*dlobatto(5,zeta)  ;
% 545
        ff(212,1) = lobatto(5,xi)*lobatto(4,eta)*lobatto(5,zeta) ;
        dff(212,1) =  dlobatto(5,xi)*lobatto(4,eta)*lobatto(5,zeta)  ;
        dff(212,2) =  lobatto(5,xi)*dlobatto(4,eta)*lobatto(5,zeta)  ;
        dff(212,3) =  lobatto(5,xi)*lobatto(4,eta)*dlobatto(5,zeta)  ;
% 552
        ff(213,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(2,zeta) ;
        dff(213,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(2,zeta)  ;
        dff(213,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(2,zeta)  ;
        dff(213,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(2,zeta)  ;
% 553
        ff(214,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(3,zeta) ;
        dff(214,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(3,zeta)  ;
        dff(214,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(3,zeta)  ;
        dff(214,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(3,zeta)  ;
% 554
        ff(215,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(4,zeta) ;
        dff(215,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(4,zeta)  ;
        dff(215,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(4,zeta)  ;
        dff(215,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(4,zeta)  ;
% 555
        ff(216,1) = lobatto(5,xi)*lobatto(5,eta)*lobatto(5,zeta) ;
        dff(216,1) =  dlobatto(5,xi)*lobatto(5,eta)*lobatto(5,zeta)  ;
        dff(216,2) =  lobatto(5,xi)*dlobatto(5,eta)*lobatto(5,zeta)  ;
        dff(216,3) =  lobatto(5,xi)*lobatto(5,eta)*dlobatto(5,zeta)  ;

   % end




            
%             ff = lobatto(k,xi)*lobatto(0,eta)*lobatto(0,zeta);
            QQ(cc,:)= [  xi eta zeta] ; 
%             ss = [ss ;  ff1 ff2 ff3 ff4 ff5 ff6 ff7 ff8 ff9 ff10 ff11 ff12 ff13 ff14 ff15 ff16 ff17 ff18 ff19 ] ; 
            ss(cc,:) = [ff'] ; 
            cc = cc  + 1; 
        end
    end
end

tri = delaunay(QQ(:,1),QQ(:,2),QQ(:,3));

%%
% stress = [ss(:,13) 1*ss(:,14) 1*ss(:,15) 1*ss(:,16)] ;
stress = ss ; 
disp = [ss(:,4) ss(:,5) ];
VTKPostProcess_new(QQ,tri,1,'Tet4','stress',stress,0*disp)
!paraview stress.vtu&
%}     