% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lobatto_v2 basis functions
% Mohsen Goudarzi September 2016 - Delft 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nv,dN]=lobatto_basis_v2(type,coord,dim,order)

if ( nargin == 2 )
    dim=1;
end
       xi = coord (1) ; 
       eta = coord(2) ;

    if length(coord) == 3 
        zeta = coord(3) ; 
    end
    
    p = order ; 
%     p = [1:12]
    nd = 4+4*(p-1)+(p-1).^2 ;

switch type
case 'Q4'
    N = zeros(nd,1) ; 
    dN = zeros(nd,2) ; 

    
% Vertex functions ( p = 1 )
        N(1,1) = lobatto_v2(0,xi)*lobatto_v2(0,eta) ;
        N(2,1) = lobatto_v2(1,xi)*lobatto_v2(0,eta) ;
        N(3,1) = lobatto_v2(1,xi)*lobatto_v2(1,eta) ;
        N(4,1) = lobatto_v2(0,xi)*lobatto_v2(1,eta) ;
    
        dN(1,1) =  dlobatto_v2(0,xi)*lobatto_v2(0,eta) ;
        dN(1,2) =  lobatto_v2(0,xi)*dlobatto_v2(0,eta) ;
        dN(2,1) =  dlobatto_v2(1,xi)*lobatto_v2(0,eta) ;
        dN(2,2) =  lobatto_v2(1,xi)*dlobatto_v2(0,eta) ;
        dN(3,1) =  dlobatto_v2(1,xi)*lobatto_v2(1,eta) ;
        dN(3,2) =  lobatto_v2(1,xi)*dlobatto_v2(1,eta) ;
        dN(4,1) =  dlobatto_v2(0,xi)*lobatto_v2(1,eta) ;
        dN(4,2) =  lobatto_v2(0,xi)*dlobatto_v2(1,eta) ;

    for p = 2 : order 
%         disp([' p = ' num2str(p) ])
        E_num = 4 ; % num new edge funcs, always 4
        B_num = (p-1)*2-1 ; % num new bubble funcs
        p_ = p-1 ; 
        nd = 4+4*(p_-1)+(p_-1).^2 ;
        nE = nd+1:nd+E_num ;
        nB = nd+E_num+1:nd+E_num+B_num ;
        
% Edge functions ( p > 1 )
        N(nE(1),1) = lobatto_v2(0,xi)*lobatto_v2(p,eta) ;
        N(nE(2),1) = lobatto_v2(1,xi)*lobatto_v2(p,eta) ;
        N(nE(3),1) = lobatto_v2(p,xi)*lobatto_v2(0,eta) ;
        N(nE(4),1) = lobatto_v2(p,xi)*lobatto_v2(1,eta) ;

        dN(nE(1),1) =  dlobatto_v2(0,xi)*lobatto_v2(p,eta) ;
        dN(nE(1),2) =  lobatto_v2 (0,xi)*dlobatto_v2(p,eta) ;
        dN(nE(2),1) =  dlobatto_v2(1,xi)*lobatto_v2(p,eta) ;
        dN(nE(2),2) =  lobatto_v2 (1,xi)*dlobatto_v2(p,eta) ;
        dN(nE(3),1) =  dlobatto_v2(p,xi)*lobatto_v2(0,eta) ;
        dN(nE(3),2) =  lobatto_v2 (p,xi)*dlobatto_v2(0,eta) ;
        dN(nE(4),1) =  dlobatto_v2(p,xi)*lobatto_v2(1,eta) ;
        dN(nE(4),2) =  lobatto_v2 (p,xi)*dlobatto_v2(1,eta) ;

% Bubble functions ( p > 1 )
        i_ = 0 ; 
        for p1 = 2 : p
            for p2 = 2 : p
                
                if p1>p-1 || p2>p-1
%                     disp([num2str(p1) ' , ' num2str(p2)])
                    i_ = i_ + 1 ; 
                    N(nB(i_),1) = lobatto_v2(p1,xi)*lobatto_v2(p2,eta) ;
                    dN(nB(i_),1) = dlobatto_v2(p1,xi)*lobatto_v2(p2,eta) ;
                    dN(nB(i_),2) = lobatto_v2(p1,xi)*dlobatto_v2(p2,eta) ;
                end
            end
        end
    end
    
    
        
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