           
nDofs = numnode*3  ; 

K = sparse(nDofs,nDofs);            
f = zeros(nDofs,1);   


% Initialize Solution vectors
u   = zeros( nDofs + 0*n_Cr , 1 );
du  = zeros( nDofs + 0*n_Cr, 1         );
ddu = zeros( nDofs + 0*n_Cr, 1         );

        
% Initialize the Nodal variables increment
du( : ) = 0.0;
