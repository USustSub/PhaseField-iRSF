%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Phase-field continuum model for faults with RSF  
% By: Mohsen Goudarzi (2023, Utrecht), Goudarzi.mohsen@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all 

! rm -r out_results/
! mkdir out_results/

L = 150000 ; 
% different length scales
for l_ = [L/500 L/1000 L/3000 L/6000]

% different distance values
for eps_ = [ 0.5 1 2 3 4]
    
% this works under linux
eval(['!sed -i ',char(39),'17s/.*/l_ =   ',num2str(l_),...
    '; /',char(39),' Input_.m']);
eval(['!sed -i ',char(39),'201s/.*/    eps_ =  ',num2str(eps_*l_),...
    '; /',char(39),' Interface_friction_phasefield_final.m']);
    
save ll.mat l_ eps_ 

% main function call
Interface_friction_phasefield_final ; 

load ll.mat l_ eps_ 

% copy results 
eval(['! cp -r out out_results/out',num2str(l_),'_',num2str(eps_)]);


end

end

