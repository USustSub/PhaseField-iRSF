LL = dlmread('LD.out');
dt = LL(itime,3) ; 

% if Time/60/60/24/365 > 35 && h__ == 1
%        dt = 0.001 ;
%        dt = min([0.2/1*L_/maxVp , 5*1e6*4]) ;
        disp(['Time step is ' num2str(dt) ' sec ' num2str(dt/60/60/24/365) ' years ' ' max Vp is ' num2str(maxVp) ])
%% Newmark constants:
beta = 2;
gama = 1.5 ; 
theta = 2;
if (beta>=(0.25*(0.5+gama)^2))&&(theta>=0.5)&&(gama>=0.5)
%     disp('Newmark Constants are OK!')
else
    error('Newmark Constants are not OK!')
end
a0 = 1/(beta*dt^2);
a1 = gama / (beta*dt);
a2 = 1 / (beta*dt) ; 
a3 = (gama/beta)-1 ;
a4 = (1/2/beta)-1 ; 
a5 = dt*((gama/2/beta)-1);

if dynamic == 0
    a0 = 0 ; 
end

% h__ = 0 ; 
% end
