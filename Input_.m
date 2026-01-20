% Inputs :: 
L = 150000 ; %mm 
D = 150000 ; %mm
pt1 = [0 0]; pt2 = [L 0];  pt3 = [L D]; pt4 = [0 D];
nnx = nn ;  nny = nn ; nne = 4 ;   nni = 4 ;   ndf = 2 ; 
b_b = 1;           % width of bricks (and mortar joints) [mm]
lin_quad = 1;       % 1 = bilinear elements (4 nodes/element); 2 = quadratic elements (8 nodes/element)
rho = 2700 ; % kg/m3

% Material properties
E_b = 36000.5e3;        % Young's modulus bricks [N/mm^2]
v_b = 0.25;          % Poisson ratio bricks [-]
ft = 1*3.19 ; 
Gf = 1*0.05 ;

% phase field length scale
l_ =   300; 

G = 30e9;
K = 50e9;

E_b = 9*K*G/(3*K+G);
v_b = (3*K-2*G)/2/(3*K+G);

% [cmat_b] = cmatrix(E_b,v_b); % plane stress constitutive matrix of continua (bricks)
C1 = E_b*(1-v_b)/(1+v_b)/(1-2*v_b);                                         % Constant for elastic constant matrix
C2 = E_b*v_b/(1+v_b)/(1-2*v_b);                                             % Constant for elastic constant matrix
C3 = E_b/2/(1+v_b);                                                       % Constant for elastic constant matrix
Cm  = [C1 C2  0;...
       C2 C1  0;...
        0  0 C3];
cmat_b = Cm  ; 

% Integration point configuration 
[rsIP_b] = integrationpoints(2,1);      % continua: 2x2 (four-node element) or 3x3 (eight-node element) integration scheme (Gauss)
[rsIP_i] = integrationpoints(2,2);      % interfaces: 2 (four-node interface) or 3 (six-node interface) integration scheme (Newton-Cotes)
% [W,Q]=quadrature(4,'GAUSS',1);
% rsIP_i = [ Q W ];    


% generate ordinary nodes 
node = square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
inc_u = 1; inc_v = nnx ; node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v); 
node = round(node*1000000000)/1000000000 ;

% 1151 951 751 551 251 51
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
densy = [30 0.49 ; 151 0.02 ; 30 0.49];
densx = [180 1];
% densy = [30 0.49 ; 151 0.02 ; 30 0.49];
% densx = [110 1];
[node] = square_node_array_ir(pt1,pt2,pt3,pt4,densx,densy);
nnx = 1 + sum(densx(:,1)) ;
nny = 1 + sum(densy(:,1)) ;
inc_u = 1;
inc_v = nnx;
node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% find initial cracked elements
xCr = [ 0   D/2
        L+100 D/2] ; 
[split_elem , tip_elem , split_elem_all] = find_cracked_elems ( node , element,  xCr) ;
nonCracked_elems = setdiff(split_elem_all , [split_elem tip_elem]) ; 


%% generate connectivity for the initial crack

% node_init = node ; element_init = element ; 
% [node,element,conint]  = correct_element_group ...
%     ( node_init , element_init ,  split_elem , nonCracked_elems , tip_elem , []);

conint = [ ] ;
for iielem = 1 : length(split_elem)
    elem_num = split_elem(iielem) ; 
    sctr = element ( elem_num , : ) ;
    conint ( iielem , : )  =  sctr ; 
end















