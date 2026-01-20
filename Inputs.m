method = 'Lobatto'; % 'Lobatto' or 'Lagrange'
Const_method = 'Penalty' ; 
Penalty_ = 1 ; 



% Dimension of the domain 
D = 1 ; 
L = 1 ;
% Four corner points
pt1 = [0 0] ;
pt2 = [D 0] ;
pt3 = [D L] ;
pt4 = [0 L] ;


% Number of nodes along two directions
nnx = 550 ;
nny = 550 ;

% Material properties
E  = 210000 ;
nu = 0.3 ;
compliance_ ;
G = 2.7 ; 
u_pres = 0.01 ; 
NR.iter      = 30;
NR.tolerance = 1e-10;

%% phase field inputs:
a = L/4; % crack length
theta = 0.2; % crack angle 
n_Cr = 5000 ; % number of penalty points along crack 
w_fac = 1e7 ; % penatly factor for imposing phase field = 1 along crack
l_ = 1/600 ; % phase field length scale
kkk = 1e-14 ; % small value
n_itr = 4 ; 
n_order = 2 ;
P_order = 1 ; % 1 is Q4 FEM... can be 1,2,3,4,5,6,7,...
fast_ass =  0 ; % if ==1, fast assembly using C++

vec_ =  1; % if ==1, fast assembly using vectorized matlab. only works for p = 1 
if P_order>1
    vec_ = 0 ; 
end

% crack geometry
xCr = [ L/2-a*cos(theta) L/2-a*sin(theta) ; L/2+a*cos(theta) L/2+a*sin(theta)];
xCr(:,2) = xCr(:,2)  ; % + .0012 ; 
    
%% mesh generation


elemType = 'Q4' ;
[node,element] = meshRectangularRegion(...
    pt1, pt2, pt3, pt4, nnx,nny,elemType);
[conn_nodes ] = find_conn_nodes (node , element ) ;


% find cracked elements
% crack_nodes = find(node(:,2)==0.5 & node(:,1)<=0.5);
if strcmp(Const_method,'Penalty') || strcmp(Const_method,'Lagrange')
    crack_elems ;
else
crack_nodes = find(node(:,2)==0.5 & node(:,1)>=0.25 & node(:,1)<=0.75);
desir_elems =conn_nodes(crack_nodes,2:end);
desir_elems = unique(desir_elems(:));
desir_elems(desir_elems==0) = [] ; 
end

if P_order>1
%     element_des = [ 5:7 32:35  ] ; 
%     element_des = [1:1:length(element)] ; 
    element_des = desir_elems ; 
%     [node_app,element_app]  =  P1toP2_func ( node , element , method );
    [node_app,element_app]  =  P1toP2_func_v3 ( node , element , method , P_order , element_des );
elseif P_order == 1 
    node_app = node ; 
    element_app = element ; 
end

    
% compute number of nodes, of elements
numnode = size(node_app,1);
numnode0 = size(node,1);
numelem = size(element,1);





%% define essential boundaries
uln = nnx*(nny-1)+1;       % upper left node number
urn = nnx*nny;             % upper right node number
lrn = nnx;                 % lower right node number
lln = 1;                   % lower left node number
cln = nnx*(nny-1)/2+1;     % node number at (0,0)

rightEdge= [ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
leftEdge = [ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';

% GET NODES ON DIRICHLET BOUNDARY AND ESSENTIAL BOUNDARY
botNodes   = unique(botEdge);     
rightNodes = unique(rightEdge);    
topNodes   = unique(topEdge);     
leftNodes  = unique(leftEdge);     

dispNodes = [leftNodes];
dispNodes = unique(dispNodes);

elem_center = nnx/2:nnx-1:length(element) ; 

figure
hold on
plot_mesh(node,element(desir_elems,:),'Q4','r-');
plot(node(:,1),node(:,2),'r.')
plot(xCr(:,1),xCr(:,2),'b-','LineWidth',2)
% plot(crack_nodes(:,1),crack_nodes(:,2),'ksq')
plot(node(topNodes,1),node(topNodes,2),'bsq')
plot(node(botNodes,1),node(botNodes,2),'bsq')
plot(node(lrn,1),node(lrn,2),'csq')
axis equal
% plot(node(split_nodes,1),node(split_nodes,2),'csq')
% plot(node(crack_nodes,1),node(crack_nodes,2),'ksq')
plot_mesh(node,element(elem_center,:),'Q4','b-');

%% load shape function
% if P_order>1
    fileName_ = ['data_p' num2str(P_order) '_nq' num2str(n_order) '.mat'];
    if isfile(fileName_)==0
        store_shape_functions ;
    end
eval(['load data_p' num2str(P_order) '_nq' num2str(n_order) '.mat N_p'])
% end