function [nmat,xyIP_i,detJ] = nmatrix(rsIP_i,coordloc,lin_quad)

xyIP_i = zeros(size(rsIP_i,1),2);
detJ = zeros(size(rsIP_i,1),2);

if (lin_quad == 1)  % four-node interfaces
    
    nmat = zeros(2,8,size(rsIP_i,1));
    
    for IP = 1:size(rsIP_i,1)
        
        r = rsIP_i(IP,1);
        x1 = coordloc(1,1);
        y1 = coordloc(1,2);
        x2 = coordloc(2,1);
        y2 = coordloc(2,2);
                
        N1 = 0.5*(1-r);
        N2 = 0.5*(1+r);        
        
                       %    x1    y1    |  x2     y2     |  x3     y3   |  x4     y4
        nmat(:,:,IP) = [    N1    0        N2    0          -N2    0       -N1    0    ;...
                            0     N1       0     N2         0      -N2     0      -N1  ];

        xyIP_i(IP,1) = N1*x1 + N2*x2;
        xyIP_i(IP,2) = N1*y1 + N2*y2;
        
        length_i = sqrt((y2-y1)^2 + (x2-x1)^2);
        
        detJ(IP,1) = (1/2)*length_i;
        detJ(IP,2) = rsIP_i(IP,2);
        
    end
elseif  (lin_quad == -1000)  % four-node interfaces
    
    nmat = zeros(1,4,size(rsIP_i,1));
    
    for IP = 1:size(rsIP_i,1)
        
        r = rsIP_i(IP,1);
        x1 = coordloc(1,1);
        y1 = coordloc(1,2);
        x2 = coordloc(2,1);
        y2 = coordloc(2,2);
                
        N1 = 0.5*(1-r);
        N2 = 0.5*(1+r);        
        
                       %    x1    y1    |  x2     y2     |  x3     y3   |  x4     y4
        nmat(:,:,IP) = [    N1      N2      -N2     -N1    ];

        xyIP_i(IP,1) = N1*x1 + N2*x2;
        xyIP_i(IP,2) = N1*y1 + N2*y2;
        
        length_i = sqrt((y2-y1)^2 + (x2-x1)^2);
        
        detJ(IP,1) = (1/2)*length_i;
        detJ(IP,2) = rsIP_i(IP,2);
        
    end
else  % six-node interfaces
    
    nmat = zeros(2,12,size(rsIP_i,1));
    
    for IP = 1:size(rsIP_i,1)
        
        r = rsIP_i(IP,1);
        x1 = coordloc(1,1);
        y1 = coordloc(1,2);
        x2 = coordloc(2,1);
        y2 = coordloc(2,2);
        x3 = coordloc(3,1);
        y3 = coordloc(3,2);
        
        N1 = 0.5*(1-r) - 0.5*(1-(r^2));
        N2 = 0.5*(1+r) - 0.5*(1-(r^2));
        N3 = 1-(r^2);
        
                       %    x1    y1    |  x2     y2     |  x3     y3   |  x4     y4    |  x5     y5   |  x6     y6
        nmat(:,:,IP) = [    N1    0        N2    0         -N2    0        -N1    0        N3    0       -N3    0;...
                            0     N1       0     N2        0     -N2        0    -N1        0    N3       0    -N3];
        
        xyIP_i(IP,1) = N1*x1 + N2*x2 + N3*x3;
        xyIP_i(IP,2) = N1*y1 + N2*y2 + N3*y3;
        
        lengthint = sqrt((y2-y1)^2 + (x2-x1)^2);
        
        detJ(IP,1) = (1/2)*lengthint;
        detJ(IP,2) = rsIP_i(IP,2);
        
    end
    
end