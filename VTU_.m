function VTU_ (xx,yy,V,output)

sxx = V(:,1);
syy = V(:,2);
sxy = V(:,3);
ux = V(:,4);
uy = V(:,5);
phi = V(:,6);


if size(yy,2) == 2 % L2 elements
    offset = [2:2:2*size(yy,1)]';     
    type = yy(:,1)*0+3;
elseif size(yy,2) == 3 % T3 element
    offset = [3:3:3*size(yy,1)]';     
    type = yy(:,1)*0+5;
elseif size(yy,2) == 4 % Q4 element
    offset = [4:4:4*size(yy,1)]';     
    type = yy(:,1)*0+9;
end

fid = fopen(output,'wt');

% preamble
fprintf(fid, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n ');
fprintf(fid, '<UnstructuredGrid> \n ');
fprintf(fid,'<Piece  NumberOfPoints=" %6.0f " NumberOfCells=" %6.0f "> \n ', size(xx,1), size(yy,1) );

% nodal data
fprintf(fid, '<Points>  \n ');
fprintf(fid, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" >   \n ');
fprintf(fid,'%8.5f %8.5f %8.5f \n', xx' );
fprintf(fid, '</DataArray>  \n ');
fprintf(fid, '</Points>  \n ');

% Field values at nodes
fprintf(fid, '<PointData  Vectors="fields"> \n ');
% ux
fprintf(fid, '<DataArray  type="Float64"  Name="ux" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', ux' );
fprintf(fid, '</DataArray>  \n ');
% uy
fprintf(fid, '<DataArray  type="Float64"  Name="uy" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', uy' );
fprintf(fid, '</DataArray>  \n ');
% phi
fprintf(fid, '<DataArray  type="Float64"  Name="phi" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', phi' );
fprintf(fid, '</DataArray>  \n ');
% % sxx
fprintf(fid, '<DataArray  type="Float64"  Name="sxx" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', sxx' );
fprintf(fid, '</DataArray>  \n ');
% sxy
fprintf(fid, '<DataArray  type="Float64"  Name="sxy" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', sxy' );
fprintf(fid, '</DataArray>  \n ');
% syy
fprintf(fid, '<DataArray  type="Float64"  Name="syy" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%8.5f \n', syy' );
fprintf(fid, '</DataArray>  \n ');



fprintf(fid, '</PointData>  \n ');

   


% connectivity data 
fprintf(fid, '<Cells>  \n ');
fprintf(fid, '<DataArray  type="Float64"  Name="connectivity"  format="ascii">  \n ');
fprintf(fid,'%8.5f %8.5f  \n', yy'-1 );
fprintf(fid, '</DataArray>  \n ');
% offsett
fprintf(fid, '<DataArray  type="Float64"  Name="offsets"  format="ascii">  \n ');
fprintf(fid,'%8.5f   \n', offset' );
fprintf(fid, '</DataArray>  \n ');
% type
fprintf(fid, '<DataArray  type="Float64"  Name="types"  format="ascii">  \n ');
fprintf(fid,'%8.5f   \n', type' );
fprintf(fid, '</DataArray>  \n ');
fprintf(fid, '</Cells>  \n ');
fprintf(fid, '</Piece>  \n ');
fprintf(fid, '</UnstructuredGrid>  \n ');
fprintf(fid, '</VTKFile>  \n ');
