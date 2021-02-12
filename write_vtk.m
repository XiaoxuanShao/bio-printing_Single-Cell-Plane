function [ ]= write_vtk(nx,ny,dx,dy,istep,data)

 format long

 %-- open output file

 fname=sprintf('time_%d.vtk',istep);
 out =fopen(fname,'w');   %w´ú±íÐ´Èë

 nz=1;

 npoint =nx*ny;

 % start writing ASCII VTK file:

 % header of VTK file

 fprintf(out,'# vtk DataFile Version 3.0\n');
 fprintf(out,'fluid_state\n');
 fprintf(out,'ASCII\n');
 fprintf(out,'DATASET RECTILINEAR_GRID\n');

 %--- coords of grid points
 fprintf(out,'DIMENSIONS %d %d %d\n',nx,ny,nz);
 
 fprintf(out, 'X_COORDINATES %d float\n', nx);
  for i = 1:nx
      fprintf(out, '%d ',i);
  end
  fprintf(out, '\n');
  
  fprintf(out, 'Y_COORDINATES %d float\n', ny);
  
  for j = 1:ny
      fprintf(out, '%d ',j);
  end
  fprintf(out, '\n');
  
   fprintf(out, 'Z_COORDINATES %d float\n', nz);
   fprintf(out, '%d ',1);
   fprintf(out, '\n');
  
 fprintf(out,'POINT_DATA %d\n',npoint);

 %--- write grid point values:


 fprintf(out,'SCALARS CON float 1\n');
 fprintf(out,'LOOKUP_TABLE default\n');

 for i = 1:nx
    for j = 1:ny
%      	ii=(i-1)*nx+j;
      	fprintf(out,'%f\n',data(i,j));
    end
 end

 fclose(out);

 end %end function