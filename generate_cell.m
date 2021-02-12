function [ncell,phis,vac] = generate_cell(Nx,Ny,ncell,R)

 format long;
 
 %--- initialize
    for icell =1:ncell
        for i=1:Nx
            for j=1:Ny
                phis(i,j,icell) = 0;
            end
        end
    end

    phis =zeros(Nx,Ny,ncell);
    
    R2 = R*R;

    xc=60;
    yc=100;

    icell =2;    
    for i=1:Nx
        for j=1:Ny
            if((i-xc)*(i-xc) + (j-yc)*(j-yc)< R2)
                phis(i,j,icell) =0.999;
            end
        end
    end


    icell=1;
    for i=Nx/2:Nx/2+10
        for j=1:Ny
            phis(i,j,icell) =0.999;
        end
    end


     %--- cell self propulsion velocity:

     vac(1) = 0;
     vac(2) = 0.8;
 end %end function