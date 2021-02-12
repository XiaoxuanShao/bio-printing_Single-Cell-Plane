%get intial wall time:
time0=clock();
format long;
%Simulation cell parameters:

 Nx   = 200;
 Ny   = 200;
 NxNy = Nx*Ny;

 dx   = 1.0;
 dy   = 1.0;

 %--- Time integration parameters:

 nstep  = 100000;
 nprint = 50;
 dtime  = 5.e-3;

 %--- Material specific Parameters:
 ncell = 2;
 R     = 30.0;
 
 gamma = 22.0;
 lamda = 7.0;
 kappa = 60;
 mu    = 40.0;
 kisa  = 1.5e3;


 pix    = 4.0  * atan(1.0);
 const1 = 30.0 / lamda^2;
 const2 = 2.0  * mu/(pix*R^2);
 const3 = 60.0 * kappa/(lamda^2*kisa);

 %---
 %prepare initial microstructure:
 %---
 
 [ncell,phis,vac] = generate_cell(Nx,Ny,ncell,R);

 %---
 %- Evolution
 %---

 for istep =1:nstep
     
     if(istep <=500 )
      for icell=1:ncell
     gamma_cell = gamma;
         
         for i=1:Nx
             for j=1:Ny
                 phi(i,j) = phis(i,j,icell);
             end
         end
         
         %----
 % calculate the laplacian and gradient terms:
 %---
         for i=1:Nx
            for j=1:Ny
                jp=j+1;
                jm=j-1;

                ip=i+1;
                im=i-1;

                jp=j+1;
                jm=j-1;

                ip=i+1;
                im=i-1;

                if(im == 0)
                    im=Nx;
                end
                
                if(ip == (Nx+1))
                    ip=1;
                end
                
                if(jm == 0)
                    jm = Ny;
                end
                
                if(jp == (Ny+1))
                    jp=1;
                end

                hne=phi(ip,j);
                hnw=phi(im,j);
                hns=phi(i,jm);
                hnn=phi(i,jp);
                hnc=phi(i,j);

                lap_phi(i,j) =(hnw + hne + hns+ hnn -4.0*hnc)/(dx*dy);

    %--gradients of phi:
                phidx(i,j) = (phi(ip,j) - phi(im,j))/dx;
                phidy(i,j) = (phi(i,jp) - phi(i,jm))/dy;
            end %for j
        end %for i

 %----
 %--- volum integrations:
        vinteg = 0.0;
        vintegx = 0.0;
        vintegy = 0.0;

        for i=1:Nx
            for j=1:Ny
                vinteg = vinteg + phi(i,j)^2;
                sum_phi = 0.0;
                for jcell =1:ncell
                    if(icell ~= jcell)
                        sum_phi = sum_phi + phis(i,j,jcell)^2;
                    end
                end
                
                vintegx = vintegx + phi(i,j)* phidx(i,j)*sum_phi;
                vintegy = vintegy + phi(i,j)*phidy(i,j)*sum_phi;
            end
        end
        
        for i=1:Nx
            for j=1:Ny
 %-- Second term:
 %-- derivative of free energy
                [dfdphi]=free_energy_v1(i,j,icell,ncell,gamma,kappa,phi,phis);
                term2 = -const1*dfdphi;

 %---
 %---Third term
                term3 = -const2*(vinteg-pix*R^2)*phi(i,j);

 %--- cell velocity:
        if(istep <= 200 )
            vac_cell = 0.0;
        else
            vac_cell = vac(icell);
        end

 %---cell velocity vector:
                vnx = vac_cell;
                vny = 0.0;

 %--- fourth term:

                vnphi = vnx*phidx(i,j) + vny*phidy(i,j);
                term4x(i,j)=vnphi;

 %-- time integration:

                phi(i,j) = phi(i,j) + dtime*(gamma_cell*lap_phi(i,j) + term2 +term3 -vnphi);
           
            end %for j
        end %for i
        
        for i=1:Nx
            for j=1:Ny

 %-- for small deviations:
                if (phi(i,j) >= 0.9999)
                    phi(i,j) = 0.9999;
                end
                
                if(phi(i,j) < 0.000)
                    phi(i,j) = 0.0;
                end
 %---
                phis(i,j,icell) = phi(i,j);
            
            end %for j
        end %for i
     end %for icell
     
     else
         dtime = 1.0e-2;
     
     icell=2;
     gamma_cell = gamma;
         
         for i=1:Nx
             for j=1:Ny
                 phi(i,j) = phis(i,j,icell);
             end
         end
         
         %----
 % calculate the laplacian and gradient terms:
 %---
         for i=1:Nx
            for j=1:Ny
                jp=j+1;
                jm=j-1;

                ip=i+1;
                im=i-1;

                jp=j+1;
                jm=j-1;

                ip=i+1;
                im=i-1;

                if(im == 0)
                    im=Nx;
                end
                
                if(ip == (Nx+1))
                    ip=1;
                end
                
                if(jm == 0)
                    jm = Ny;
                end
                
                if(jp == (Ny+1))
                    jp=1;
                end

                hne=phi(ip,j);
                hnw=phi(im,j);
                hns=phi(i,jm);
                hnn=phi(i,jp);
                hnc=phi(i,j);

                lap_phi(i,j) =(hnw + hne + hns+ hnn -4.0*hnc)/(dx*dy);

    %--gradients of phi:
                phidx(i,j) = (phi(ip,j) - phi(im,j))/dx;
                phidy(i,j) = (phi(i,jp) - phi(i,jm))/dy;
            end %for j
        end %for i

 %----
 %--- volum integrations:
        vinteg = 0.0;
        vintegx = 0.0;
        vintegy = 0.0;

        for i=1:Nx
            for j=1:Ny
                vinteg = vinteg + phi(i,j)^2;
                sum_phi = 0.0;
                for jcell =1:ncell
                    if(icell ~= jcell)
                        sum_phi = sum_phi + phis(i,j,jcell)^2;
                    end
                end
                
                vintegx = vintegx + phi(i,j)* phidx(i,j)*sum_phi;
                vintegy = vintegy + phi(i,j)*phidy(i,j)*sum_phi;
            end
        end
        
        for i=1:Nx
            for j=1:Ny
 %-- Second term:
 %-- derivative of free energy
                [dfdphi]=free_energy_v1(i,j,icell,ncell,gamma,kappa,phi,phis);
                term2 = -const1*dfdphi;

 %---
 %---Third term
                term3 = -const2*(vinteg-pix*R^2)*phi(i,j);

 %--- cell velocity:
        if(istep <= 200 )
            vac_cell = 0.0;
        else
            vac_cell = vac(icell);
        end

 %---cell velocity vector:
                vnx = vac_cell;
                vny = 0.0;

 %--- fourth term:

                vnphi = vnx*phidx(i,j) + vny*phidy(i,j);
                term4x(i,j)=vnphi;

 %-- time integration:

                phi(i,j) = phi(i,j) + dtime*(gamma_cell*lap_phi(i,j) + term2 +term3 -vnphi);
           
            end %for j
        end %for i
        
        for i=1:Nx
            for j=1:Ny

 %-- for small deviations:
                if (phi(i,j) >= 0.9999)
                    phi(i,j) = 0.9999;
                end
                
                if(phi(i,j) < 0.000)
                    phi(i,j) = 0.0;
                end
 %---
                phis(i,j,icell) = phi(i,j);
            
            end %for j
        end %for i
 end %else

 %---- print results

    if((mod(istep,nprint) == 0) || (istep == 1) ) %mod() 取余数
        fprintf('done step: %5d\n',istep);

 %--- write vtk file:
        for i=1:Nx
            for j=1:Ny
                phi1(i,j) =0.0;
            end
        end
        
        for icell=1:ncell
            for i=1:Nx
                for j=1:Ny
                    phi1(i,j) = phi1(i,j) + (phis(i,j,icell)); %本来这里有^2
                end % for j
            end % for i
        end %icell
    
    write_vtk_grid_values_dlchen(Nx,Ny,dx,dy,istep,phi1);
    end %end if


 end %istep

 %--- calculate compute time:

 compute_time = etime(clock(),time0);
 fprintf('Compute Time: %10d\n',compute_time);