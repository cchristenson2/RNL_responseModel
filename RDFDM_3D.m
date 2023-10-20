%{
Forward evaluation of the Reaction-Diffusion model in 3D for modeling treatment response
Author: Chase - 11/2/2021

dN/dt = grad*(D*grad(N)) + kp*N*(1-N/theta)

Inputs:
    - initial: starting point for forward evaulation (sy,sx,sz)
    - D: scalar diffusivity
    - kp: proliferation scalar or map (sy,sx,sz)
    - dx, dy, dt: FDM grid spacing variables
    - tspan: times to extract cell simulations
    - bcs: Boundary condition map for the domain; (sy,sx,sz,3)
    - theta: carrying capacity for cells

Outputs:
    - sim: Cell maps at times of tspan
%}

function Sim = RDFDM_3D(initial, D, kp, dx, dy, dz, dt, tspan, bcs, theta)
    [sy,sx,sz] = size(initial);
    t_idx = 1 + tspan/dt;
    nt = t_idx(end);
    
    %Check if kp is scalar
    if(numel(kp)==1)
        kp = kp.*ones(sy,sx,sz);
    end
    
    %Initialize FDM grid and plug in initial values
    N = zeros(sy,sx,sz,nt);
    N(:,:,:,1) = initial;
    
    thresh = theta*0.01; %Small threshold to determine if cells are present
    
    
    for t = 1:nt
        for z = 1:sz
            for y = 1:sy
                for x = 1:sx
                    if(bcs(y,x,z,1)~=2)
                        %Adaptive proliferation
                        %If tumor moves into location with kp=0, but is above a threshold
                        %Kp from nearest non-zero neighbors is added to the voxel
                        if(kp(y,x,z) == 0 && N(y,x,z) > thresh)
                            try
                                nn = [kp(y+1,x-1,z),kp(y+1,x,z),kp(y+1,x+1,z),kp(y,x+1,z),kp(y-1,x+1,z),kp(y-1,x,z),kp(y-1,x-1,z),kp(y,x-1,z)];
                                nn = nn(nn>0);
                                kp(y,x,z) = sum(nn,'all')/numel(nn);
                                if(isnan(kp(y,x,z)))
                                    kp(y,x,z)=0;
                                end
                            catch
                                %Default low proliferation rate
                                kp(y,x,z) = 1e-3;
                            end
                        end
            
                        %FDM on Y Boundaries
                        if(bcs(y,x,z,1)==0)
                            inv_y = D*(N(y+1,x,z,t)-2*N(y,x,z,t)+N(y-1,x,z,t))/(dy^2);
                        elseif(bcs(y,x,z,1)==1)
                            inv_y = D*(-2*N(y,x,z,t)+2*N(y-1,x,z,t))/(dy^2);
                        elseif(bcs(y,x,z,1)==-1)
                            inv_y = D*(-2*N(y,x,z,t)+2*N(y+1,x,z,t))/(dy^2);
                        end
            
                        %FDM on X Boundaries
                        if(bcs(y,x,z,2)==0)
                            inv_x = D*(N(y,x+1,z,t)-2*N(y,x,z,t)+N(y,x-1,z,t))/(dx^2);
                        elseif(bcs(y,x,z,2)==1)
                            inv_x = D*(-2*N(y,x,z,t)+2*N(y,x-1,z,t))/(dx^2);
                        elseif(bcs(y,x,z,2)==-1)
                            inv_x = D*(-2*N(y,x,z,t)+2*N(y,x+1,z,t))/(dx^2);
                        end
            
            
                        %FDM on Z Boundaries
                        if(bcs(y,x,z,3)==0)
                            inv_z = D*(N(y,x,z+1,t)-2*N(y,x,z,t)+N(y,x,z-1,t))/(dz^2);
                        elseif(bcs(y,x,z,3)==1)
                            inv_z = D*(-2*N(y,x,z,t)+2*N(y,x,z-1,t))/(dz^2);
                        elseif(bcs(y,x,z,3)==-1)
                            inv_z = D*(-2*N(y,x,z,t)+2*N(y,x,z+1,t))/(dz^2);
                        end
            
                        invasion = inv_y + inv_x + inv_z;
                        prolif    = N(y,x,z,t)*kp(y,x,z)*(1-(N(y,x,z,t)/theta));
                        

                        N(y,x,z,t+1) = N(y,x,z,t) + dt*(invasion + prolif);
                    end
                end
            end
        end
    end
    
    Sim = N(:,:,:,t_idx);
end