function bcs = BuildBoundaries(brain)
    [sy,sx,sz] = size(brain);
    
    bcs = zeros(sy,sx,sz,3);
    for k = 1:sz
        for j = 1:sy
            for i = 1:sx
                
                boundary = [0,0,0]; %[y x z], exists on -1 to 1
                %Check y
                if(j == 1)
                    boundary(1)=-1;
                elseif(j==sy)
                    boundary(1)=1;
                end
                if(boundary(1)==0)
                    try %negative y-check
                        test = brain(j-1,i,k);
                        if(test==0)
                           boundary(1) = -1; 
                        end
                    catch
                        boundary(1) = -1; 
                    end
                end
                if(boundary(1)==0)
                    try %positive y-check
                        test = brain(j+1,i,k);
                        if(test==0)
                           boundary(1) = 1; 
                        end
                    catch
                        boundary(1) = 1; 
                    end
                end
                
                %Check x
                if(i == 1)
                    boundary(2)=-1;
                elseif(i==sx)
                    boundary(2)=1;
                end
                if(boundary(2)==0)
                    try %negative x-check
                        test = brain(j,i-1,k);
                        if(test==0)
                           boundary(2) = -1; 
                        end
                    catch
                        boundary(2) = -1; 
                    end
                end
                if(boundary(2)==0)
                    try %positive x-check
                        test = brain(j,i+1,k);
                        if(test==0)
                           boundary(2) = 1; 
                        end
                    catch
                        boundary(2) = 1; 
                    end
                end
                
                %Check z
                if(k == 1)
                    boundary(3)=-1;
                elseif(k==sz)
                    boundary(3)=1;
                end
                if(boundary(3)==0)
                    try %negative z-check
                        test = brain(j,i,k-1);
                        if(test==0)
                           boundary(3) = -1; 
                        end
                    catch
                        boundary(3) = -1; 
                    end
                end
                if(boundary(3)==0)
                    try %positive z-check
                        test = brain(j,i,k+1);
                        if(test==0)
                           boundary(3) = 1; 
                        end
                    catch
                        boundary(3) = 1; 
                    end
                end
                
                bcs(j,i,k,:) = boundary;
            end
        end
    end
    

end