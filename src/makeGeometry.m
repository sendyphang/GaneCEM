classdef makeGeometry
    % this is a ststic class to make geometry
    
    properties
    end
    
    methods(Static)
        
        %% makeRectangle
        % logical value of snap_on_grid -> true if you want the rectangular
        % to be perfectly snapped on the grid, false otherwise
        function rectangle(mygrid,material,x0,y0,xspan,yspan,snap_on_grid)
            
            xi = (x0 - xspan/2.)/mygrid.dl ; yi = (y0 - yspan/2.)/mygrid.dl; % bottom left point
            
            Xspan = xspan/mygrid.dl; Yspan = yspan/mygrid.dl;
            
            %%
            if snap_on_grid 
                if ( abs(Xspan-round(Xspan))>1e-7 || abs(Yspan-round(Yspan))>1e-7 ) 
                    error('makeGeometry::no_nodes, not an integer');
                end
            end
            
            xi=round(xi); yi=round(yi);
            Xspan=round(Xspan); Yspan=round(Yspan);
            
            %%
            if xi == 0
                xi = 1;
            end
            
            if yi == 0
                yi = 1;
            end
            
            Xspan = Xspan - 1;
            Yspan = Yspan - 1;
                       
%             fprintf('\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Yspan : %i   yi : %i',Yspan,yi); pause;
            
            if ((xi+Xspan>mygrid.NX) || (yi+Yspan>mygrid.NY) )
                error('makeGeometry::too long')
            end
            
            %%
            % carefull that matlab has (0,0) on top left corner 
            % first dimension is row whilst second is coloumn
            for i = xi:xi+Xspan     % associate material ID to the matTag in grid object
                for j = yi:(yi+Yspan) %j = mygrid.NY+1-yi:-1:mygrid.NY+1-(yi+Yspan)
                    
                    mygrid.matTag(j,i) = material.ID;
                    mygrid.matCat(j,i) = material.Cat; 
                    mygrid.matSus0(j,i) = material.susceptibility_0; 
                    mygrid.matCond0(j,i) = material.conductivity_0; % real material conductivity 
                    
                    mygrid.matOmgPlasma(j,i) = material.omg_plasma;
                    mygrid.matDamping(j,i) = material.damping;
                    %... add later on
                    
                end
                
            end 
            
            %% to list the type material
            if material.Cat == 1
                mygrid.matListing(1) = true;
            elseif material.Cat == 2
                mygrid.matListing(material.Cat) = true;
            else
               %% if you want to have different material put here 
            end
            
            %%
            fprintf('//----------------------------------------------------------------------------------------------------');
            fprintf('\nmakeRectangle : OK!');
            fprintf('\n\tfor material : %s',material.Properties.RowNames{1});
            fprintf('\n\tMaterial Catagory : %i',material.Cat);
            fprintf('\n\tcentred at : ( %e m x %e m )',x0,y0);  
            fprintf('\n\tspan : ( %e m x %e m )\n\n',xspan,yspan);  
        end
        
        
        %% makeCircle
        function circle( mygrid, material, centre_x, centre_y, radius, angle_start, angle_end)
            
            x_0 = centre_x;
            y_0 = centre_y;
            
            thei = angle_start; 
            thee = angle_end;
            
            if thei < 0. 
                thei = 360. + thei;
            end
            if thei >= 360. 
                thei = thei - 360.;
            end
            if thee < 0. 
                thee = 360. + thee;
            end
            if thei >= 360.  
                thee = thee - 360.;
            end
            %%
            
            %%
            % carefull that matlab has (0,0) on top left corner 
            % first dimension is row whilst second is coloumn
            for i = 1:mygrid.NX     % associate material ID to the matTag in grid object
                
                xnorm = (((i)+0.5)*mygrid.dl) - x_0;
                
                for j = 1:mygrid.NY %j = mygrid.NY+1-yi:-1:mygrid.NY+1-(yi+Yspan)
                    
                    ynorm = (((j)+0.5)*mygrid.dl) - y_0;
                    rad = sqrt( (xnorm*xnorm) + (ynorm*ynorm) );
                    
                    the = mathFunc.getAngle(xnorm,ynorm);
                    
                    if rad <= radius

                        if thei < thee

                            if the >= thei && the <= thee
                                mygrid.matTag(j,i) = material.ID;
                                mygrid.matCat(j,i) = material.Cat; 
                                mygrid.matSus0(j,i) = material.susceptibility_0; 
                                mygrid.matCond0(j,i) = material.conductivity_0; % real material conductivity 
                                
                                mygrid.matOmgPlasma(j,i) = material.omg_plasma;
                                mygrid.matDamping(j,i) = material.damping;
                                %... add later on
                            end

                        end
                        if thei > thee

                            if the <= thee || the >= thei 
                                mygrid.matTag(j,i) = material.ID;
                                mygrid.matCat(j,i) = material.Cat; 
                                mygrid.matSus0(j,i) = material.susceptibility_0; 
                                mygrid.matCond0(j,i) = material.conductivity_0; % real material conductivity 
                                
                                mygrid.matOmgPlasma(j,i) = material.omg_plasma;
                                mygrid.matDamping(j,i) = material.damping;
                                %... add later on
                            end

                        end
                        if thei == thee 

                            mygrid.matTag(j,i) = material.ID;
                            mygrid.matCat(j,i) = material.Cat; 
                            mygrid.matSus0(j,i) = material.susceptibility_0; 
                            mygrid.matCond0(j,i) = material.conductivity_0; % real material conductivity 
                            
                            mygrid.matOmgPlasma(j,i) = material.omg_plasma;
                            mygrid.matDamping(j,i) = material.damping;
                                %... add later on

                        end
                    end
                                        
                end
                
            end 
            
            %% to list the type material
            if material.Cat == 1
                mygrid.matListing(1) = true;
            elseif material.Cat == 2
                mygrid.matListing(material.Cat) = true;
            else
               %% if you want to have different material put here 
            end
            
            %%
            fprintf('//----------------------------------------------------------------------------------------------------');
            fprintf('\nmakeCircle : OK!');
            fprintf('\n\tfor material : %s',material.Properties.RowNames{1});
            fprintf('\n\tMaterial Catagory : %i',material.Cat);
            fprintf('\n\tcentred at : ( %e m x %e m )',x_0,y_0);  
            fprintf('\n\tRadius : %e m ',radius);  
            fprintf('\n\tAngle : ( %f  -- %f  )\n\n',thei,thee); 
        end
        
         
        %% make .... 
        % if you want to make different kind of geometries (say, gratings?)
        % write here
        
        
    end
    
end
