classdef source
    %% source functions
    
    properties
    end
    
    methods(Static)
        
        %% heaviside function 
        function dipole_heaviside(grid,A,x0,y0)
            %% check computational dimension 
            NX = grid.NX; NY = grid.NY; polarisation = grid.pol_type;

%             X0 = round((x0-grid.dl/2.)/grid.dl); Y0 = round((y0-grid.dl/2.)/grid.dl);
            X0 = ceil((x0)/grid.dl); Y0 = ceil((y0)/grid.dl);
            
            if X0>NX || X0<0 || Y0>NY || Y0<0
                error('SOURCE::outside computational windwow');
            end
            
            if X0 ==0
                X0 = 1;
            end
            if Y0 ==0
                Y0 = 1;
            end
 
            %% normalized amplitude against the material at which it's excited
            if strcmp(polarisation , 'Hz') 
                A = A/sqrt(2.)/sqrt(grid.matSus0(Y0,X0)+1.);
            else
                A = sqrt(grid.matSus0(Y0,X0)+1.)*A/sqrt(2.);
            end
            
            %% do source
            if strcmp(polarisation , 'Hz') 
                grid.i_z(Y0,X0) = grid.i_z(Y0,X0) + A;
            else
                grid.V_z(Y0,X0) = grid.V_z(Y0,X0) + A;
            end
 
        end
            
        %% Gaussian pulse
        % grid to get access to the grid field
        % A = the amplitude
        % x0 is where the point source 
        % idT = current time step
        % tmax = when it reached the maximum value
        % r = truncation in ratio
        % h = the width from the middle to the trunctaion
        % stat = return 1 if on 0 if off
        function stat = dipole_gauss(grid,A,x0,y0,idT,tmax,r,h)
            
            %% check computational dimension 
            NX = grid.NX; NY = grid.NY; polarisation = grid.pol_type;
            X0 = ceil((x0)/grid.dl); Y0 = ceil((y0)/grid.dl);
%             X0 = round((x0-grid.dl/2.)/grid.dl); Y0 = round((y0-grid.dl/2.)/grid.dl);
            
            if X0>NX || X0<0 || Y0>NY || Y0<0
                error('SOURCE::outside computational windwow');
            end
            
            if X0 ==0
                X0 = 1;
            end
            if Y0 ==0
                Y0 = 1;
            end
            
            
            
            %% normalized amplitude against the material at which it's excited
            if strcmp(polarisation , 'Hz') 
                A = A/sqrt(2.)/sqrt(grid.matSus0(Y0,X0)+1.);
            else
                A = sqrt(grid.matSus0(Y0,X0)+1.)*A/sqrt(2.);
            end

            % Calculate the value of the source 
            gauss = A*exp((log(r)*(idT*grid.dt - tmax)^2)/(h*h));

            if (abs(gauss/A)<(0.25*r) && (idT*grid.dt - tmax)>0. )
                gauss = 0; %turn source off
            end

            %% do source
            if strcmp(polarisation , 'Hz') 
                grid.i_z(Y0,X0) = grid.i_z(Y0,X0) + gauss;
            else
                grid.V_z(Y0,X0) = grid.V_z(Y0,X0) + gauss;
            end
            
            %% return status
            if (abs(gauss/A)<(0.25*r) && (idT*grid.dt - tmax)>0. )
                stat =  0; %source off
            
            else
                stat = 1; %source on
            end
            
        end
        
        %% Gaussian wavepackage surface package
        % other parameters see gaussian pulse
        % f0 middle frequency
        function stat = dipole_gaussWavepacket(grid,A,x0,y0,idT,tmax,r,h,f0)  
                        
            %% check computational dimension 
            NX = grid.NX; NY = grid.NY; polarisation = grid.pol_type;

            X0 = ceil((x0-grid.dl/2.)/grid.dl); Y0 = round((y0-grid.dl/2.)/grid.dl);
            
            if X0>NX || X0<0 || Y0>NY || Y0<0
                error('SOURCE::outside computational windwow');
            end
            
            if X0 ==0
                X0 = 1;
            end
            if Y0 ==0
                Y0 = 1;
            end
            
            %% normalized amplitude against the material at which it's excited
            if strcmp(polarisation , 'Hz') 
                A = A/sqrt(2.)/sqrt(grid.matSus0(Y0,X0)+1.);
            else
                A = sqrt(grid.matSus0(Y0,X0)+1.)*A/sqrt(2.);
            end

            % Calculate the value of the source 
            gauss = A*exp((log(r)*(idT*grid.dt - tmax)^2)/(h*h));
            gaussWP = gauss*sin(2*pi*f0*idT*grid.dt);

            if (abs(gauss/A)<(0.25*r) && (idT*grid.dt - tmax)>0. )
                gauss = 0; %turn source off
            end

            %% do source
            if strcmp(polarisation , 'Hz') 
                grid.i_z(Y0,X0) = grid.i_z(Y0,X0) + gaussWP;
            else
                grid.V_z(Y0,X0) = grid.V_z(Y0,X0) + gaussWP;
            end
            
            %% return status
            if (abs(gauss/A)<(0.25*r) && (idT*grid.dt - tmax)>0. )
                stat =  0; %source off
            
            else
                stat = 1; %source on
            end
            
        end
        
        %% Continuous wave
        % A for amplitude
        % f0 middle frequency
        function dipole_CW(grid,A,x0,y0,idT,f0)   
            %% check computational dimension 
            NX = grid.NX; NY = grid.NY; polarisation = grid.pol_type;

            X0 = round((x0-grid.dl/2.)/grid.dl); Y0 = round((y0-grid.dl/2.)/grid.dl);
            
            if X0>NX || X0<0 || Y0>NY || Y0<0
                error('SOURCE::outside computational windwow');
            end
            
            if X0 ==0
                X0 = 1;
            end
            if Y0 ==0
                Y0 = 1;
            end
            
            %% normalized amplitude against the material at which it's excited
            if strcmp(polarisation , 'Hz') 
                A = A/sqrt(2.)/sqrt(grid.matSus0(Y0,X0)+1.);
            else
                A = sqrt(grid.matSus0(Y0,X0)+1.)*A/sqrt(2.);
            end

            % Calculate the value of the source 
            sine = A*sin(2*pi*f0*idT*grid.dt);

            %% do source
            if strcmp(polarisation , 'Hz') 
                grid.i_z(Y0,X0) = grid.i_z(Y0,X0) + sine;
            else
                grid.V_z(Y0,X0) = grid.V_z(Y0,X0) + sine;
            end

        end
        
        
        
        
    end
    
end

