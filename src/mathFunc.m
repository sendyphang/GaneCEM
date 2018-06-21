

classdef mathFunc
    
    properties
    end
    
    methods(Static)
        
        function angle = getAngle(x,y)
            %% calculate angle
            % given double x
            % given double y
            % return angle in degree from 0 - 360
            
            %%
            theta = atan(y/x);
            
            theta = theta*180./pi;

            if x < 0. 
                theta = theta + 180.;
            else
                if y < 0.
                    theta = theta + 360.;
                else
                    theta = theta;
                end
            end
            
            angle = theta;     
        end
        
        %% if you want to add more function add here!
        
        
        
    end
    
    
end