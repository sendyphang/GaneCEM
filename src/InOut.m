classdef InOut
    %% an IO class function 
    % static function for IO
    
    properties
    end
    
    methods (Static)
        
        %% for plot refractive index distribution on computational space
        function plotMaterial(grid,what_to_plot)
            
            if strcmp(what_to_plot , 'n')
                [x,y] = meshgrid(((0:grid.NX-1)+0.5)*grid.dl,((0:grid.NY-1)+0.5)*grid.dl);
                x = reshape(x,[numel(x),1]);
                y = reshape(y,[numel(y),1]);
%                 pcolor(((0:grid.NX-1))*grid.dl,((0:grid.NY-1))*grid.dl,sqrt(grid.matSus0+1));
                imagesc(x,y,sqrt(grid.matSus0+1));
                xlabel('x / m'); ylabel('y / m'); zlabel('Refractive index');
                axis('image'); view(0,90)
                set(gca,'Ydir','normal')
            else
                %....
            end
            
        end
        
        %% monitor point 
        function [value] = monitor_point(grid,what_to_observe,x0,y0)
            X = round(x0/grid.dl); Y = round(y0/grid.dl);

            value = what_to_observe(Y,X);
            
        end
        
        %% line point 
        function [value] = monitor_line(grid,what_to_observe,x0,xe,y0,ye)
            X0 = round(x0/grid.dl); Y0 = round(y0/grid.dl);
            Xe = round(xe/grid.dl); Ye = round(ye/grid.dl);
            
            if Xe<X0
                dum = Xe;
                Xe = X0;
                X0 = dum;
            end
              
            if Ye<Y0
                dum = Ye;
                Ye = Y0;
                Y0 = dum;
            end
            
            if X0 <= 0
                X0 = 1;
            end
            
            if Y0 <= 0
                Y0 = 1;
            end
            
            if Xe >= grid.NX
                Xe = grid.NX;
            end
            
            if Ye >= grid.NY
                Ye = grid.NY;
            end
            
            if Ye-Y0>0 && Xe-X0>0
                error('InOut::Line Monitor Both size > 0');
            end
            
            
            value = what_to_observe(Y0:Ye,X0:Xe);  
        end
        
        %% plot nodes for debug or demonstration
        function plotTLMnodes(grid,sz)
            [x,y] = meshgrid(((0:grid.NX-1)+0.5)*grid.dl,((0:grid.NY-1)+0.5)*grid.dl);
            x = reshape(x,[numel(x),1]);
            y = reshape(y,[numel(y),1]);
            scatter(x,y,sz,'+');
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            xlabel('x / m'); ylabel('y / m'); 
            axis('image');
                        
        end
        
        %% write a grid value to plot figure for debug
        function writeText(grid,x,y,what_to_plot)
            text(x,y,string(what_to_plot));
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            drawnow;
            xlabel('x / m'); ylabel('y / m'); 
            axis('image'); view(0,90);
            set(gca,'Ydir','normal')
        end
        
        %% if you want to add more function write down here
        
        
    end
    
end

