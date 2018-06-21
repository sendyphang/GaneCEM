
%% to handle different boundary condition 
function ZTLM2Dkernel_BoundaryHandling(grid, west, east, south, north )
    
    %% check computational dimension 
    NX = grid.NX; NY = grid.NY; polarisation = grid.pol_type;
    
    %% locate the boundary material -- the refractive index
    n_west = sqrt(grid.matSus0(:,1) + 1.);
    n_east = sqrt(grid.matSus0(:,end) + 1.);
    n_south = sqrt(grid.matSus0(1,:) + 1.);
    n_north = sqrt(grid.matSus0(end,:) + 1.);
    
    
    %% List of boundary conditions 
    pec = -1.;

    %% For series nodes - Hz polarisation
    if strcmp(polarisation , 'Hz') 
        
        r4(1:NY,1)=0.; % For left --> V4i = r4 * V4r */
        r5(1:NY,1)=0; % For right --> V5i = r5 * V5r */
        r2(1,1:NX)=0; % For bottom --> V2i = r2 * V2r */
        r3(1,1:NX)=0; % For top --> V3i = r3 * V3r */
        
        if grid.gpu_yes_no
            r4 = gpuArray(r4); r5 = gpuArray(r5); 
            r2 = gpuArray(r2); r3 = gpuArray(r3);
        end

        %% If PEC 
        if strcmp(west , 'PEC') 
            r4(1:end,1) = pec;
        end
        if strcmp(east , 'PEC')
            r5(1:end,1) = pec;
        end
        if strcmp(south , 'PEC')
            r2(1,1:end) = pec;
        end
        if strcmp(north , 'PEC')
            r3(1,1:end) = pec;
        end
        
        
        %% if matched boundary condition - MBC
        if strcmp(west , 'MBC') 
            r4(:,1) =  (1. - n_west) ./(1. + n_west);
        end
        if strcmp(east , 'MBC')
            r5(:,1) =  (1. - n_east) ./(1. + n_east);
        end
        if strcmp(south , 'MBC')
            r2(1,:) =  (1. - n_south) ./(1. + n_south);
        end
        if strcmp(north , 'MBC')
            r3(1,:) =  (1. - n_north)  ./(1. + n_north);
        end
        
        %% if periodic boundary - will be added later
        
        
        %% Perform boundary condition 
        % Left side
        if ~strcmp(west , 'P_LR') % if not periodic        
            grid.V4ir(1:end,1) = r4(1:end,1) .* grid.V4ir(1:end,1);
        end
        % Right side
        if ~strcmp(east , 'P_LR') % if not periodic        
            grid.V5ir(1:end,end) = r5(1:end,1) .* grid.V5ir(1:end,end);
        end
        % bottom side
        if ~strcmp(south , 'P_BT') % if not periodic        
            grid.V2ir(1,1:end) = r2(1,1:end) .* grid.V2ir(1,1:end);
        end
        % top side
        if ~strcmp(north , 'P_BT') % if not periodic        
            grid.V3ir(end,1:end) = r3(1,1:end) .* grid.V3ir(end,1:end);
        end
        
    end
    
    %% For shunt nodes - Ez polarization
    if strcmp(polarisation , 'Ez') 
        
        r10(1:NY,1)=0.; % For left --> V10i = r10 * V10r */
        r11(1:NY,1)=0; % For right --> V11i = r11 * V11r */
        r8(1,1:NX)=0; % For bottom --> V8i = r8 * V8r */
        r9(1,1:NX)=0; % For top --> V9i = r9 * V9r */
        
        if grid.gpu_yes_no
            r10 = gpuArray(r10); r11 = gpuArray(r11); 
            r8 = gpuArray(r8); r9 = gpuArray(r9);
        end

        %% If PEC 
        if strcmp(west , 'PEC') 
            r10(1:end,1) = pec;
        end
        if strcmp(east , 'PEC')
            r11(1:end,1) = pec;
        end
        if strcmp(south , 'PEC')
            r8(1,1:end) = pec;
        end
        if strcmp(north , 'PEC')
            r9(1,1:end) = pec;
        end
        
        %% if matched boundary condition - MBC
        if strcmp(west , 'MBC') 
            r10(1:end,1) =  (1. - n_west) ./(1. + n_west);
        end
        if strcmp(east , 'MBC')
            r11(1:end,1) =  (1. - n_east) ./(1. + n_east);
        end
        if strcmp(south , 'MBC')
            r8(1,1:end) =   (1. - n_south)  ./(1. + n_south);
        end
        if strcmp(north , 'MBC')
            r9(1,1:end) =   (1. - n_north)   ./(1. + n_north);
        end
        %% if periodic boundary - will be added later
        
        
        %% Perform boundary condition 
        % Left side
        if ~strcmp(west , 'P_LR') % if not periodic        
            grid.V10ir(1:end,1) = r10(1:end,1) .* grid.V10ir(1:end,1);
        end
        % Right side
        if ~strcmp(east , 'P_LR') % if not periodic        
            grid.V11ir(1:end,end) = r11(1:end,1) .* grid.V11ir(1:end,end);
        end
        % bottom side
        if ~strcmp(south , 'P_BT') % if not periodic        
            grid.V8ir(1,1:end) = r8(1,1:end) .* grid.V8ir(1,1:end);
        end
        % top side
        if ~strcmp(north , 'P_BT') % if not periodic        
            grid.V9ir(end,1:end) = r9(1,1:end) .* grid.V9ir(end,1:end);
        end
        
    end
    
    

    

    
    
end