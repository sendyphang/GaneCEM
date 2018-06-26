function ZTLM2Dkernel_doConnection(grid)
    %% check polarisation 
    polarisation = grid.pol_type;
    
    %% connection    
    if strcmp(polarisation , 'Hz') 
        
        %% gather all component again because will get through swapping   
        if grid.gpu_yes_no
            grid.V2ir=gather(grid.V2ir);
            grid.V3ir=gather(grid.V3ir);
            grid.V4ir=gather(grid.V4ir);
            grid.V5ir=gather(grid.V5ir);
        end
        
        %% connection
        tempV5ir = grid.V5ir(:,1:end-1);
        grid.V5ir(:,1:end-1) = grid.V4ir(:,2:end);
        grid.V4ir(:,2:end) = tempV5ir;
    
        tempV3ir = grid.V3ir(1:end-1,:);
        grid.V3ir(1:end-1,:) = grid.V2ir(2:end,:);
        grid.V2ir(2:end,:) = tempV3ir;
        
        %% put to gpu again 
        if grid.gpu_yes_no
            
            grid.V2ir = gpuArray(grid.V2ir); grid.V3ir = gpuArray(grid.V3ir); 
            grid.V4ir = gpuArray(grid.V4ir); grid.V5ir = gpuArray(grid.V5ir);
            
        end
        
    else
        
        %% gather all component again because will get through swapping   
        if grid.gpu_yes_no
            grid.V8ir=gather(grid.V8ir);
            grid.V9ir=gather(grid.V9ir);
            grid.V10ir=gather(grid.V10ir);
            grid.V11ir=gather(grid.V11ir);
        end
        
        %% connection
        tempV11ir = grid.V11ir(:,1:end-1);
        grid.V11ir(:,1:end-1) = grid.V10ir(:,2:end);
        grid.V10ir(:,2:end) = tempV11ir;
    
        tempV9ir = grid.V9ir(1:end-1,:);
        grid.V9ir(1:end-1,:) = grid.V8ir(2:end,:);
        grid.V8ir(2:end,:) = tempV9ir;
        
        %% put to gpu again 
        if grid.gpu_yes_no
            
            grid.V8ir = gpuArray(grid.V8ir); grid.V9ir = gpuArray(grid.V9ir); 
            grid.V10ir = gpuArray(grid.V10ir); grid.V11ir = gpuArray(grid.V11ir);
            
        end
        
    end
end
