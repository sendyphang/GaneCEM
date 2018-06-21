classdef makeMesh < handle
    
    properties
        
        %% as usual
        Lx=0;  % length x 
        Ly=0;  % length y
        
        dl=0; dt=0; % obvious
        
        pol_type; % series nodes -> Hz polarisation whilst shunt node is -> Ez polarisation
        
        NX=0; NY=0; % number of mesh point
        
        %% gpu setting false by default
        gpu_yes_no = false;
        
        %% material tag for each node 1 by default
        matListing; % to list if there is specific kind of material 
        
        matTag; % material tag for each node 1 by default
        matCat;
        matSus0; 
        matCond0;
        matOmgPlasma;
        matDamping;
        
        
        %...  you can add later on for different material type
        
        %% for H polarised
        V_x; V_y; i_z;   
        V2ir; V3ir; V4ir; V5ir;% will be use as V4i and V4r respectively
   
        Sex; Sey;
        Aex1;Aey1;Aex2;Aey2; % for drude model
        
        %% for E polarised
        i_x; i_y; V_z;       
        V8ir; V9ir; V10ir; V11ir; % will be use as V4i and V4r respectively
        
        Sez;
        %
        
    end
     
    methods
        
        function obj=makeMesh(x, y, d_l, polarisation) %constructor
            obj.Lx = x; 
            obj.Ly = y; 
            
            obj.dl = d_l;
            obj.dt = d_l/(sqrt(2)*physC.c0);
            
            obj.pol_type = polarisation;
            %%
            obj.NX = round(x/d_l);
                if abs(round(obj.NX)-obj.NX)>1e-7
                    error('ERR:Discretisation not an integer');
                end
                obj.NX=round(obj.NX);
            
            obj.NY = round(y/d_l);
                if abs(round(obj.NY)-obj.NY)>1e-7
                    error('ERR:Discretisation not an integer');
                end
                obj.NY=round(obj.NY);
            
            %%
            % carefull that matlab has (0,0) on top left corner 
            % first dimension is row whilst second is coloumn
            if strcmp(polarisation , 'Hz')
                
                obj.V_x(1:obj.NY,1:obj.NX)=0; obj.V_y(1:obj.NY,1:obj.NX)=0; obj.i_z(1:obj.NY,1:obj.NX)=0;
                
                obj.V2ir(1:obj.NY,1:obj.NX)=0; obj.V3ir(1:obj.NY,1:obj.NX)=0;
                obj.V4ir(1:obj.NY,1:obj.NX)=0; obj.V5ir(1:obj.NY,1:obj.NX)=0;
                
                obj.Sex(1:obj.NY,1:obj.NX)=0; obj.Sey(1:obj.NY,1:obj.NX)=0;
                
                obj.Aex1(1:obj.NY,1:obj.NX)=0;obj.Aey1(1:obj.NY,1:obj.NX)=0;
                obj.Aex2(1:obj.NY,1:obj.NX)=0;obj.Aey2(1:obj.NY,1:obj.NX)=0;
                
            elseif strcmp(polarisation , 'Ez')
                
                obj.i_x(1:obj.NY,1:obj.NX)=0; obj.i_y(1:obj.NY,1:obj.NX)=0; obj.V_z(1:obj.NY,1:obj.NX)=0;
                
                obj.V8ir(1:obj.NY,1:obj.NX)=0; obj.V9ir(1:obj.NY,1:obj.NX)=0;
                obj.V10ir(1:obj.NY,1:obj.NX)=0; obj.V11ir(1:obj.NY,1:obj.NX)=0;
                
                obj.Sez(1:obj.NY,1:obj.NX)=0; 
                    
            else
                error('Polarisation must be either Hz or Ez'); 
            end
            
            %% material related
            obj.matListing(1)=true; % for type 1 dielectric by default is true
            
            obj.matTag(1:obj.NY,1:obj.NX)=1;
            obj.matCat(1:obj.NY,1:obj.NX)=1;
            obj.matSus0(1:obj.NY,1:obj.NX)=0;
            obj.matCond0(1:obj.NY,1:obj.NX)=0;
            
            obj.matListing(2)=false;  % for type 2 plasmonic by default is false
            obj.matOmgPlasma(1:obj.NY,1:obj.NX)=0;
            obj.matDamping(1:obj.NY,1:obj.NX)=0;
            
            %%
            fprintf('//----------------------------------------------------------------------------------------------------');
            fprintf('\nMesh generator : OK!');
            fprintf('\n\tComputational window : ( %e m x %e m ) ',x , y);
            fprintf('\n\tDiscretisation : %e m',d_l);
            fprintf('\n\tTime discretisation parameter : %e s',obj.dt);
            fprintf('\n\tTotal number of node points : ( %i x %i = %i ) \n\n',obj.NX,obj.NY,obj.NX*obj.NY);
            fprintf('\nRECOMMENDATION - for homogeneous structure');
            fprintf('\n\tMax simulation frequency : %e Hz or min free-space wavelength : %e m',0.05/obj.dt,20*obj.dl);
            fprintf('\n\tNyquist criterion : %e Hz \n\n',0.5/obj.dt);

        end
        
        
        function GPU_parallelisation(obj)
            
            obj.gpu_yes_no = true;

            %% make array parallel
            obj.matTag=gpuArray(obj.matTag);
            obj.matCat=gpuArray(obj.matCat);
            obj.matSus0=gpuArray(obj.matSus0);
            obj.matCond0=gpuArray(obj.matCond0);
            
            obj.matOmgPlasma=gpuArray(obj.matOmgPlasma);
            obj.matDamping=gpuArray(obj.matDamping);
            
            if strcmp(obj.pol_type , 'Hz')
                
                obj.V_x = gpuArray(obj.V_x); obj.V_y = gpuArray(obj.V_y); obj.i_z = gpuArray(obj.i_z);
                
                obj.V2ir = gpuArray(obj.V2ir); obj.V3ir = gpuArray(obj.V3ir); 
                obj.V4ir = gpuArray(obj.V4ir); obj.V5ir = gpuArray(obj.V5ir);
                 
                obj.Sex  = gpuArray(obj.Sex); obj.Sey = gpuArray(obj.Sey);
                
                obj.Aex1  = gpuArray(obj.Aex1); obj.Aey1 = gpuArray(obj.Aey1);
                obj.Aex2  = gpuArray(obj.Aex2); obj.Aey2 = gpuArray(obj.Aey2);
                
            elseif strcmp(obj.pol_type , 'Ez')
                
                obj.i_x = gpuArray(obj.i_x); obj.i_y = gpuArray(obj.i_y); obj.V_z = gpuArray(obj.V_z);  
                
                obj.V8ir = gpuArray(obj.V8ir); obj.V9ir = gpuArray(obj.V9ir);
                obj.V10ir = gpuArray(obj.V10ir); obj.V11ir = gpuArray(obj.V11ir);
                 
                obj.Sez = gpuArray(obj.Sez);
                    
            else
                error('Polarisation must be either Hz or Ez'); 
            end

            %%
            fprintf('//----------------------------------------------------------------------------------------------------');
            fprintf('\nGrid is now in the GPU... \n\n');
            
        end
        
        
        
    end
    
    
end