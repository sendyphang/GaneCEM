
clc; clear; clf 
reset(gpuDevice(1));
addpath('./src')

load matLibrary.mat;

%% computational parameters 
d_l = 1;

length_x = 251*d_l;
length_y = 251*d_l;

polarisation = 'Ez'; % 'Hz'; % Hz means the Ex, Ey, Hz != 0 and Hx, Hy, Vz = 0  while Ez polarization is the opposite   

NT = 300; % number of timestep

% source parameters
source_centre_x = 0.5*length_x; source_centre_y = 0.5*length_y;
amplitude = 1; 


%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);


%% build geometry
makeGeometry.rectangle(grid,matLibrary({'Vacuum'},:),0.5*length_x,0.5*length_y,length_x,length_y,true);

% InOut.plotMaterial(grid,'n'); hold on; % if you want to plot the material refractive index uncomment 
% InOut.plotTLMnodes(grid,4000);

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

% % monitor points position
% mon1x = 0.8*length_x; mon1y = 0.8*length_y;

%% for write function for debug or demonstration 
[x,y] = meshgrid(((0:grid.NX-1)+0.5)*grid.dl,((0:grid.NY-1)+0.5)*grid.dl);
% x = reshape(x,[numel(x),1]);
% y = reshape(y,[numel(y),1]);

%% allocate memory for animation
animation(grid.NY,grid.NX,round(NT)/2)=0;
animation = gpuArray(animation);
% point_monitor(1,NT)=0; 

tic
%% run TLM 
zzz = waitbar(0,'Please wait...');
TT = 1; 
for T = 1:NT
    if  mod(T,5) == 1
        waitbar(T / NT) ; %% update waitbar
    end
        
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    
    %% source
    %{
    if T == 1
        grid.V8ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V9ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V10ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V11ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
%         source.dipole_heaviside(grid,amplitude,source_centre_x,source_centre_y);
    end
    %}
    
    source.dipole_gauss(grid,amplitude,source_centre_x,source_centre_y,T,10.*grid.dt,0.001,7*grid.dt);
    
    %% scattering and connection
    ZTLM2Dkernel_doScattering(grid);
        
    %% connection
    ZTLM2Dkernel_doConnection(grid);
    
    %% boundary 
%     ZTLM2Dkernel_BoundaryHandling(grid, 'PEC', 'PEC' , 'PEC' , 'PEC' );
    ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data 
    if mod(T,2)==1
        if strcmp(polarisation , 'Hz')
            if ~grid.gpu_yes_no
                animation(:,:,TT) = (grid.V_z);
                TT = TT +1;
            else 
                animation(:,:,TT) = gather(grid.V_z);
                TT = TT +1;
            end

        else 
            if ~grid.gpu_yes_no
                animation(:,:,TT) = (grid.V_z);
                TT = TT +1;
            else
                animation(:,:,TT) = gather(grid.V_z);
                TT = TT +1;
            end
        end
    end
    
end
close(zzz); %% close waitbar
toc

%% if want animation uncomment 
%
figure(3); clf
[~,~,NTT] = size(animation);
for T = 1:NTT-2
    set(gcf,'color','w');
    surf(x,y,(animation(:,:,T))); 
    view(0,90);  
    title_string =  ['at T = ' (num2str((T-1)))];
    title(title_string)    
    axis('image'); 
    colormap jet;
%     shading flat; 
    shading('interp');
    set(gca,'Ydir','normal')
    
    caxis([-0.1 0.1])
    M(T) = getframe(gcf); T = T +1; 
    
    clf
end
%}

%% 
writerObj = VideoWriter('animation.avi');
writerObj.FrameRate = 5; 
open(writerObj);
writeVideo(writerObj,M);
close(writerObj);












