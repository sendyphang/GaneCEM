
clc; clear; 

addpath('./src')

load matLibrary.mat;

%% circular Silicon cavity
radius = 0.54e-6;


%% computational parameters 
length_x = 4.*radius ;
length_y = 4.*radius ;

d_l = radius/50;

polarisation = 'Ez';   

run_time = 0.5e-12; % real time 

% source parameters
source_centre_x = 0.5*length_x-0.8*radius; source_centre_y = 0.5*length_y;
amplitude = 1; 

f_op = 336.85e12; 

%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);

NT = run_time/grid.dt; % number of timestep 

%% build geometry
makeGeometry.circle(grid,matLibrary({'Silicon'},:),0.5*length_x,0.5*length_y,radius,0,360);

InOut.plotMaterial(grid,'n'); pause;% if you want to plot the material refractive index uncomment 

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

% monitor points position
mon1x = source_centre_x; mon1y = source_centre_y ;


%% allocate memory for animation
% animation(1:grid.NY,1:grid.NX,1:300) = 0;
id_animation =1;

tic
%% run TLM 
animation(grid.NY,grid.NX,round(NT/50))=0;
zzz = waitbar(0,'Please wait...');
for T = 1:NT
    waitbar(T / NT) ; %% update waitbar
    
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid); 
    
    %% source
    source_status = source.dipole_gaussWavepacket(grid,amplitude,source_centre_x,source_centre_y,T,10.1/f_op,0.001,10/f_op,f_op);

    %% scattering and connection
    ZTLM2Dkernel_doScatConn(grid);
    %% boundary 
    ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data
    %
        if (id_animation<1000)&&(mod(T,50)==0)
            if strcmp(polarisation , 'Hz')
                if ~grid.gpu_yes_no
                    animation(:,:,id_animation) = (grid.i_z);
                else 
                    animation(:,:,id_animation) = gather(grid.i_z);
                end
            else 
                if ~grid.gpu_yes_no
                    animation(:,:,id_animation) = (grid.V_z);
                else
                    animation(:,:,id_animation) = gather(grid.V_z);
                end
            end
            id_animation = id_animation +1;
        end
    %}
    
    
    if(~source_status)
        if strcmp(polarisation , 'Hz')
            point_monitor(1,T)= InOut.monitor_point(grid,grid.i_z,mon1x,mon1y);
        else 
            point_monitor(1,T)= InOut.monitor_point(grid,grid.V_z,mon1x,mon1y);
        end
    end
    
end
close(zzz); %% close waitbar
toc

%% spectra
figure(1)
hold on;
dt = d_l/(sqrt(2)*physC.c0);
freq = (0:length(point_monitor)-1)./(length(point_monitor)*dt);
plot(freq,abs(fft(point_monitor)));


%% if want animation uncomment 
%
figure(2)
[nr,nc,np] = size(animation);
for T = 1:np
    pcolor((animation(:,:,T))); view(0,90);  axis('image'); 
    title_string =  ['at T = ' (num2str((T-1)))];
        title(title_string)
        
%     shading flat; 
    shading('interp');
    colormap jet;
    caxis([-0.1 0.1]);
    getframe(gcf); 
        
end
%}














