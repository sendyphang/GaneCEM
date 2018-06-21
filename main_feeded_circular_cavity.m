
clc; clear;

addpath('./src')

reset(gpuDevice(1));

load matLibrary.mat;

%% circular Silicon cavity
radius = 0.54e-6;

%% slab waveguide 
width = 0.15e-6;
gap = 0.02e-6;

%% computational parameters 
length_x = 4.*radius + 2.*width;
length_y = 4.*radius ;

d_l = width/10;

polarisation = 'Ez'; % for plasmonic support Hz polarisation   

run_time = 0.5e-12; % real time 

% source parameters
source_centre_x = 0.5*length_x-radius-gap-0.5*width; source_centre_y = 2*d_l;
amplitude = 1; 

f_op = 336.85e12; 

%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);

NT = run_time/grid.dt; % number of timestep 

%% build geometry
makeGeometry.rectangle(grid,matLibrary({'Silicon'},:),0.5*length_x-radius-gap-0.5*width,0.5*length_y,width,length_y,true);
makeGeometry.circle(grid,matLibrary({'Silicon'},:),0.5*length_x,0.5*length_y,radius,0,360);
makeGeometry.rectangle(grid,matLibrary({'Silicon'},:),0.5*length_x+radius+gap+0.5*width,0.5*length_y,width,length_y,true);

InOut.plotMaterial(grid,'n'); % if you want to plot the material refractive index uncomment 

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

% monitor points position
mon1x = 0.5*length_x-radius-gap-0.5*width; mon1y = 10*d_l ;
mon2x = 0.5*length_x-radius-gap-0.5*width; mon2y = length_y-10*d_l ;
mon3x = 0.5*length_x+radius+gap+0.5*width; mon3y = 10*d_l ;
mon4x = 0.5*length_x+radius+gap+0.5*width; mon4y = length_y-10*d_l ;


%% allocate memory for animation
% animation(1:grid.NY,1:grid.NX,1:300) = 0;
id_animation =1;
source_status=false;

tic
%% run TLM 
zzz = waitbar(0,'Please wait...');
point_monitor1(1,NT)=0; 
point_monitor2(1,NT)=0; 
point_monitor3(1,NT)=0; 
point_monitor4(1,NT)=0; 
animation(grid.NY,grid.NX,round(NT/50))=0;
for T = 1:NT
    waitbar(T / NT) ; %% update waitbar
    
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    
    %% source
    source_status = source.dipole_gaussWavepacket(grid,amplitude,source_centre_x,source_centre_y,T,5.1/f_op,0.001,5/f_op,f_op);
%     source.dipole_CW(grid,amplitude,source_centre_x,source_centre_y,T,f_op) 
    
    %% scattering and connection
    ZTLM2Dkernel_doScatConn(grid)
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
            point_monitor1(1,T)= InOut.monitor_point(grid,grid.i_z,mon1x,mon1y);
            point_monitor2(1,T)= InOut.monitor_point(grid,grid.i_z,mon2x,mon2y);
            point_monitor3(1,T)= InOut.monitor_point(grid,grid.i_z,mon3x,mon3y);
            point_monitor4(1,T)= InOut.monitor_point(grid,grid.i_z,mon4x,mon4y);
        else 
            point_monitor1(1,T)= InOut.monitor_point(grid,grid.V_z,mon1x,mon1y);
            point_monitor2(1,T)= InOut.monitor_point(grid,grid.V_z,mon2x,mon2y);
            point_monitor3(1,T)= InOut.monitor_point(grid,grid.V_z,mon3x,mon3y);
            point_monitor4(1,T)= InOut.monitor_point(grid,grid.V_z,mon4x,mon4y);
        end
    end
    
end
close(zzz); %% close waitbar
toc

%% spectra
figure(2)
hold on;
dt = d_l/(sqrt(2)*physC.c0);
freq = (0:length(point_monitor1)-1)./(length(point_monitor1)*dt); hold on
plot(freq,abs(fft(point_monitor1)));
plot(freq,abs(fft(point_monitor2)));
plot(freq,abs(fft(point_monitor3)));
plot(freq,abs(fft(point_monitor4)));

%% if want animation uncomment 
figure(3)
for T = 1:id_animation
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














