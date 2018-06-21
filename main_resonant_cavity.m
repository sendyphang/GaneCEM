
clc; clear; clf
reset(gpuDevice(1));
addpath('./src')


%% example problem from Microwave Engineering by David Pozar 4th ed
% page 116 example 3.1

load matLibrary.mat;

%% computational parameters 
length_x = 1.07e-2;
length_y = 0.43e-2;

d_l = 0.005e-2;

polarisation = 'Hz';%'Ez';  % Hz means the Ex, Ey, Hz != 0 and Hx, Hy, Vz = 0  while Ez polarization is the opposite   

NT = 40000; % number of timestep

% source parameters
source_centre_x = 0.2*length_x; source_centre_y = 0.2*length_y;
amplitude = 1; 


%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);


%% build geometry
makeGeometry.rectangle(grid,matLibrary({'Teflon'},:),0.5*length_x,0.5*length_y,length_x,length_y,true);

InOut.plotMaterial(grid,'n'); pause;% if you want to plot the material refractive index uncomment 

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

% monitor points position
mon1x = 0.8*length_x; mon1y = 0.8*length_y;


%% allocate memory for animation
% animation(1:grid.NY,1:grid.NX,1:NT) = 0;
% animation = gpuArray(animation);

tic
%% run TLM 
point_monitor(1,NT)=0; 
animation(grid.NY,grid.NX,round(NT/50))=0;
zzz = waitbar(0,'Please wait...');
for T = 1:NT
    waitbar(T / NT) ; %% update waitbar
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    
    %% source
    if T == 1
        source.dipole_heaviside(grid,amplitude,source_centre_x,source_centre_y);
    end
    
    
    %% scattering and connection
    ZTLM2Dkernel_doScatConn(grid)
    %% boundary 
    ZTLM2Dkernel_BoundaryHandling(grid, 'PEC', 'PEC' , 'PEC' , 'PEC' );
%     ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data 
    if (T>0 && T<NT)&&(mod(T,50)==0)
        if strcmp(polarisation , 'Hz')
            if ~grid.gpu_yes_no
                animation(:,:,T) = (grid.i_z);
                point_monitor(1,T)= InOut.monitor_point(grid,grid.i_z,mon1x,mon1y);
            else 
                animation(:,:,T) = gather(grid.i_z);
                point_monitor(1,T)= gather(InOut.monitor_point(grid,grid.i_z,mon1x,mon1y));
            end

        else strcmp(obj.pol_type , 'Ez')
            if ~grid.gpu_yes_no
                animation(:,:,T) = (grid.V_z);
                point_monitor(1,T)= InOut.monitor_point(grid,grid.V_z,mon1x,mon1y);
            else
                animation(:,:,T) = gather(grid.V_z);
                point_monitor(1,T)= gather(InOut.monitor_point(grid,grid.V_z,mon1x,mon1y));
            end
    %         point_monitor(1,T)= InOut.monitor_point(grid,grid.V_z,mon1x,mon1y);
        end
    end
    
end
close(zzz); %% close waitbar
toc

%% spectra
figure(2)
dt = d_l/(sqrt(2)*physC.c0);
freq = (0:length(point_monitor)-1)./(length(point_monitor)*dt);
plot(freq,abs(fft(point_monitor)));


%% if want animation uncomment 
%
figure(3)
for T = 1:NT
    surf((animation(:,:,T))); view(0,90);  axis('image'); shading flat; %shading('interp');
    caxis([-0.01 0.01])
    getframe(gcf); 
    title_string =  ['at T = ' (num2str((T-1)))];
        title(title_string)    
end
%}














