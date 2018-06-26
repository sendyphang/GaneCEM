clc; clear;
addpath('./src')
reset(gpuDevice(1));

load matLibrary.mat;

%% plasmonic nanowire
radius = 40e-9;


%% computational parameters 
length_x = 6.*radius ;
length_y = 6.*radius ;

d_l = radius/60;

polarisation = 'Hz'; % for plasmonic support Hz polarisation   

NT = 10000; % number of timestep

% source parameters
source_centre_x = 0.5*length_x; source_centre_y = 0.5*length_y-1.1*radius;
amplitude = 1; 


%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);


%% build geometry
makeGeometry.circle(grid,matLibrary({'Gold_lossless'},:),0.5*length_x,0.5*length_y,radius,0,360);

% InOut.plotMaterial(grid,'n'); pause;% if you want to plot the material refractive index uncomment 

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

%% for write function for debug or demonstration 
[x,y] = meshgrid(((0:grid.NX-1)+0.5)*grid.dl,((0:grid.NY-1)+0.5)*grid.dl);
% x = reshape(x,[numel(x),1]);
% y = reshape(y,[numel(y),1]);

% monitor points position
mon1x = 0.5*length_x; mon1y = 0.5*length_y - 0.9*radius;

%% allocate memory for animation
% point_monitor(1,NT)=0; 
%     point_monitor = gpuArray(point_monitor);
animation(grid.NY,grid.NX,round(NT/10))=0;
    animation = gpuArray(animation);
id_animation =1;

tic
%% run TLM 
zzz = waitbar(0,'Please wait...');

for T = 1:NT
    waitbar(T / NT) ; %% update waitbar
    
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    
    %% source
%     for xcen = d_l:d_l:length_x
%         status = source.dipole_gaussWavepacket(grid,amplitude,xcen,source_centre_y,T,502*grid.dt,0.001,500*grid.dt,1336e12);
%     end
    status = source.dipole_gaussWavepacket(grid,amplitude,source_centre_x,source_centre_y,T,7/1543.8e12,0.001,5/1543.8e12,1543.8e12);
%     if T == 1
%         source.dipole_heaviside(grid,amplitude,source_centre_x,source_centre_y);
%     end
    
    
    %% scattering and connection
    ZTLM2Dkernel_doScatConn(grid)
    %% boundary 
    ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data 
    if (mod(T,10)==1)
        if strcmp(polarisation , 'Hz')
            if ~grid.gpu_yes_no
                animation(:,:,id_animation) = (grid.i_z);
            else 
                animation(:,:,id_animation) = gather(grid.i_z);
            end
        else strcmp(obj.pol_type , 'Ez')
            if ~grid.gpu_yes_no
                animation(:,:,id_animation) = (grid.V_z);
            else
                animation(:,:,id_animation) = gather(grid.V_z);
            end
        end
        id_animation = id_animation +1;
    end
    
    %{
    if strcmp(polarisation , 'Hz')
        point_monitor(1,T)= InOut.monitor_point(grid,grid.i_z,mon1x,mon1y);
    else strcmp(obj.pol_type , 'Ez')
        point_monitor(1,T)= InOut.monitor_point(grid,grid.V_z,mon1x,mon1y);
    end
    %}
    
end
close(zzz); %% close waitbar
toc

%% spectra
%{
figure(1)
hold on;
dt = d_l/(sqrt(2)*physC.c0);
freq = (0:length(point_monitor)-1)./(length(point_monitor)*dt);
plot(freq,abs(fft(point_monitor)));
%}

%% if want animation uncomment 
%
figure(3); clf
    the = linspace(0,2*pi,100);
    cx = 0.5*length_x+radius*cos(the);
    cy = 0.5*length_y+radius*sin(the);
[~,~,NTT] = size(animation);
for T = 1:NTT-2
    set(gcf,'color','w');
    plot3(cx*1e9,cy*1e9,cy.*0+12,'color','[0.93 0.69 0.13]'); hold on
    surf(x*1e9,y*1e9,(animation(:,:,T))); 
        view(0,90);  
        xlabel('x / nm'); ylabel('y / nm'); 
        title_string =  ['at T = ' (num2str((T-1)))];
        title(title_string)    
        axis('image'); 
        colormap jet;
    %     shading flat; 
        shading('interp');
        set(gca,'Ydir','normal')

        caxis([-0.01 0.01])
    M(T) = getframe(gcf); T = T +1; 
    
    clf
end














