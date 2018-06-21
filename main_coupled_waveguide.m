
clc; clear

addpath('./src')
reset(gpuDevice(1));

load matLibrary.mat;


%% slab waveguide 
length = 10e-6; 
width = 0.15e-6;
gap = 0.02e-6;

pad = 1.5*width;
%% computational parameters 
length_x = length;
length_y = pad+width+gap+width+pad ;

d_l = width/10;

polarisation = 'Ez'; % for plasmonic support Hz polarisation   

run_time = 0.5e-12; % real time 

% structure 
centre_x_1 = 0.5*length_x ; centre_y_1 = pad + 0.5*width;
centre_x_2 = 0.5*length_x ; centre_y_2 = centre_y_1 + gap + width;

% source parameters
source_centre_x = 2*d_l; source_centre_y = centre_y_1;
amplitude = 1; 

f_op = 336.85e12; 

%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);

NT = round(run_time/grid.dt); % number of timestep 

%% build geometry
makeGeometry.rectangle(grid,matLibrary({'Silicon'},:),centre_x_1,centre_y_1,length,width,true);
makeGeometry.rectangle(grid,matLibrary({'Silicon'},:),centre_x_2,centre_y_2,length,width,true);

%InOut.plotMaterial(grid,'n'); pause;% if you want to plot the material refractive index uncomment 

%% send the mesh to GPU - comment if run in CPU
grid.GPU_parallelisation(); 

% monitor points position
monx0_l1 = 0; monxe_l1 = length_x; mony0_l1 = centre_y_1 ; monye_l1 = centre_y_1 ;
monx0_l2 = 0; monxe_l2 = length_x; mony0_l2 = centre_y_2 ; monye_l2 = centre_y_2 ;

%% allocate memory for animation
% animation(1:grid.NY,1:grid.NX,1:300) = 0;
id_animation =1;
source_status=false;

tic
%% run TLM 
animation(grid.NY,grid.NX,round(NT/4))=0;
zzz = waitbar(0,'Please wait...');
for T = 1:NT
    waitbar(T / NT) ; %% update waitbar
    
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    
    %% source
    %source_status = source.dipole_gaussWavepacket(grid,amplitude,source_centre_x,source_centre_y,T,5.1/f_op,0.001,5/f_op,f_op);
     source.dipole_CW(grid,amplitude,source_centre_x,source_centre_y,T,f_op) 
    
    %% scattering and connection
    ZTLM2Dkernel_doScatConn(grid)
    %% boundary 
    ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data
    %
        if (id_animation<3000)&&(mod(T,4)==0)
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
    
    if (T==NT-1)
        line_monitor1(1,:) = InOut.monitor_line(grid,grid.V_z,monx0_l1,monxe_l1,mony0_l1,mony0_l1);
        line_monitor2(1,:) = InOut.monitor_line(grid,grid.V_z,monx0_l2,monxe_l2,mony0_l2,mony0_l2);
    end
    
end
close(zzz); %% close waitbar
toc

%% spectra
figure(1)
hold on;
[nR, nC] = size(line_monitor1);
beta = 2*pi.*((0:nC-1)./(nC*d_l));
plot(beta,abs(fft(line_monitor1)));
plot(beta,abs(fft(line_monitor2)));

%% if want animation uncomment 
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














