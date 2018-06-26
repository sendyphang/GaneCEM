
clc; clear; clf 
reset(gpuDevice(1));
addpath('./src')



%% example problem from Microwave Engineering by David Pozar 4th ed
% page 116 example 3.1

load matLibrary.mat;

%% computational parameters 
d_l = 1;

length_x = 5*d_l;
length_y = 5*d_l;

polarisation = 'Ez'; % 'Hz'; % Hz means the Ex, Ey, Hz != 0 and Hx, Hy, Vz = 0  while Ez polarization is the opposite   

NT = 10; % number of timestep

% source parameters
source_centre_x = 0.5*length_x; source_centre_y = 0.5*length_y;
amplitude = sqrt(2); 


%% make grid
grid = makeMesh(length_x,length_y,d_l,polarisation);


%% build geometry
makeGeometry.rectangle(grid,matLibrary({'Vacuum'},:),0.5*length_x,0.5*length_y,length_x,length_y,true);

InOut.plotMaterial(grid,'n'); hold on; % if you want to plot the material refractive index uncomment 
InOut.plotTLMnodes(grid,4000);

%% send the mesh to GPU - comment if run in CPU
% grid.GPU_parallelisation(); 

% monitor points position
mon1x = 0.8*length_x; mon1y = 0.8*length_y;

%% for write function for debug or demonstration 
[x,y] = meshgrid(((0:grid.NX-1)+0.5)*grid.dl,((0:grid.NY-1)+0.5)*grid.dl);
x = reshape(x,[numel(x),1]);
y = reshape(y,[numel(y),1]);

%% allocate memory for animation
animation(grid.NY,grid.NX,round(NT))=0;
% animation = gpuArray(animation);
point_monitor(1,NT)=0; 

tic
%% run TLM 
% zzz = waitbar(0,'Please wait...');
ii = 1; 
for T = 1:NT
%     waitbar(T / NT) ; %% update waitbar
    set(gcf,'color','w');
    % 
    if(T>1)
        clf
        InOut.plotTLMnodes(grid,4000); hold on
        for rr=1:numel(x)
            if(abs(grid.V8ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)-0.4*grid.dl,string(grid.V8ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.1*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)-0.45*grid.dl 0];
                mArrow3(p2,p1,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V9ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)+0.4*grid.dl,string(grid.V9ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)+0.05*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)+0.4*grid.dl 0];
                mArrow3(p2,p1,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V10ir(rr))>1e-10)
                text(x(rr)-0.4*grid.dl,y(rr)+0.07*grid.dl,string(grid.V10ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)-0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p2,p1,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V11ir(rr))>1e-10)
                text(x(rr)+0.1*grid.dl,y(rr)+0.07*grid.dl,string(grid.V11ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)+0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)+0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p2,p1,'color','grey','linewidth',1,'stemWidth',0.01)
            end
        end
        
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            title(['Incomming at T = ' (num2str((T-1)))],'FontSize',13) 
            drawnow;
            xlabel('x / m'); ylabel('y / m'); 
            axis('image'); view(0,90);
            set(gca,'Ydir','normal')
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
%             saveas(gcf,['in at T = ' (num2str((T-1))) ' .png'])
            M(ii) = getframe(gcf); ii = ii +1;
            
    end
    %}
    
    %% calc field
    ZTLM2Dkernel_doFieldCalc(grid);
    %{ 
    clf
    InOut.plotTLMnodes(grid,4000);
    text(x+0.02*grid.dl,y+0.07*grid.dl,string(grid.V_z),'Color','red','FontSize',10);
        axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
        title(['at T = ' (num2str((T-1)))],'FontSize',13) 
        drawnow;
        xlabel('x / m'); ylabel('y / m'); 
        axis('image'); view(0,90);
        set(gca,'Ydir','normal')
        axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
    %}
    
    %% scattering and connection
    ZTLM2Dkernel_doScattering(grid)
    
    % 
    if (T>1)
        clf
        InOut.plotTLMnodes(grid,4000); hold on
        for rr=1:numel(x)
            if(abs(grid.V8ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)-0.4*grid.dl,string(grid.V8ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.1*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)-0.45*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V9ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)+0.4*grid.dl,string(grid.V9ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)+0.05*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)+0.4*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V10ir(rr))>1e-10)
                text(x(rr)-0.4*grid.dl,y(rr)+0.07*grid.dl,string(grid.V10ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)-0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V11ir(rr))>1e-10)
                text(x(rr)+0.1*grid.dl,y(rr)+0.07*grid.dl,string(grid.V11ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)+0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)+0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
        end
        
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            title(['Scatter at T = ' (num2str((T-1)))],'FontSize',13) 
            drawnow;
            xlabel('x / m'); ylabel('y / m'); 
            axis('image'); view(0,90);
            set(gca,'Ydir','normal')
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
%             saveas(gcf,['out at T = ' (num2str((T-1))) ' .png'])
            M(ii) = getframe(gcf); ii = ii +1; 
    end
    %}
    
    
    
    %%
    %% source
    if T == 1
        grid.V8ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V9ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V10ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
        grid.V11ir(ceil(source_centre_y/grid.dl),ceil(source_centre_x/grid.dl)) = 1;
%         source.dipole_heaviside(grid,amplitude,source_centre_x,source_centre_y);
        % 
        clf
        InOut.plotTLMnodes(grid,4000); hold on
        for rr=1:numel(x)
            if(abs(grid.V8ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)-0.4*grid.dl,string(grid.V8ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.1*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)-0.45*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V9ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)+0.4*grid.dl,string(grid.V9ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)+0.05*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)+0.4*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V10ir(rr))>1e-10)
                text(x(rr)-0.4*grid.dl,y(rr)+0.07*grid.dl,string(grid.V10ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)-0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V11ir(rr))>1e-10)
                text(x(rr)+0.1*grid.dl,y(rr)+0.07*grid.dl,string(grid.V11ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)+0.05*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)+0.4*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p1,p2,'color','grey','linewidth',1,'stemWidth',0.01)
            end
        end
        
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            title(['Initial condition at T = ' (num2str((T-1)))],'FontSize',13) 
            drawnow;
            xlabel('x / m'); ylabel('y / m'); 
            axis('image'); view(0,90);
            set(gca,'Ydir','normal')
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
%             saveas(gcf,['source at T = ' (num2str((T-1))) ' .png']);
            M(ii) = getframe(gcf); ii = ii +1;
        %}
        
    end
    
    %%
    
    ZTLM2Dkernel_doConnection(grid);
    % 
        clf
        InOut.plotTLMnodes(grid,4000); hold on
%         text(x+0.02*grid.dl,y-0.4*grid.dl,string(grid.V8ir),'Color','k','FontSize',8);
%         text(x+0.02*grid.dl,y+0.4*grid.dl,string(grid.V9ir),'Color','k','FontSize',8);
%         text(x-0.4*grid.dl,y+0.07*grid.dl,string(grid.V10ir),'Color','k','FontSize',8);
%         text(x+0.1*grid.dl,y+0.07*grid.dl,string(grid.V11ir),'Color','k','FontSize',8);
%         
%         clf
        InOut.plotTLMnodes(grid,4000); hold on
        for rr=1:numel(x)
            if(abs(grid.V8ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)-0.4*grid.dl,string(grid.V8ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)-0.4*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)-0.7*grid.dl 0];
                mArrow3(p2,p1,'color','red','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V9ir(rr))>1e-10)
                text(x(rr)+0.02*grid.dl,y(rr)+0.4*grid.dl,string(grid.V9ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.05*grid.dl y(rr)+0.4*grid.dl 0];
                p2 = [x(rr)-0.05*grid.dl y(rr)+0.7*grid.dl 0];
                mArrow3(p2,p1,'color','red','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V10ir(rr))>1e-10)
                text(x(rr)-0.4*grid.dl,y(rr)+0.07*grid.dl,string(grid.V10ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)-0.4*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)-0.7*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p2,p1,'color','red','linewidth',1,'stemWidth',0.01)
            end
            
            if(abs(grid.V11ir(rr))>1e-10)
                text(x(rr)+0.1*grid.dl,y(rr)+0.07*grid.dl,string(grid.V11ir(rr)),'Color','k','FontSize',8);
                p1 = [x(rr)+0.4*grid.dl y(rr)-0.05*grid.dl 0];
                p2 = [x(rr)+0.7*grid.dl y(rr)-0.05*grid.dl 0];
                mArrow3(p2,p1,'color','red','linewidth',1,'stemWidth',0.01)
            end
        end
        
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
            title(['Connection at T = ' (num2str((T-1)))],'FontSize',13) 
            drawnow;
            xlabel('x / m'); ylabel('y / m'); 
            axis('image'); view(0,90);
            set(gca,'Ydir','normal')
            axis([0 grid.NY*grid.dl 0 grid.NX*grid.dl])
%             saveas(gcf,['conn at T = ' (num2str((T-1))) ' .png'])
            M(ii) = getframe(gcf); ii = ii +1;
    %}
        
    
    %% boundary 
%     ZTLM2Dkernel_BoundaryHandling(grid, 'PEC', 'PEC' , 'PEC' , 'PEC' );
    ZTLM2Dkernel_BoundaryHandling(grid, 'MBC', 'MBC' , 'MBC' , 'MBC' );
    
    %% saving data 
    if (T>0 && T<NT)
        if strcmp(polarisation , 'Hz')
            if ~grid.gpu_yes_no
                animation(:,:,T) = (grid.i_z);
                point_monitor(1,T)= InOut.monitor_point(grid,grid.i_z,mon1x,mon1y);
            else 
                animation(:,:,T) = gather(grid.i_z);
                point_monitor(1,T)= gather(InOut.monitor_point(grid,grid.i_z,mon1x,mon1y));
            end

        else strcmp(polarisation , 'Ez')
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
% close(zzz); %% close waitbar
toc

%% spectra
% figure(2)
% dt = d_l/(sqrt(2)*physC.c0);
% freq = (0:length(point_monitor)-1)./(length(point_monitor)*dt);
% plot(freq,abs(fft(point_monitor)));


%% if want animation uncomment 
%
figure(3)
for T = 1:NT
    
    imagesc(x,y,(animation(:,:,T))); view(0,90);  axis('image'); shading flat; %shading('interp');
    set(gca,'Ydir','normal')
    
    caxis([-2 2])
    getframe(gcf); 
    title_string =  ['at T = ' (num2str((T-1)))];
    
        title(title_string)    
        clf
end
%}

%% 
writerObj = VideoWriter('animation.avi');
writerObj.FrameRate = 1; 
open(writerObj);
writeVideo(writerObj,M);
close(writerObj);












