function ZTLM2Dkernel_doFieldCalc(grid)
    
    %% check polarisation 
    polarisation = grid.pol_type;
    
    
    %% do field calculation
    %% I have set material catagory 1 as simple dielectric which is frequency and intensity independent
    if (grid.matListing(1))
        X1 = find(grid.matCat==1);  % find which node is material catagory 1 % BE AWARE, of indexing in matlab see find function for 2d matrix in matlab
    end
    
    if numel(X1)>0
        if strcmp(polarisation , 'Hz')  
            [grid.V_x(X1),...
             grid.V_y(X1),...
             grid.i_z(X1),...
             grid.Sex(X1),...
             grid.Sey(X1)] = arrayfun(@dFC_Hz_SimpleDielectric,...
                                        grid.matSus0(X1),...  % suscepetibility of dielectric
                                        grid.matCond0(X1) * grid.dl * physC.NU0 / sqrt(2.),...  % normalised conductivity
                                        grid.V2ir(X1),...
                                        grid.V3ir(X1),...
                                        grid.V4ir(X1),...
                                        grid.V5ir(X1),...
                                        grid.Sex(X1),...
                                        grid.Sey(X1));

        else
            [grid.i_x(X1),...
             grid.i_y(X1),...
             grid.V_z(X1),...
             grid.Sez(X1)] = arrayfun(@dFC_Ez_SimpleDielectric,...
                                        grid.matSus0(X1),...  % suscepetibility of dielectric
                                        grid.matCond0(X1) * grid.dl * physC.NU0 * sqrt(2.),...  % normalised conductivity
                                        grid.V8ir(X1),...
                                        grid.V9ir(X1),...
                                        grid.V10ir(X1),...
                                        grid.V11ir(X1),...
                                        grid.Sez(X1));
        end
    end
    
    
    
    %% I have set material catagory 2 as Drude model which is suitable for plasmonics
    X2=[]; 
    if (grid.matListing(2))
        X2 = find(grid.matCat==2);  % find which node is material catagory 1 % BE AWARE, of indexing in matlab see find function for 2d matrix in matlab
    end  
    
    if (numel(X2)>0)
        if strcmp(polarisation , 'Hz')  
            [grid.V_x(X2),...
             grid.V_y(X2),...
             grid.i_z(X2),...
             grid.Sex(X2),...
             grid.Sey(X2),...
             grid.Aex1(X2),...
             grid.Aey1(X2),...
             grid.Aex2(X2),...
             grid.Aey2(X2)] = arrayfun(@dFC_Hz_Drude,...
                                        grid.matOmgPlasma(X2)*grid.dt,...  % normalised plasma frequency in radial
                                        grid.matDamping(X2) *grid.dt,...  % normalised dumping in rad
                                        grid.V2ir(X2),...
                                        grid.V3ir(X2),...
                                        grid.V4ir(X2),...
                                        grid.V5ir(X2),...
                                        grid.Sex(X2),...
                                        grid.Sey(X2),...
                                        grid.Aex1(X2),...
                                        grid.Aey1(X2),...
                                        grid.Aex2(X2),...
                                        grid.Aey2(X2));
        end
    end
    
    
    
    %% if you want to model different material set the function header here
    % ...
    
end

%% function implementation to calculate the voltage scattering in case of simple dielectric - SERIES nodes
function [vx,vy,iz,sex,sey] = dFC_Hz_SimpleDielectric(chi_e,...
                                                      g_e,...
                                                      V2ir,...
                                                      V3ir,...
                                                      V4ir,...
                                                      V5ir,...
                                                      Sex,...
                                                      Sey)
    %% Calculate K_0
    K_0 = g_e + 2.*(1.+chi_e);
    
%     chi_e
% %     fprintf('xxxxxxxxxxxxxxxxxxx chi_e : %f' , chi_e);
%     pause;
                                    
    %% Calculate Vx_i etc - here V_2ir etc act as the incident
    Vxi = V2ir + V3ir; % Vxi = v2_i + v3_i 
    Vyi = V4ir + V5ir; % Vyi = v4_i + v5_i 
    izi = V2ir + V5ir - V3ir - V4ir ; % izi = v2_i - v3_i - v4_i + V5_i 
                                    
    %% Calculate V_x , V_y , i_z 
    vx = 2.*Vxi/K_0 + Sex/K_0; % V_x = 2*Vxi/K0 +z^-1(Sex/K0) 
    vy = 2.*Vyi/K_0 + Sey/K_0; % V_y = 2*Vyi/K0 +z^-1(Sey/K0) 
    iz = izi/-2.; % i_z = izi/-2 

    %% return sex , sey 
    sex = 2.*Vxi - vx*(2.*(1.-chi_e) + g_e); % Sex = 2Vxi + Vx(2(1-chi_e) + g_e) 
    sey = 2.*Vyi - vy*(2.*(1.-chi_e) + g_e); % Sey = 2Vyi + Vy(2(1-chi_e) + g_e)

end

%% function implementation to calculate the voltage scattering in case of drude model - Shunt nodes
function [vx,vy,iz,sex,sey,aex1,aey1,aex2,aey2] = dFC_Hz_Drude(omg_p_normalised,...
                                                               damping_normalised,...
                                                               V2ir,...
                                                               V3ir,...
                                                               V4ir,...
                                                               V5ir,...
                                                               Sex,...
                                                               Sey,...
                                                               Aex1,...
                                                               Aey1,...
                                                               Aex2,...
                                                               Aey2)
    %% some constants
    Ke1 = omg_p_normalised*omg_p_normalised;
    Ke2 = 4 + 2*damping_normalised;
    Ke3 = 4 - 2*damping_normalised;
    
    chi_e0 = Ke1/Ke2;
    chi_e1 = 0.;
    K_0 = 2*(1+chi_e0);
    
    b0 = (chi_e0)*(1+8/Ke2);
    b1 = (chi_e0)*(-1-Ke3/Ke2);
    b2 = -chi_e0;
    a1 = -8/Ke2;
    a2 = Ke3/Ke2;
    
                                
    %% Calculate Vx_i etc - here V_2ir etc act as the incident
    Vxi = V2ir + V3ir; % Vxi = v2_i + v3_i 
    Vyi = V4ir + V5ir; % Vyi = v4_i + v5_i 
    izi = V2ir + V5ir - V3ir - V4ir ; % izi = v2_i - v3_i - v4_i + V5_i 
    
    %% Calculate V_x , V_y , i_z 
    vx = (2*Vxi + Sex)/K_0;
    vy = (2*Vyi + Sey)/K_0;
    iz = izi/-2.; % i_z = izi/-2 
    
    %% 1st delay accumulator
    sedx =  -2*b0*vx + Aex1;
    sedy =  -2*b0*vy + Aey1;
    
    sex = 2*Vxi + 2*vx*(chi_e1 - 1) + sedx; 
    sey = 2*Vyi + 2*vy*(chi_e1 - 1) + sedy;
                                          
    %% 2nd delay accumulator              
    aex1 = -2*b1*vx - a1*sedx + Aex2;   
    aey1 = -2*b1*vy - a1*sedy + Aey2;
    
    %% 3rd delay accumulator
    aex2 = -2*b2*vx - a2*sedx;
    aey2 = -2*b2*vy - a2*sedy;
                                          
                                          
end
                                                  
                                                  


%% function implementation to calculate the voltage scattering in case of simple dielectric - SHUNT nodes
function [ix,iy,vz,sez] = dFC_Ez_SimpleDielectric(chi_e,...
                                                      g_e,...
                                                      V8ir,...
                                                      V9ir,...
                                                      V10ir,...
                                                      V11ir,...
                                                      Sez)
    %% Calculate K_0
    K_0 = g_e + 4.*(1.+chi_e);
    
    %% Calculate ix_i etc - here i_2ir etc act as the incident
    ixi = V8ir - V9ir; % ixi = v8_i - v9_i 
    iyi = V10ir - V11ir; % iyi = v10_i - v11_i
    Vzi = V8ir + V9ir + V10ir + V11ir ; % Vzi = v8_i + v9_i + v10_i + V11_i 
    
    %% Calculate i_x , i_y , V_z 
    ix = 1.*ixi; %
    iy = -1.*iyi; % i_y = -iyi
    vz = 2.*Vzi/K_0 +  Sez/K_0  ; % V_z = 2*Vzi/K0  + z^-1(Sez/K0)

    %% return calculate Sez */
    sez = 2.*Vzi - vz*(4.*(1.-chi_e) + g_e); % Sez = 2Vzi - Vz(4(1-chi_e) + g_e)

end










