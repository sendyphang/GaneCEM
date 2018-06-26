function ZTLM2Dkernel_doScattering(grid)
    %% check polarisation 
    polarisation = grid.pol_type;
    
    %% scattering
    if strcmp(polarisation , 'Hz') 
        
        temp_V2ir = grid.V2ir;
        grid.V2ir = grid.V_x + grid.i_z - grid.V3ir; % V_2ir acts as the reflected V2r = Vx + iz - V3i */
        grid.V3ir = grid.V_x - grid.i_z - temp_V2ir; % V3r = Vx - iz - v2i -> use the temporary saved one */

        temp_V4ir = grid.V4ir;
        grid.V4ir = grid.V_y - grid.i_z - grid.V5ir; % V4r = Vy - iz - V5i */
        grid.V5ir = grid.V_y + grid.i_z - temp_V4ir; % V5r = Vy + iz - V4i */
        
    else
        
        temp_V8ir = grid.V8ir;
        grid.V8ir = grid.V_z - grid.i_x - grid.V9ir; % V_8ir acts as the reflected V8r = Vz - ix - V9i */
        grid.V9ir = grid.V_z + grid.i_x - temp_V8ir; % V9r = Vz + ix - v8i -> use the temporary saved one */

        temp_V10ir = grid.V10ir;
        grid.V10ir = grid.V_z + grid.i_y - grid.V11ir; % V10r = Vz + iy - V11i */
        grid.V11ir = grid.V_z - grid.i_y - temp_V10ir; % V11r = Vz - iy - V10i */
    end
end