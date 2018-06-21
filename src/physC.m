
%% just some constants
classdef physC
    
    properties(Constant)
        EPS0=8.854187817e-12;
        MU0=4e-7*pi;

        NU0 = sqrt(physC.MU0/physC.EPS0); 
        c0=1/sqrt(physC.MU0*physC.EPS0);
    end
    
end