function f = colebrookRN(Re,epsilon,D)
% __________________________________________________________
% This function computes the friction factor, using
% the Colebrook-White equation. Friction factor is
% computing with the Newton-Rhapson method.
%
% Inputs:
% Re = Reynolds number (-)
% epsilon = Pipe Rugosity (mm)
% D = Pipe diameter (mm)
% 
% Output:
% f = friction factor (-)
%
% Example:
% f = colebrookRN(1.97E5,0.015,76.2)
%
% Reference:
% Ladino M., E.O, Garcia U., C.A. y Garcia V., M.C. (2019). 
%     Darcy-Weisbach resistance coefficient determination 
%     using Newton-Raphson approach for android 4.0. 
%     Tecnura, 23(60), 52-58. 
%     DOI: https://doi.org/10.14483/22487638.14929
% 
% Gabriel Ruiz
% 2023
% ___________________________________________________________

narginchk(3,3)
erel = epsilon/D;
nMaxIter = 180;

if Re >= 2300

    % Computing the inital value from transient flow
    flprandtl = epsilon/(3.7 * D); % using Nikuradse
    C = flprandtl;
    ld = -0.86 * log(C);
    fseed = (1/ld)^2;
    f = fseed;
    i = 0;

    while i < nMaxIter
        flprandtl = (erel/3.7);
        flkarman = 2.51/(Re*sqrt(f));

        f_f = (1/sqrt(f)) + (2*log10(flprandtl+flkarman));
        fprime  = (-0.5 * f^-1.5) +....
                  ((2*((-2.51/(2*Re))*f^-1.5)*log10(epsilon))/...
                  (flprandtl+flkarman));

        fn  = f - (f_f/fprime);
        i = i + 1;
        %fprintf('Iter: %i, f = %5.4f\n',i,fn);

        if abs(f-fn) > 0.0001
            f = fn;
        else
            f = round(f*1000)/1000;
            break;
        end
    end
else
	f = 64/Re;  % using Poiselle
end