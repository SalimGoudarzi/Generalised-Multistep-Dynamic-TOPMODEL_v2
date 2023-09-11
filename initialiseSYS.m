%**************************************************************************
function [Sx,Su,dtw]=initialiseSYS(Nc)
%--------------------------------------------------------------------------
%                            initialise variables
%--------------------------------------------------------------------------
Sx  = zeros(Nc,1) + 1e-16;
%maximum subsurface flow per unit area
dtw = zeros(Nc,1) + 1e-16;
%unsaturated zone initial storage 
Su  = zeros(Nc,1) + 1e-16;
