% ---------------------------------------------------------------------
% BetaSolenoid.m   Calculates a line of best fit for a pdb struct
% ---------------------------------------------------------------------
%  Input is a struct defining a PDB structure
%  Output is a single line of best fit, defined by point of origin
%  and direction
%  [point, direction] = BetaSolenoid(gfl)

 function [point, direction] = LineBestFit(gfl,printAxis)
coords=CoordsGenerator(gfl, 0);
% Plot a circle for every atom
hold on;
% Create line of best fit using Least squared method
[x0, a] = ls3dline(coords);
%Plot a vector as axis, direction is given by best fit plane


%Assign final values
point = transpose(x0);
direction = transpose(a);
%Plot Axis (Optional)
PlotAxis(point, direction, gfl, coords,printAxis)
%End