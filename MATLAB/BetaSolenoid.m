% ---------------------------------------------------------------------
% BetaSolenoid.m   Calculates a line of best fit for a pdb struct
% ---------------------------------------------------------------------
%  Input is a struct defining a PDB structure
%  Output is a single line of best fit, defined by point of origin
%  and direction
%  [point, direction] = BetaSolenoid(gfl)

 function [point, direction] = BetaSolenoid(gfl)
 
coords=CoordsGenerator(gfl, 0);
% Plot a circle for every atom
scatter3(coords(:,1),coords(:,2), coords(:,3)); 
hold on;
% Create line of best fit using Least squared method
[x0, a] = ls3dline(coords);

%Plot a vector as axis, direction is given by best fit plane
%Plot Axis (Optional)
PlotAxis(x0, a, coords)
%Assign final values
point = x0;
direction = a;
%End