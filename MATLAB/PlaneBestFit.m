% ---------------------------------------------------------------------
% PlaneBestFit.m  Calculates a axis using plane-of-best fit for pdb struct
% ---------------------------------------------------------------------
%  Input is a struct defining a PDB structure
%  Output is a single line that acts as normal to the plane
%  of best fit. The line is defined by point of origin and direction 
%  [point, direction] = BetaSolenoid(gfl)

function [point, direction] = PlaneBestFit(gfl,printAxis)
coords=CoordsGenerator(gfl, 0);
%Calculate point of best fit and normal to plane
[x0, a] = lsplane(coords);
point = transpose(x0);
direction = transpose(a);
%Plot Axis (Optional)
PlotAxis(point, direction, gfl, coords,printAxis)
%Assign final values

%End