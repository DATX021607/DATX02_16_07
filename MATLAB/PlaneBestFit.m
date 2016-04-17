% ---------------------------------------------------------------------
% PlaneBestFit.m  Calculates a axis using plane-of-best fit for pdb struct
% ---------------------------------------------------------------------
%  Input is a struct defining a PDB structure
%  Output is a single line that acts as normal to the plane
%  of best fit. The line is defined by point of origin and direction 
%  [point, direction] = BetaSolenoid(gfl)

function [point, direction] = PlaneBestFit(gfl)

coords=CoordsGenerator(gfl, 0);
%Calculate point of best fit and normal to plane
[x0, a] = lsplane(coords);
point = transpose(x0);
normal = transpose(a);

%Plot Axis (Optional)
PlotAxis(x0, a, coords)
%Assign final values
hold on
axis([-250 250 -250 250 -100 150])
point = x0;
direction = a;
%End