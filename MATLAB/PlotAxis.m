% ---------------------------------------------------------------------
% PlotAxis.m   Plots Axis and Scatter plot of all atoms in pdb struct
% ---------------------------------------------------------------------
%  Input is a point of origin and direction of the line and
%  all plots for scatter plot.
%  Output -

function PlotAxis(x0, a, gfl, coords, print)
    %Scatter coords
    if print ~= 0
        if coords == 0
            coords = CoordsGenerator(gfl,0);
            scatter3(coords(:,1),coords(:,2), coords(:,3)); hold on;
        else
           scatter3(coords(:,1),coords(:,2), coords(:,3)); hold on;
        end
        %Scaling the Axis
        multiplier = 50;
        x0 = x0;
        a  = multiplier*a;
        %Ploting the axis in both directions from the origin point
        hold on
        quiver3(x0(:,1),x0(:,2),x0(:,3),a(:,1),a(:,2),a(:,3), 'Color', 'green', 'LineWidth',3 )
        b = -a;
        quiver3(x0(:,1),x0(:,2),x0(:,3),b(:,1),b(:,2),b(:,3), 'Color', 'green', 'LineWidth',3 )
        hold on
        axis([-250 250 -250 250 -200 150])
    end