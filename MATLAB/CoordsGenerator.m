function [coords]=CoordsGenerator(gfl, pos)
    if pos == 0
        pos = 1:(numel(gfl.Model(1).Atom()));
    end
    %Scan the coordinates of these atoms to respective coordinate array 
    X=[gfl.Model(1).Atom(pos).X];
    Y=[gfl.Model(1).Atom(pos).Y];
    Z=[gfl.Model(1).Atom(pos).Z];

    % Transpose and Genereate a single coords array
    TX=transpose(X);
    TY=transpose(Y);
    TZ=transpose(Z);
    coords=[TX TY TZ];