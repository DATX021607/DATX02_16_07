function[coordinatesFinished] = createCoordinatesArray(coords, pos, Restraints)

        XRestraint = Restraints(1);
        YRestraint = Restraints(2);
        ZLowRestraint = Restraints(3);
        ZHighRestraint = Restraints(4);
        numberOfLayers = 0;
        extraplot = [coords pos(:)]
        extraplot = sortrows(extraplot, 3);

        for x = 1:(length(extraplot))
           coordinatesAA1 = extraplot(x, :);
           saved = 0;
           for y1 = 1:numberOfLayers

               allCoordinatesLayer2 = coordinatesArray{y1};
               coordinatesAA2 = allCoordinatesLayer2(size(allCoordinatesLayer2, 1), :);
               aminoAcid1 = [coordinatesAA2(1); coordinatesAA2(2); coordinatesAA2(3)];
               aminoAcid2 = [coordinatesAA1(1); coordinatesAA1(2); coordinatesAA1(3)];


               distance = norm(aminoAcid1-aminoAcid2);
               Zdiff = abs(coordinatesAA1(3) - coordinatesAA2(3));
               Xdiff = abs(coordinatesAA1(1) - coordinatesAA2(1));
               Ydiff = abs(coordinatesAA1(2) - coordinatesAA2(2));

                if(Xdiff < XRestraint && Ydiff < YRestraint && Zdiff < ZHighRestraint && Zdiff > ZLowRestraint) %2RNM needs 7, 11, 15, 3.7
                      coordinatesArray{y1} = [allCoordinatesLayer2; coordinatesAA1];
                      saved = 1;
                end

           end
           if(saved == 0)
              numberOfLayers = numberOfLayers + 1;
              coordinatesArray{numberOfLayers} = coordinatesAA1;
           end
        end
        coordinatesFinished = coordinatesArray;