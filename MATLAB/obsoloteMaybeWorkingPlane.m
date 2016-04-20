
 normals = [];
        lowestSize = min(size(coordinatesArray{1},1), size(coordinatesArray{2},1));
        lowestSize = min(lowestSize, size(coordinatesArray{3},1))
        for q = 2:lowestSize-1
            P1 = coordinatesArray{1}(q,1:3);
            P2 = coordinatesArray{2}(q,1:3);%1+SequenceLength
            P3 = coordinatesArray{3}(q,1:3); %1+SequenceLength*2
            normals = [normals; cross(P1-P3, P1-P2)];
            zplane = CalculatePlane(P1,P2,P3);
            %hold on
            %ezmesh(zplane, [0,28,25,50])
        end
        normalsX = normals(:,1)
        normalsX = mean(normalsX);
        normalsY = normals(:,2);
        normalsY = mean(normalsY);
        normalsZ = normals(:,3);
        normalsZ = mean(normalsZ);
        normal = [normalsX normalsY normalsZ];
        meanPointAxis = [mean(extraplot(1:3,1)), mean(extraplot(1:3,2)), mean(extraplot(1:3,3))];
        pointMeanNormal = [[normalsX normalsY normalsZ]; meanPointAxis];
