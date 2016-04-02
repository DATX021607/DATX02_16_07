function [crossproduct] = CalculateNormal (P1, P2, P3)
    crossproduct = cross(P1-P3, P1-P2);
    crossproduct = crossproduct*10;
   if(crossproduct(3) < 0)
        crossproduct = crossproduct .* -1;
    end
    Point1 = transpose(P1(:));
    Point2 = transpose(P2(:));
    Point3 = transpose(P3(:));
    points = [Point1;Point2;Point3];
    meanpoint = [mean(points(:,1)), mean(points(:,2)), mean(points(:,3))];
    pointAromatic = [crossproduct; meanpoint ]
    quiver3(pointAromatic(2,1), pointAromatic(2,2), pointAromatic(2,3),pointAromatic(1,1), pointAromatic(1,2), pointAromatic(1,3))
    
    hold on
    fill3(points(:,1), points(:,2), points(:,3), 'o')
    
    hold on
end
    