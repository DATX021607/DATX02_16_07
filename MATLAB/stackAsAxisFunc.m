function[point1, finishedvector]  = stackAsAxisFunc(coordinatesArray, gfl,printAxis)
    %if(current_proteine == '2n0aModel1Beta')
    %    stack1 = 1;
    %    stack2 = 5;
    %    stack3 = 14;
    %else
    if(strcmp(gfl.Header.idCode,'2LMN') == 1)
        stack1 = 2;
        stack2 = 4;
        stack3 = 7;
    else
        stack1 = 1;
        stack2 = 2;
        stack3 = 3;
    end

    %Generates Vector 1
    point1 = coordinatesArray{stack1}(1,1:3)
    point2 = coordinatesArray{stack1}(size(coordinatesArray{stack1},1),1:3)
    vector1 = [point2-point1;point1]
    %Generates Vector 2
    point3 = coordinatesArray{stack2}(1,1:3)
    point4 = coordinatesArray{stack2}(size(coordinatesArray{stack2},1),1:3)
    vector2 = [point4-point3;point3]
    %Generator Vector 3
    point5 = coordinatesArray{stack3}(1,1:3)
    point6 = coordinatesArray{stack3}(size(coordinatesArray{stack3},1),1:3)
    vector3 = [point6-point5;point5]

    vectorsX = [vector1(1,1);vector2(1,1);vector3(1,1)]
    vectorsY = [vector1(1,2);vector2(1,2);vector3(1,2)]
    vectorsZ = [vector1(1,3);vector2(1,3);vector3(1,3)]
    %Creates mean vector
    meanvectorX = mean(vectorsX)
    meanvectorY = mean(vectorsY)
    meanvectorZ = mean(vectorsZ)

    finishedvector = [meanvectorX meanvectorY meanvectorZ]

    PlotAxis(point1,finishedvector,gfl,0,printAxis);






