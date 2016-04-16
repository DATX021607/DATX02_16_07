% gfl = pdbread('2e8d.pdb');
%
% pos = find((strcmp{gfl.Model.Atom(:).AtomName},atom_name))
%
% X=[gfl.Model.Atom(pos).X]
% Generate a PDB file (example from MatLab help)
% gfl = getpdb('1GFL','TOFILE','1gfl.pdb')
% Read the PDB file


%--------------------------------
%Available Functions
%[point, direction] = PlaneBestFit(gfl)
%[point, direction] = BetaSolenoid(gfl)
%[coords]=CoordsGenerator(gfl)
%PlotAxis(x0, a, coords)
%-------------------------------

%Protein to be used
figures = [];
coordinatesArray = {};
normalsJan = [];
vectors = [];
%proteins = ['2LMQ'; '2LMP'; '2LMO'; '2LMN'; '2M4J']
%LEU på 2M5N
%LEU på 2KIB
%2NNT endast 1 stack VAL används ASN
%'2E8D' endast 2 stacks VAL används ILE
%'2RNM'; Specialstorlek
%'2N1E' : Växer i Y-led istället för Z-led
%;'2LNQ' : Sned i kanterna, eventuellt endast använda betasheets
%'2BEG'; : Måste ha max 1 i x-diff för att inte ta fel kopplingar
%'2LMN' : Måste ha max 1 i y-diff för att inte ta fel kopplingar
  
proteins = ['2LMP';'2KJ3';'2LMQ';'2M4J';'2MXU';'2MPZ';'2MVX']
for p = 1:length(proteins)
    %current_proteine = proteins(p,:);
    current_proteine = '1Z2F'
    %Check if the current protein is already scanned into matlab
    %If yes then skip this part
    if((strcmp({gfl.Header.idCode},current_proteine)) == 0) 
        filename = strcat(current_proteine, '.pdb');
        if isempty(dir(filename)) == 1 
            test = 0
            gfl = getpdb(current_proteine,'ToFile',filename);
        else
            gfl = pdbread(filename);
        end
        BetaSolenoid(gfl)
    end
%%    
    atom_name='CA';
    
    res_name_axis='VAL';
    %For Aromatic Ring Calcualtion
    aromatic_res = ['PHE';'TYR';'TRP';'HIS']
    %aromatic_res = ['TYR'];
    aromatic_res = cellstr(aromatic_res)

    %gfl.Model(1).Atom(:)

    %res_name_axis = 'LEU'
    seqNr = 40;
    searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name_axis));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
    special = 1;
    %ST = gfl.Model(1).Atom(:)(2) == 'CA' && gfl.Model(1).Atom(:)(4) == 'ILE' 
    % Search for the couple "Atom name"
    pos = find(searchterm);
    search2 = [];
    search40 = [];
    search70 = [];
    search95 = [];
    Sequence = [];
    for i = pos
        if (special == 1)
            if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1)
                
%                 if(gfl.Model(1).Atom(i).resSeq == 70)
%                    search70 = [search70 i]
%                    
%                 elseif(gfl.Model(1).Atom(i).resSeq == 95)
%                     search95 = [search95 i]
%                 elseif((gfl.Model(1).Atom(i).resSeq == 40))
%                     search40 = [search40 i]
%                 end
                    
                search2 = [search2 i];
                
               Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
            end
        else
            
           if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1 ) 
               search2 = [search2 i];
               Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
           end
        end
    end
    
    pos = search2;

    %coords = sortrows(coords,3);
    %pos = search40
    
%     for fp = 1:length(pos)
% 
% 
%Generating Coords from gfl 
coords=CoordsGenerator(gfl);

%     if(length(coordinatesArray) < 1)
%         coordinatesArray{1} = coords(fp,:);
%     else
%         coordinatesArray{1} = [coordinatesArray{1}; coords(fp,:)]
%     end
%     end
%     pos = search70;
%     for fp = 1:length(pos)
%         
% 
%     X=[gfl.Model(1).Atom(pos).X];
%     Y=[gfl.Model(1).Atom(pos).Y];
%     Z=[gfl.Model(1).Atom(pos).Z]; 
% 
%     TX=transpose(X);
%     TY=transpose(Y);
%     TZ=transpose(Z);
% 
%     coords=[TX TY TZ];
%     if(length(coordinatesArray) < 2)
%         coordinatesArray{2} = coords(fp,:);
%     else
%         coordinatesArray{2} = [coordinatesArray{2}; coords(fp,:)]
%     end
%     end
%     pos = search95;
%     for fp = 1:length(pos)
%         
% 
%     X=[gfl.Model(1).Atom(pos).X];
%     Y=[gfl.Model(1).Atom(pos).Y];
%     Z=[gfl.Model(1).Atom(pos).Z]; 
% 
%     TX=transpose(X);
%     TY=transpose(Y);
%     TZ=transpose(Z);
% 
%     coords=[TX TY TZ];
%         %coordinatesArray{3} = [coordinatesArray{3}; coords(pos,:)]
%     if(length(coordinatesArray) < 3)
%         coordinatesArray{3} = coords(fp,:);
%     else
%         coordinatesArray{3} = [coordinatesArray{3}; coords(fp,:)]
%     end
%     end
    numberOfLayers = 0;
    extraplot = [coords pos(:)];
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
           %if(distance < 6)
            %   if(Zdiff > 3.7)
            if(Xdiff < 4.5 && Ydiff < 4 && Zdiff < 15 && Zdiff > 1.5) %2RNM needs 7, 11, 15, 3.7
                  coordinatesArray{y1} = [allCoordinatesLayer2; coordinatesAA1];
                  saved = 1;
            end
           %end
       end
       if(saved == 0)
          numberOfLayers = numberOfLayers + 1;
          coordinatesArray{numberOfLayers} = coordinatesAA1;
       end
    end



    %plot3(X,Y,Z, '.')


    % Extract the coordinates of the Atoms matching the search criteria
    % cmap = colormap(parula(length(Z)));
    % scatter3(X,Y,Z,50, cmap)

    % line(X,Y,Z)

    hold on

    max1 = max(coordinatesArray{1}(:,1));
    max2 = max(coordinatesArray{2}(:,1));
    max3 = max(coordinatesArray{3}(:,1));
    maxtot = max([max1,max2,max3]);
    B = linspace(0, maxtot + 5);
    max1 = max(coordinatesArray{1}(:,2));
    max2 = max(coordinatesArray{2}(:,2));
    max3 = max(coordinatesArray{3}(:,2));
    maxtot2 = max([max1,max2,max3]);

    min1 = min(coordinatesArray{1}(:,2));
    min2 = min(coordinatesArray{2}(:,2));
    min3 = min(coordinatesArray{3}(:,2));
    mintot = min([min1,min2,min3]);
    C = linspace(mintot-5,maxtot2+10);
    coordinates = {};
    points = zeros(numberOfLayers,3)
    %for x = 1:numberOfLists
    %    coordinates{x} = test{x}(1,:)
    %    points(x,:) = coordinates{x}(1:3)

    %     coordinates2{x} = test{x}(size(test{x},1),:)
    %     points2(x,:) = coordinates2{x}(1:3)
    %     
    %     coordinates3{x} = test{x}(size(test{x},1)/2,:)
    %     points3(x,:) = coordinates3{x}(1:3)

    %end
    normals = [];
    SequenceLength = length(unique(Sequence));
    lowestSize = min(size(coordinatesArray{1},1), size(coordinatesArray{2},1));
    lowestSize = min(lowestSize, size(coordinatesArray{3},1))
    
%     for q = 2:size(lowestSize)-1
%         P1 = coordinatesArray{1}(q,1:3);
%         P2 = coordinatesArray{1+SequenceLength}(q,1:3);%1+SequenceLength
%         P3 = coordinatesArray{1+SequenceLength*2}(q,1:3); %1+SequenceLength*2
%         normals = [normals; cross(P1-P3, P1-P2)];
%         zplane = CalculatePlane(P1,P2,P3);
%         hold on
%         ezmesh(zplane, [0,28,25,50])
%     end

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
    
    f = figure;
    set(f,'name',current_proteine,'numbertitle','off');
    hold on;
    %plot3(pointMeanNormal(:,1),pointMeanNormal(:,2),pointMeanNormal(:,3), 'Color', 'red', 'LineWidth',4 )
    hold on;
    
    if(length(current_proteine) > 5)
        
        if(current_proteine == '2n0aModel1Beta')
            stack1 = 1;
            stack2 = 5;
            stack3 = 14;
        else
            stack1 = 1;
            stack2 = 2;
            stack3 = 3;
        end
    else
        stack1 = 1;
        stack2 = 2;
        stack3 = 3;
    end
    
    %Use Plane of Best Fit method to find axis
    [x0, a] = PlaneBestFit(gfl);
    
    
point1 = coordinatesArray{stack1}(1,1:3)
point2 = coordinatesArray{stack1}(size(coordinatesArray{stack1},1),1:3)

vector1 = [point2-point1;point1]

point3 = coordinatesArray{stack2}(1,1:3)
point4 = coordinatesArray{stack2}(size(coordinatesArray{stack2},1),1:3)

vector2 = [point4-point3;point3]

point5 = coordinatesArray{stack3}(1,1:3)
point6 = coordinatesArray{stack3}(size(coordinatesArray{stack3},1),1:3)

vector3 = [point6-point5;point5]

vectorsX = [vector1(1,1);vector2(1,1);vector3(1,1)]
vectorsY = [vector1(1,2);vector2(1,2);vector3(1,2)]
vectorsZ = [vector1(1,3);vector2(1,3);vector3(1,3)]

meanvectorX = mean(vectorsX)
meanvectorY = mean(vectorsY)
meanvectorZ = mean(vectorsZ)

finishedvector = [meanvectorX meanvectorY meanvectorZ]
finishedvector = finishedvector.*3;

vectors = [vectors; finishedvector];
quiver3(coordinatesArray{stack1}(1,1), coordinatesArray{stack1}(1,2),coordinatesArray{stack1}(1,3),finishedvector(1), finishedvector(2), finishedvector(3),'Color','red', 'LineWidth',4)




figures = [figures f];
hold off;
end
savefig(figures,'figures.fig')
 %%   %-------------------------------------------------------------
        %FOR PRINTING
        %hold on
    for i=1:length(aromatic_res)
        curr_res=aromatic_res(i,:);
        result = [];
        %Result is Cell matrix with ChainID : Angle
       

        result = [result AromaticRings(curr_res, gfl, normal); ];
        if (length(result)>0)
            %FOR PRINTING
            %hold off
            %Data extracion
            % Create Table Titles
            title = ['Chain ID  ';'Sequence #';'Angle     '];
            %Combines function data witl table titles
            celltitle = transpose(cellstr(title));
            dataToWrite = vertcat(celltitle, result);
            %Creat tabe name as protein - Residue
            tabName=strcat(current_proteine, ' - ',curr_res);
            %Writes Data (only works if excel is available atm)
            xlswrite('Angles',dataToWrite,tabName{1});
        end
    end
    f = figure;
    figures = [figures f];
    
    
%end
savefig(figures,'figures.fig')