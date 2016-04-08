% gfl = pdbread('2e8d.pdb');
%
% pos = find((strcmp{gfl.Model.Atom(:).AtomName},atom_name))
%
% X=[gfl.Model.Atom(pos).X]
% Generate a PDB file (example from MatLab help)
% gfl = getpdb('1GFL','TOFILE','1gfl.pdb')
% Read the PDB file

%Protein to be used
figures = [];
coordinatesArray = {};
normalsJan = [];
vectors = [];

%proteins = ['2LMQ'; '2LMP'; '2LMO'; '2LMN'; '2M4J']
%LEU p� 2M5N
%LEU p� 2KIB
%2NNT endast 1 stack VAL anv�nds ASN
%'2E8D' endast 2 stacks VAL anv�nds ILE
%'2RNM'; Specialstorlek
%'2N1E' : V�xer i Y-led ist�llet f�r Z-led
%;'2LNQ' : Sned i kanterna, eventuellt endast anv�nda betasheets
%'2BEG'; : M�ste ha max 1 i x-diff f�r att inte ta fel kopplingar
%'2LMN' : M�ste ha max 1 i y-diff f�r att inte ta fel kopplingar

%proteins = ['2LMP';'2KJ3';'2LMQ';'2M4J';'2MXU';'2MPZ';'2MVX']
%for p = 1:length(proteins)
    %current_proteine = proteins(p,:);
    current_proteine = '2m5k'
    filename = strcat(current_proteine, '.pdb');
    if isempty(dir(filename)) == 1 
        test = 0
        gfl = getpdb(current_proteine,'ToFile',filename);
    else
        gfl = pdbread(filename);
    end
%%    
    atom_name='CA';
    
    res_name_axis='LEU';
    %For Aromatic Ring Calcualtion
    aromatic_res = ['PHE';'TYR';'TRP';'HIS']
    %aromatic_res = ['TYR'];
    aromatic_res = cellstr(aromatic_res)

    %gfl.Model(:).Atom(:)
    searchterm = [];
    %res_name_axis = 'LEU'
    seqNr = 40;
    for s = 1:114
        searchterm = [searchterm (strcmp({gfl.Model(s).Atom(:).resName},res_name_axis))];% & (strcmp({gfl.Model(:).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(:).Atom(:).resSeq),'31')));
    end
    special = 1;
    %ST = gfl.Model(:).Atom(:)(2) == 'CA' && gfl.Model(:).Atom(:)(4) == 'ILE' 
    % Search for the couple "Atom name"
    pos = find(searchterm);
    search2 = [];
    search40 = [];
    search70 = [];
    search95 = [];
    Sequence = [];
    atoms = [];
    for p=1:114
        atoms = [atoms gfl.Model(p).Atom(:)];
    end
        for i = pos
            if (special == 1)
                if((strcmp({atoms(i).AtomName},atom_name)) == 1)

    %                 if(gfl.Model(:).Atom(i).resSeq == 70)
    %                    search70 = [search70 i]
    %                    
    %                 elseif(gfl.Model(:).Atom(i).resSeq == 95)
    %                     search95 = [search95 i]
    %                 elseif((gfl.Model(:).Atom(i).resSeq == 40))
    %                     search40 = [search40 i]
    %                 end

                    search2 = [search2 i];

                   %Sequence = [Sequence gfl.Model(p).Atom(i).resSeq];
                end
            else

               if((strcmp({gfl.Model(p).Atom(i).AtomName},atom_name)) == 1 ) 
                   search2 = [search2 i];
                   Sequence = [Sequence gfl.Model(p).Atom(i).resSeq];
               end
            end
        end
    pos = search2;

    %coords = sortrows(coords,3);
    %pos = search40
    
%     for fp = 1:length(pos)
% 
% 

    X=[atoms(pos).X];
    Y=[atoms(pos).Y];
    Z=[atoms(pos).Z]; 
    TX=transpose(X);
    TY=transpose(Y);
    TZ=transpose(Z);

    coords=[TX TY TZ];
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
%     X=[gfl.Model(:).Atom(pos).X];
%     Y=[gfl.Model(:).Atom(pos).Y];
%     Z=[gfl.Model(:).Atom(pos).Z]; 
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
%     X=[gfl.Model(:).Atom(pos).X];
%     Y=[gfl.Model(:).Atom(pos).Y];
%     Z=[gfl.Model(:).Atom(pos).Z]; 
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
    numberOfLayers = 0
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

%     for q = 2:lowestSize-1
%         P1 = coordinatesArray{1}(q,1:3);
%         P2 = coordinatesArray{2}(q,1:3);%1+SequenceLength
%         P3 = coordinatesArray{3}(q,1:3); %1+SequenceLength*2
%         normals = [normals; cross(P1-P3, P1-P2)];
%         zplane = CalculatePlane(P1,P2,P3);
%         %hold on
%         %ezmesh(zplane, [0,28,25,50])
%     end
%     normalsX = normals(:,1)
%     normalsX = mean(normalsX);
%     normalsY = normals(:,2);
%     normalsY = mean(normalsY);
%     normalsZ = normals(:,3);
%     normalsZ = mean(normalsZ);
%     normal = [normalsX normalsY normalsZ];
%     meanPointAxis = [mean(extraplot(1:3,1)), mean(extraplot(1:3,2)), mean(extraplot(1:3,3))];
%     pointMeanNormal = [[normalsX normalsY normalsZ]; meanPointAxis];
%     
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
    atomsJan = [];
    for p = 1:114
       atomsJan = [atomsJan gfl.Model(p).Atom()];
       pos = 1:(numel(gfl.Model(p).Atom()));
    end
     
        %Scan the coordinates of these atoms to respective coordinate array 
    X=[atomsJan.X];
    Y=[atomsJan.Y];
    Z=[atomsJan.Z];
    % Transpose and Genereate a single coords array
    TX=transpose(X);
    TY=transpose(Y);
    TZ=transpose(Z);
    coords=[TX TY TZ];
    % Plot a circle for every atom
    scatter3(coords(:,1),coords(:,2), coords(:,3)); hold on;
    % Create best-fit plane using Least squared method
    [x0, a] = lsplane(coords);
    point = transpose(x0);
    normal = transpose(a);
    %Calculate plane values
    d = -point*normal';
    [xx,yy]=ndgrid(-30:80,-30:80);
    z = (-normal(1)*xx - normal(2)*yy - d)/normal(3);
    % Plot plane
    surf(xx,yy,z);  hold on;
    normalsJan = [normalsJan; normal];
    %Plot a vector as axis, direction is given by best fit plane
    multiplier = 50;
    x0 = transpose(x0);
    a  = multiplier* transpose(a);
    hold on
    quiver3(x0(:,1),x0(:,2),x0(:,3),a(:,1),a(:,2),a(:,3), 'Color', 'green', 'LineWidth',3 )
    b = -a
    quiver3(x0(:,1),x0(:,2),x0(:,3),b(:,1),b(:,2),b(:,3), 'Color', 'green', 'LineWidth',3 )
    hold on
    
    hold on
    axis([-60 80 -60 80 -400 200])
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
%end
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