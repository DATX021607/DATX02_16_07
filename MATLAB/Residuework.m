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
%proteins = ['2LMQ'; '2LMP'; '2LMO'; '2LMN'; '2M4J']
proteins = ['2LMP']
%for p = 1:length(proteins)
    %current_proteine = proteins(p,:);
    current_proteine = '2LMP'
    filename = strcat(current_proteine, '.pdb');
    if isempty(dir(filename)) == 1 
        test = 0
        gfl = getpdb(current_proteine,'ToFile',filename);
    else
        gfl = pdbread(filename);
    end
    
    atom_name='CA';
    
    res_name_axis='ILE';
    %For Aromatic Ring Calcualtion
    aromatic_res = ['PHE';'TYR';'TRP';'HIS']
    %aromatic_res = ['TYR'];
    aromatic_res = cellstr(aromatic_res)

    %gfl.Model(1).Atom(:)
    if(current_proteine == '2LMP')
        res_name_axis = 'LEU'
        seqNr = 34;
        searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name_axis));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
        special = 1;
    end
    %ST = gfl.Model(1).Atom(:)(2) == 'CA' && gfl.Model(1).Atom(:)(4) == 'ILE' 
    % Search for the couple "Atom name"
    pos = find(searchterm);
    search2 = [];
    Sequence = [];
    for i = pos
        if (special == 1)
            if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1 && gfl.Model(1).Atom(i).resSeq == seqNr) 
                search2 = [search2 i];
                Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
            end
        else
            
           if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1) 
               search2 = [search2 i];
               Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
           end
        end
    end
    pos = search2;

    X=[gfl.Model(1).Atom(pos).X];
    Y=[gfl.Model(1).Atom(pos).Y];
    Z=[gfl.Model(1).Atom(pos).Z]; 

    TX=transpose(X);
    TY=transpose(Y);
    TZ=transpose(Z);

    coords=[TX TY TZ];



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
           if(distance < 6)
               if(Zdiff > 3.7)
                  coordinatesArray{y1} = [allCoordinatesLayer2; coordinatesAA1];
                  saved = 1;
               end
           end
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
    for q = 2:size(coordinatesArray{1},1)-1
        P1 = coordinatesArray{1}(q,1:3);
        P2 = coordinatesArray{1+SequenceLength}(q,1:3);%1+SequenceLength
        P3 = coordinatesArray{1+SequenceLength*2}(q,1:3); %1+SequenceLength*2
        normals = [normals; cross(P1-P3, P1-P2)];
        zplane = CalculatePlane(P1,P2,P3);
        hold on
        ezmesh(zplane, [0,28,25,50])
    end
    normalsX = normals(:,1);
    normalsX = mean(normalsX);
    normalsY = normals(:,2);
    normalsY = mean(normalsY);
    normalsZ = normals(:,3);
    normalsZ = mean(normalsZ);
    normal = [normalsX normalsY normalsZ];
    meanPointAxis = [mean(extraplot(1:3,1)), mean(extraplot(1:3,2)), mean(extraplot(1:3,3))];
    pointMeanNormal = [[normalsX normalsY normalsZ]; meanPointAxis];
    plot3(pointMeanNormal(:,1),pointMeanNormal(:,2),pointMeanNormal(:,3), 'Color', 'red', 'LineWidth',4 )

    axis([-250 250 -250 250 50 150])
    figure;
%end
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