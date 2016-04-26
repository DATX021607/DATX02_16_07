function[finishedvector, f] = StericZippers(gfl, Restraints, current_proteine)

% gfl = pdbread('2e8d.pdb');
%
% pos = find((strcmp{gfl.Model.Atom(:).AtomName},atom_name))
%
% X=[gfl.Model.Atom(pos).X]
% Generate a PDB file (example from MatLab help)
% gfl = getpdb('1GFL','TOFILE','1gfl.pdb')
% Read the PDB file
ValSteric = {'4XFN';'4ZNN';'4NP8';'4ONK';'4OLR'};
GlySteric = {'3NHC';'3NHD'};
IleSteric = {'4ROP';'3NVF'};
LeuSteric = {'2OMP';'2OMQ'};
AsnSteric = {'3FVA';'3FTL'};
PheSteric = {'3OW9'};
SerSteric = {'3LOZ'};
MetSteric = {'3NVE'};
%Protein to be used
figures = [];
coordinatesArray = {};
normalsJan = [];
vectors = [];

   
    atom_name='CA';
    
    if  any(ismember(ValSteric, current_proteine)) == 1
        res_name_axis ='VAL';
    elseif any(ismember(GlySteric, current_proteine)) == 1
        res_name_axis = 'GLY';
    elseif any(ismember(AsnSteric,current_proteine)) == 1
        res_name_axis = 'ASN';
    elseif any(ismember(PheSteric, current_proteine)) == 1
        res_name_axis = 'PHE';
    elseif any(ismember(SerSteric, current_proteine)) == 1
        res_name_axis = 'SER';
    elseif any(ismember(MetSteric, current_proteine)) == 1
        res_name_axis = 'MET';
    elseif any(ismember(IleSteric, current_proteine)) == 1
        res_name_axis = 'ILE';
    elseif any(ismember(LeuSteric, current_proteine)) == 1
        res_name_axis = 'LEU';
    end

    %For Aromatic Ring Calcualtion
    aromatic_res = ['PHE';'TYR';'TRP';'HIS']
    %aromatic_res = ['TYR'];
    aromatic_res = cellstr(aromatic_res)

    %gfl.Model(:).Atom(:)
    searchterm = [];
    %res_name_axis = 'LEU'
    seqNr = 40;
    for s = 1:size(gfl.Model(:),1)
        searchterm = [searchterm (strcmp({gfl.Model(s).Atom(:).resName},res_name_axis))];% & (strcmp({gfl.Model(:).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(:).Atom(:).resSeq),'31')));
    end
    special = 1;
    %ST = gfl.Model(:).Atom(:)(2) == 'CA' && gfl.Model(:).Atom(:)(4) == 'ILE' 
    % Search for the couple "Atom name"
    pos = find(searchterm);
    search2 = [];
    Sequence = [];
    atoms = [];
    for p=1:size(gfl.Model(:),1)
        atoms = [atoms gfl.Model(p).Atom(:)];
    end
        for i = pos
                if((strcmp({atoms(i).AtomName},atom_name)) == 1)
                    search2 = [search2 i];
                end
        end
    pos = search2;
    X=[atoms(pos).X];
    Y=[atoms(pos).Y];
    Z=[atoms(pos).Z]; 
    TX=transpose(X);
    TY=transpose(Y);
    TZ=transpose(Z);

    coords=[TX TY TZ];

    coordinatesArray = createCoordinatesArray(coords, pos, Restraints);

    hold on

    

    f = figure;
    set(f,'name',current_proteine,'numbertitle','off');
    hold on;
    %plot3(pointMeanNormal(:,1),pointMeanNormal(:,4),pointMeanNormal(:,3), 'Color', 'red', 'LineWidth',4 )
    hold on;
    
   
    stack1 = 1;
    stack2 = 2;
    stack3 = 3;

    atomsJan = [];
    for p = 1:size(gfl.Model(:),1)
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
    axis([-60 80 -60 80 -300 100])
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