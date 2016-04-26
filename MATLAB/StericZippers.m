function[point1,finishedvector] = StericZippers(gfl, Restraints, current_proteine, printAxis)

% gfl = pdbread('2e8d.pdb');
%
% pos = find((strcmp{gfl.Model.Atom(:).AtomName},atom_name))
%
% X=[gfl.Model.Atom(pos).X]
% Generate a PDB file (example from MatLab help)
% gfl = getpdb('1GFL','TOFILE','1gfl.pdb')
% Read the PDB file
% Which static zippers use what Amino Acid
ValSteric = {'4XFN';'4ZNN';'4NP8';'4ONK';'4OLR'};
GlySteric = {'3NHC';'3NHD'};
IleSteric = {'4R0P';'3NVF'};
LeuSteric = {'2OMP';'2OMQ'};
AsnSteric = {'3FVA';'3FTL'};
PheSteric = {'3OW9'};
SerSteric = {'3LOZ'};
MetSteric = {'3NVE'};

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
aromatic_res = ['PHE';'TYR';'TRP';'HIS'];
%aromatic_res = ['TYR'];
aromatic_res = cellstr(aromatic_res);

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

stack1 = 1;
stack2 = 2;
stack3 = 3;
testflag = 0
%Plot all atoms

point1 = coordinatesArray{stack1}(1,1:3);
point2 = coordinatesArray{stack1}(size(coordinatesArray{stack1},1),1:3);
vector1 = [point2-point1;point1];

point3 = coordinatesArray{stack2}(1,1:3);
point4 = coordinatesArray{stack2}(size(coordinatesArray{stack2},1),1:3);
vector2 = [point4-point3;point3];

point5 = coordinatesArray{stack3}(1,1:3);
point6 = coordinatesArray{stack3}(size(coordinatesArray{stack3},1),1:3);
vector3 = [point6-point5;point5];

vectorsX = [vector1(1,1);vector2(1,1);vector3(1,1)];gfl
vectorsY = [vector1(1,2);vector2(1,2);vector3(1,2)];
vectorsZ = [vector1(1,3);vector2(1,3);vector3(1,3)];

meanvectorX = mean(vectorsX);
meanvectorY = mean(vectorsY);
meanvectorZ = mean(vectorsZ);

finishedvector = [meanvectorX meanvectorY meanvectorZ];
vectors = [vectors; finishedvector];

PlotAxis(point1,finishedvector,gfl,0,printAxis);
%end