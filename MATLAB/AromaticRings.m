% gfl = pdbread('2e8d.pdb');
%
% atom_name = ' CA ';
%
% pos = find((strcmp{gfl.Model.Atom(:).AtomName},atom_name))
%
% X=[gfl.Model.Atom(pos).X]
% Generate a PDB file (example from MatLab help)
% gfl = getpdb('1GFL','TOFILE','1gfl.pdb')
% Read the PDB file
gfl = pdbread('2lmq.pdb')
%% Define the Atom Name
atom_name='CG';
atom_name2 = 'CE1';
atom_name3 = 'CE2';
%
res_name='PHE';
%gfl.Model(1).Atom(:)
searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
%ST = gfl.Model(1).Atom(:)(2) == 'CA' && gfl.Model(1).Atom(:)(4) == 'ILE' 
% Search for the couple "Atom name"
pos = find(searchterm);
search2 = [];
Sequence = [];
CGArray = [];
CE1Array = [];
CE2Array = [];
for i = pos
   if(gfl.Model(1).Atom(i).resSeq == 19)
       if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1) 
           CGArray = [CGArray i];
       elseif((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name2)) == 1) 
           CE1Array = [CE1Array i];
       elseif((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name3)) == 1) 
           CE2Array = [CE2Array i];
       end
   end
end

posCG = CGArray;
posCE1 = CE1Array;
posCE2 = CE2Array;

XCG = [gfl.Model(1).Atom(posCG(5)).X];
YCG = [gfl.Model(1).Atom(posCG(5)).Y];
ZCG = [gfl.Model(1).Atom(posCG(5)).Z];

XCE1 = [gfl.Model(1).Atom(posCE1(5)).X];
YCE1 = [gfl.Model(1).Atom(posCE1(5)).Y];
ZCE1 = [gfl.Model(1).Atom(posCE1(5)).Z];

XCE2 = [gfl.Model(1).Atom(posCE2(5)).X];
YCE2 = [gfl.Model(1).Atom(posCE2(5)).Y];
ZCE2 = [gfl.Model(1).Atom(posCE2(5)).Z];

P1 = [XCG YCG ZCG];
P2 = [XCE1 YCE1 ZCE1];
P3 = [XCE2 YCE2 ZCE2];

normal = cross(P1-P3, P1-P2);
normalAtom = [normalsX normalsY normalsZ];
