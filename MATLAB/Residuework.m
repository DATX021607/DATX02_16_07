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
% Define the Atom Name
atom_name='CA';
%
res_name='ILE';
%gfl.Model(1).Atom(:)
searchterm = (strcmp({gfl.Model(1).Atom(:).AtomName},atom_name));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
%ST = gfl.Model(1).Atom(:)(2) == 'CA' && gfl.Model(1).Atom(:)(4) == 'ILE' 
% Search for the couple "Atom name"
pos = find(searchterm);
search2 = [];
for i = pos
    
   if((strcmp({gfl.Model(1).Atom(i).resName},res_name)) == 1 && gfl.Model(1).Atom(i).resSeq == 32) 
       search2 = [search2 i];
   end
end
pos = search2;
X=[gfl.Model(1).Atom(pos).X];
Y=[gfl.Model(1).Atom(pos).Y];
Z=[gfl.Model(1).Atom(pos).Z];
save = zeros(88,88);
m1 = [];
m2 = [];
m3 = [];
for i = 1:(length(pos)-1)
    for j = (i+1):(length(pos))
       
        ai = [X(i); Y(i); Z(i)];
        aj = [X(j); Y(j); Z(j)];
        
        distance = norm(ai-aj);
        Zdiff = abs(Z(i) - Z(j));
        
        
        hold on
        if(distance < 5.27 && Zdiff > 3.7)
            Xs = [X(i) X(j)];
            Ys = [Y(i) Y(j)];
            Zs = [Z(i) Z(j)];

            saved(i,j) = 1;
            plot3(Xs, Ys, Zs,'o');
        
        end
       
    end
end

%plot3(X,Y,Z, '.')

%%
% Extract the coordinates of the Atoms matching the search criteria
% cmap = colormap(parula(length(Z)));
% scatter3(X,Y,Z,50, cmap)
% line(X,Y,Z)