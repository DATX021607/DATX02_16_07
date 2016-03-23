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
gfl = pdbread('2lmq.pdb');

%% Define the Atom Name

atom_name='CA';
%
res_name='ILE';
%gfl.Model(1).Atom(:)
searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
%ST = gfl.Model(1).Atom(:)(2) == 'CA' && gfl.Model(1).Atom(:)(4) == 'ILE' 
% Search for the couple "Atom name"
pos = find(searchterm);
search2 = [];
Sequence = [];
for i = pos
    
   if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1) 
       search2 = [search2 i];
       Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
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



numberOfLists = 0;
extraplot = [coords pos(:)];
extraplot = sortrows(extraplot, 3);
for x = 1:(length(extraplot))
   xs = extraplot(x, :);
   saved = 0;
   for y1 = 1:numberOfLists
       
       y2 = coordinatesArray{y1};
       ys = y2(size(y2, 1), :);
       ai = [ys(1); ys(2); ys(3)];
       aj = [xs(1); xs(2); xs(3)];
        
       
       distance = norm(ai-aj);
       Zdiff = abs(xs(3) - ys(3));
       if(distance < 5.5)
           if(Zdiff > 3.7)
              coordinatesArray{y1} = [y2; xs];
              saved = 1;
           end
       end
   end
   if(saved == 0)
      numberOfLists = numberOfLists + 1;
      coordinatesArray{numberOfLists} = xs;
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
points = zeros(numberOfLists,3)
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
for i = 2:size(coordinatesArray{1},1)-1
    P1 = coordinatesArray{1}(i,1:3);
    P2 = coordinatesArray{1+SequenceLength}(i,1:3);
    P3 = coordinatesArray{1+SequenceLength*2}(i,1:3);
    normals = [normals; cross(P1-P3, P1-P2)]; 
    zplane = CalculatePlane(P1,P2,P3);
    hold on
    %ezmesh(zplane, [0,28,25,50])
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
hold on
%ezmesh(zplane, [0,maxtot,mintot,maxtot2])
%plot(normal)
hold on
%plot3(coords(:,1),coords(:,2), coords(:,3), 'o')
hold on
axis([-250 250 -250 250 50 150])

hold on
result = [];
result = [result AromaticRings('PHE', gfl, normal)];
for i = 1:length(result)
    result(i) = radtodeg(acos(result(i)));
end
hold off
