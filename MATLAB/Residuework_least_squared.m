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
    current_proteine = '2LMQ'
    filename = strcat(current_proteine, '.pdb');
    if isempty(dir(filename)) == 1 
        test = 0
        gfl = getpdb(current_proteine,'ToFile',filename);
    else
        gfl = pdbread(filename);
    end
%%   
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
    pos = 1:(numel(gfl.Model(1).Atom()));
    %Scan the coordinates of these atoms to respective coordinate array 
    X=[gfl.Model(1).Atom(pos).X];
    Y=[gfl.Model(1).Atom(pos).Y];
    Z=[gfl.Model(1).Atom(pos).Z];

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