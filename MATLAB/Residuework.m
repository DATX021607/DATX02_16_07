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
useBeta = 1;
usePlanes = 1;
stericZipper = 0;
XRestraint = 4.5;
YRestraint = 4;
ZLowRestraint = 1.5;
ZHighRestraint = 15;
%Protein to be used
figures = [];
coordinatesArray = {};
normalsJan = [];
vectors = [];
%proteins = ['2LMQ'; '2LMP'; '2LMO'; '2LMN'; '2M4J']

%'2LMN' : Måste ha max 1 i y-diff för att inte ta fel kopplingar
%'2N1E'; '2LNQ'; '2BEG'; '2LMN'
proteins = ['2LMP';'2KJ3';'2LMQ';'2M4J';'2MXU';'2MPZ';'2MVX'; '2KIB';'2N1E';'2M5N';'2NLQ']
%for p = 1:length(proteins)
    current_proteine = '2LMN'
    %current_proteine = proteins(p,:);
    %Check if the current protein is already scanned into matlab
    %If yes then skip this part
    %if((strcmp({gfl.Header.idCode},current_proteine)) == 0) 
     filename = strcat(current_proteine, '.pdb');
        if isempty(dir(filename)) == 1 
            test = 0
            gfl = getpdb(current_proteine,'ToFile',filename);
        else
            gfl = pdbread(filename);
        end
       % BetaSolenoid(gfl)
    %end
   
    atom_name='CA';
    if ((strcmp(current_proteine, '2KIB') == 1))
        res_name_axis = 'ILE';
        YRestraint = 15;
        ZHighRestraint = 9;
        ZLowRestraint = 0;
        XRestraint = 8;
    %elseif (strcmp(current_proteine, '2NNT') == 1)
    %   res_name_axis = 'ASN';
    %  usePlanes = 1;
    elseif (strcmp(current_proteine, '2N1E') == 1)
        res_name_axis = 'VAL'
        YRestraint = 15;
        ZHighRestraint = 4.6;
        ZLowRestraint = 0;
        XRestraint = 0.7;
        useBeta = 1;
    elseif (strcmp(current_proteine, '2M5N') == 1)
        res_name_axis = 'LEU';
        YRestraint = 15;
        ZHighRestraint = 4.6;
        ZLowRestraint = 0;
        XRestraint = 0.7;
    elseif (strcmp(current_proteine, '2LNQ') == 1)
        res_name_axis = 'ILE';
    else
        res_name_axis='VAL';
    end
    %For Aromatic Ring Calcualtion
    aromatic_res = ['PHE';'TYR';'TRP';'HIS']
     
    f = figure;
    set(f,'name',current_proteine,'numbertitle','off');
    aromatic_res = cellstr(aromatic_res)
    searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name_axis));% & (strcmp({gfl.Model(1).Atom(:).resName},res_name))  & (strcmp(int2str(gfl.Model(1).Atom(:).resSeq),'31')));
    special = 1;
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
    if useBeta == 1
        [dontcare, normal] = BetaSolenoid(gfl);
        normalT = transpose(normal);
        normalsJan = [normalsJan;normalT];
    end
    if usePlanes == 1
        [dontcare, normal] = PlaneBestFit(gfl);
        normalT = transpose(normal);
        normalsJan = [normalsJan;normalT];
    end
    
        %Generating Coords from gfl 
        coords=CoordsGenerator(gfl, pos);
        %Create the list of arrays containing the stacks
        coordinatesArray = {};
        Restraints = [XRestraint, YRestraint, ZLowRestraint, ZHighRestraint];
        coordinatesArray = createCoordinatesArray(coords, pos, Restraints);
        hold on;
        %plot3(pointMeanNormal(:,1),pointMeanNormal(:,2),pointMeanNormal(:,3), 'Color', 'red', 'LineWidth',4 )
        hold on;


        [normal] = manylines(coordinatesArray, gfl);

        vectors = [vectors;normal];

        figures = [figures f];
        hold off;
        XRestraint = 4.5;
        YRestraint = 4;
        ZLowRestraint = 1.5;
        ZHighRestraint = 15;
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