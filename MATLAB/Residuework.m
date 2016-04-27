%--------------------------------
% Available Functions
% [point, direction] = PlaneBestFit(gfl)
% [point, direction] = BetaSolenoid(gfl)
% [coords]=CoordsGenerator(gfl)
% PlotAxis(x0, a, gfl,  coords)
%-------------------------------
% Steric Zippers and which way they grow
% X: 4OLR, 4ONK, 3NVE, 2OMP
% Y: 4NP8, 2OMQ, 3FTL, 3FVA, 3OW9, 3NHC, 4ZNN, 4R0P
% Z: 3NVF, 3NHD
%-------------------------------
% Initiate variables
figures = [];
coordinatesArray = {};
normal= [];
vectors = [];
%-------------------------------
% Initiate Flags
steric       = 0;
useLine      = 0;
usePlane    = 0;
useStackAxis = 0;
%-------------------------------
%Print figures for 
printAromatic = 0;
printAxis = 1;
%-------------------------------
% Default constraints for Stack as Axis methd
XRestraint = 4.5;
YRestraint = 4;
ZLowRestraint = 1.5;
ZHighRestraint = 15;
%-------------------------------
% Protein to be used
% proteins = ['2LMQ'; '2LMP'; '2LMO'; '2LMN'; '2M4J']
%Problematiska Zippers : '3NVG';'3NVH';'4XFN';'3LOZ'
stericZippers = {'4ZNN';'3NHC';'3NHD';'4R0P';
                 '3NVF';'3OW9';'3LOZ';'2OMP';
                 '3NVE';'3FVA';'3FTL';'2OMQ';
                 '4NP8';'4ONK';'4OLR'};
%'2LMN' : Måste ha max 1 i y-diff för att inte ta fel kopplingar
%'2N1E'; '2LNQ'; '2BEG'; '2LMN' ;'3LOZ'

proteins = ['2LMP';'2KJ3';'2LMQ';'2M4J';
            '2MXU';'2MPZ';'2MVX';'2KIB';
            '2N1E';'2M5N';'2LNQ';'2NNT';
            '4ZNN';'3NHC';'3NHD';
            '4R0P';'3NVF';'3OW9';
            '2OMP';'3NVE';'3FVA';'3FTL';
            '2OMQ';'4NP8';'4ONK';'4OLR'];
%endast zippers 
%proteins = ['4ZNN';'3NHC';'3NHD';'4R0P';'3NVF';'3OW9';'2OMP';'3NVE';'3FVA';'3FTL';'2OMQ';'4NP8';'4ONK';'4OLR'];
% Configuration regarding which method to use
% Use Plane of best fit
planeofBestFit = ['2E8D';'2MVX';'2MPZ';'2LMP';
                  '2M4J';'2LMQ';'2KJ3'];
% Use Line of best fit
lineofBestFit  = ['2N1E';'2KIB';'2MXU';'2RNM'];
% Use Stack as Axis
stackAsAxis    = ['2MXU';'2BEG';'2M5N';'2LNQ';'2NNT'];
% Is Beta Solenoid
betaSolenoid = ['4S37';'2N3D';'4IHG';'1TYU'];
% Tell program to use Alpha Carbon atoms in calculations.
atom_name='CA';
% For Aromatic Ring Calcualtion
aromatic_res = ['PHE';'TYR';'TRP';'HIS'];
for p = 1:length(proteins)
    %------------------
    % Reset flags
    steric = 0;
    useLine      = 0;
    usePlane    = 0;
    useStackAxis = 0;
    %------------------
    % Reset Constraints
    XRestraint = 4.5;
    YRestraint = 4;
    ZLowRestraint = 1.5;
    ZHighRestraint = 15;
    %------------------
    %Crashade på 2M5N, kollar varför
    current_proteine = proteins(p,:);
    % Check if the current proteine is a Steric zipper, if yes use pdb1
    % , i.e. Biological assembly. Also set Steric to 1.
    if ismember(current_proteine,stericZippers)
        filename = strcat(current_proteine, '.pdb');
        gfl = pdbread(filename);
        steric = 1;
    else
        % Check if .pdb exist, otherwise fetch it. Sets Steric to 0.
        filename = strcat(current_proteine, '.pdb');
        if isempty(dir(filename)) == 1 
            gfl = getpdb(current_proteine,'ToFile',filename);
        else
            gfl = pdbread(filename);
        end
        steric = 0;
    end
    if steric ~= 1
        % If member of lineofBestFit use Line of Best Fit method
        if ismember(current_proteine,lineofBestFit,'rows') || ismember(current_proteine,betaSolenoid,'rows')
            useLine = 1;
        % If member of planeofBestFit use Plane of Best Fit method
        elseif ismember(current_proteine,planeofBestFit,'rows')
            usePlane = 1;
        % If member of stackASAxis use stackAsAxis Fit method
        elseif ismember(current_proteine,stackAsAxis,'rows')
            useStackAxis = 1;
        % If not member of any of these arrays, break and tell user to
        % specify method
        elseif useLine ~= 1 || usePlane ~= 1 || useStackAxis ~=1
            error('Unclassefied protine, please specify method') 
        end
    end
    
    % Protein specific constraints
    if ((strcmp(current_proteine, '2KIB') == 1))
        res_name_axis = 'ILE';
        YRestraint = 15;
        ZHighRestraint = 9;
        ZLowRestraint = 0;
        XRestraint = 8;
        useLine = 1;
    elseif (strcmp(current_proteine, '2NNT') == 1)
        res_name_axis = 'THR';
        usePlane = 1;
    elseif (strcmp(current_proteine, '2N1E') == 1)
        res_name_axis = 'VAL';
        YRestraint = 15;
        ZHighRestraint = 4.6;
        ZLowRestraint = 0;
        XRestraint = 0.7;
        useLine = 1;
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

    % Create a figure with name of the current proteine
    if printAxis ~= 0
        f = figure;
        set(f,'name',current_proteine,'numbertitle','off');
    end
    aromatic_res = cellstr(aromatic_res);
    searchterm = (strcmp({gfl.Model(1).Atom(:).resName},res_name_axis));
    special = 1;
    % Search for the couple "Atom name"
    pos = find(searchterm);
    % Reset search lists
    search2 = [];
    Sequence = [];
    % Extract atoms from positions find in pos
    for i = pos
        if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1)
           search2 = [search2 i];
           Sequence = [Sequence gfl.Model(1).Atom(i).resSeq];
        end
    end
    
    pos = search2;
    % Use Line of best Fit method
    if useLine == 1
        [dontcare, normal] = LineBestFit(gfl, printAxis);
        normalT = transpose(normal);
        normal = [normalsJan;normalT];
    end
    % Use Plane of best Fit method
    if usePlane == 1
        [dontcare, normal] = PlaneBestFit(gfl,printAxis);
        normalT = transpose(normal);
        normal = [normalsJan;normalT];
    end
    % Use Stack as Axis method
    if useStackAxis == 1
        if steric == 1
           % Apply Special Restraints
           Restraints = [XRestraint, YRestraint, ZLowRestraint, ZHighRestraint];
           % Call function StericZippers to calculate normal
           [x0,normal] = StericZippers(gfl, Restraints, current_proteine, printAxis);
           
        else
            % Generating Coords from gfl 
            coords=CoordsGenerator(gfl, pos);
            % Create the list of arrays containing the stacks
            coordinatesArray = {};
            Restraints = [XRestraint, YRestraint, ZLowRestraint, ZHighRestraint];
            coordinatesArray = createCoordinatesArray(coords, pos, Restraints)
            [x0,normal] = stackAsAxisFunc(coordinatesArray, gfl, printAxis);
        end
    end
    if strcmp(current_proteine, '2M5N') || strcmp(current_proteine,'2LNQ') || strcmp(current_proteine, '2NNT')
        normal = transpose(normal);
    end
        vectors = [vectors;normal];

        figures = [figures f];
        hold off;

%% 
        % FOR PRINTING
        hold on
    for i=1:length(aromatic_res)
        curr_res=aromatic_res(i,:);
        result = [];
        % Result is Cell matrix with ChainID : Angle
        result = [result AromaticRings(curr_res, gfl, normal ,steric, current_proteine, printAromatic) ];
        if (length(result)>0)
            %FOR PRINTING
            % Data extracion
            % Create Table Titles
            % Leave the padding
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
    figures = [figures f];
end