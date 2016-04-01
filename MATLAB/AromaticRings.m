function [result] = AromaticRings (ResName, gfl, normal)

    if(strcmp(ResName, 'PHE') == 1 || strcmp(ResName, 'TYR') == 1)
        atom_name='CG';
        atom_name2 = 'CE1';
        atom_name3 = 'CE2';
    elseif(strcmp(ResName, 'TRP') == 1)
        atom_name = 'CD2';
        atom_name2 = 'CZ2';
        atom_name3 = 'CZ3';
    elseif(strcmp(ResName, 'HIS') == 1)
        ResName
        atom_name = 'CG';
        atom_name2 = 'CD2';
        atom_name3 = 'CE1';    
    end

    %gfl.Model(1).Atom(:)
    searchterm = (strcmp({gfl.Model(1).Atom(:).resName},ResName));
    % Search for the couple "Atom name"
    pos = find(searchterm);
    search2  = [];
    Sequence = [];
    CGArray  = [];
    CE1Array = [];
    CE2Array = [];
    chainID  = [];
    seqNum = [];
    resultV  = [];
    result   = [];
    for i = pos

       if((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name)) == 1) 
           CGArray = [CGArray i];
           %Stores ChainID for each atom extracted
           chainID = [chainID gfl.Model(1).Atom(i).chainID]
           seqNum = [seqNum gfl.Model(1).Atom(i).resSeq]
       elseif((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name2)) == 1) 
           CE1Array = [CE1Array i];
       elseif((strcmp({gfl.Model(1).Atom(i).AtomName},atom_name3)) == 1) 
           CE2Array = [CE2Array i];
       end

    end

    posCG = CGArray;
    posCE1 = CE1Array;
    posCE2 = CE2Array;
    for i = 1:length(CGArray)
        XCG = [gfl.Model(1).Atom(posCG(i)).X];
        YCG = [gfl.Model(1).Atom(posCG(i)).Y];
        ZCG = [gfl.Model(1).Atom(posCG(i)).Z];

        XCE1 = [gfl.Model(1).Atom(posCE1(i)).X];
        YCE1 = [gfl.Model(1).Atom(posCE1(i)).Y];
        ZCE1 = [gfl.Model(1).Atom(posCE1(i)).Z];

        XCE2 = [gfl.Model(1).Atom(posCE2(i)).X];
        YCE2 = [gfl.Model(1).Atom(posCE2(i)).Y];
        ZCE2 = [gfl.Model(1).Atom(posCE2(i)).Z];

        P1Aro = [XCG YCG ZCG];
        P2Aro = [XCE1 YCE1 ZCE1];
        P3Aro = [XCE2 YCE2 ZCE2];

        calculatedNormal = CalculateNormal(P1Aro,P2Aro,P3Aro);
        
        resultV = [resultV dot(calculatedNormal, normal)/(norm(calculatedNormal)*norm(normal))]

    end
    if (length(resultV)>0)
        for i = 1:length(resultV)
            result(i) = radtodeg(acos(resultV(i)));
        end 
        C1 = num2cell(transpose(result));
        C2 = cellstr (transpose(chainID));
        C3 = num2cell(transpose(seqNum));
        result = horzcat(C2,C3,C1)
    end
end

