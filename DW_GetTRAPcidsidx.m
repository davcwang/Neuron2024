function [TRAPcids, TRAPidx, age, ExptNo, N, OID, OO, OOforA, GOVA, depthcutoffs] = DW_GetTRAPcidsidx(Expt,NPXSpikes)

% ExptNo goes from 0801000 to 0911000 to 1304000, where last 3 digits are for idx (not cid)
% The TRAP ids for APC/PPC specifically is for OlfacFigures (where it's run individually by APC/PPC)
% The TRAP ids for an entire expt are for spont, CCGs, etc

if contains(Expt,'8-1') && contains(Expt,'PPC') && ~contains(Expt,'8-10')
    TRAPcids = [374]; % 8-1-PPC
    age = 6;
    ExptNo = 0801000;

    depthcutoffs = [0 300 550 1500];

    GOVA = "GO";

    N = 14; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex"; "GO"; "MO"; "GO ventrum"; "GO milk"; "Maternal ur"; "Male ur"; "GO AF"; "VA"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 14 9 10 13 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-3') && contains(Expt,'PPC')
    TRAPcids = [91]; % 8-3-PPC
    age = 9;
    ExptNo = 0803000;

    depthcutoffs = [0 50 250 1150];

    GOVA = "GO";

    N = 14; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex"; "GO"; "MO"; "GO ventrum"; "GO milk"; "Maternal ur"; "Male ur"; "Bedding/nesting"; "VA"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 14 9 10 13 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-4') && contains(Expt,'PPC')
    TRAPcids = [306 664]; % 8-4-PPC
    age = 11;
    ExptNo = 0804000;

    depthcutoffs = [0 89 253 1420];

    GOVA = "GO";

    N = 14; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex"; "GO"; "MO"; "GO ventrum"; "GO milk"; "Maternal ur"; "Male ur"; "Bedding/nesting"; "VA"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 14 9 10 13 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-7') && contains(Expt,'PPC')
    TRAPcids = [35]; % 8-7-PPC
    age = 7;
    ExptNo = 0807000;

    depthcutoffs = [0 250 600 1200];

    GOVA = "VA";

    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "VA"; "GO"; "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 9 13 10 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-8') && contains(Expt,'PPC')
    TRAPcids = [11]; % 8-8-PPC
    age = 8;
    ExptNo = 0808000;

    depthcutoffs = [0 100 300 1200];

    GOVA = "VA";
    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "VA"; "GO"; "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 9 13 10 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-10') && contains(Expt,'APC')
    TRAPcids = [943]; % 8-10-APC
    age = 10;
    ExptNo = 0810000;

    depthcutoffs = [0 50 125 850];

    GOVA = "VA";

    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "VA"; "GO"; "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 9 13 10 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'8-10') && contains(Expt,'PPC')
    TRAPcids = [10 18 411 32 923]; % 8-10-PPC
    age = 10;
    ExptNo = 0810000;

    depthcutoffs = [0 150 450 1500];

    GOVA = "VA";

    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti";  "VA"; "GO"; "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [8 1 2 3 4 5 6 7 9 13 10 11 12]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-3') % && contains(Expt,'APC')
    TRAPcids = [710 381 110 403]; % 9-3-APC
    age = 9;
    ExptNo = 0903000;

    depthcutoffs = [0 150 350 1250];

    GOVA = "GO";

    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex"; "MO"; "GO ventrum"; "GO milk"; "Maternal ur"; "Male ur"; "Male ventrum"; "VA ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [7 1 2 3 4 5 6 13 8 12 9 10 11]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-4') % && contains(Expt,'APC')
    TRAPcids = [22 63]; % 9-4-APC
    age = 8;
    ExptNo = 0904000;

    depthcutoffs = [0 196 493 1213];

    GOVA = "VA";

    N = 13; %number of unique odors in experiment
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "GO ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"; "VA ventrum"]; %this is the original non-permuted order of odor IDs 
    OO = [7 1 2 3 4 5 6 13 8 12 9 10 11]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-6')
    TRAPcids = [1011 2042 2374 2104]; % 9-6
    age = 8;
    ExptNo = 0906000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 100 275 2031];
        TRAPcids = [11];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 95 299 1100];
        TRAPcids = [42 374 104];
    end

    GOVA = "VA";

    N = 14;
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex"; "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"; "GO ventrum"; "GO milk"];
    OO = [7 1 2 3 4 5 6 8 13 12 9 14 10 11]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-7')
    TRAPcids = [1370 1774 1089 2092]; % 9-7
    age = 10;
    ExptNo = 0907000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 121 254 1714];
        TRAPcids = [370 774 89];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 1 2 3]; % this looks like OT, exclude
        % depthcutoffs = [0 100 300 1000];
        TRAPcids = [92];
    end

    GOVA = "VA";

    N = 14;
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"; "GO ventrum"; "GO milk"];
    OO = [7 1 2 3 4 5 6 8 13 12 9 14 10 11]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-8')
    TRAPcids = [10274 10455 10313]; % 9-8
    age = 10;
    ExptNo = 0908000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 150 350 1500];
        TRAPcids = [274 455 313];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 250 450 1000];
    end

    GOVA = "VA";

    N = 14;
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "VA ventrum"; "VA milk"; "Maternal ur"; "Male ur"; "Male ventrum"; "GO ventrum"; "GO milk"];
    OO = [7 1 2 3 4 5 6 8 13 12 9 14 10 11]; %Odor order (OO) vector; write elements in desired order

elseif contains(Expt,'9-9')
    TRAPcids = [1234]; % 9-9
    age = 6;
    ExptNo = 0909000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 100 200 1200];
        TRAPcids = [234];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 400 1000 1000];
    end

    GOVA = "GO";

    N = 14;
    OO = [7 1 2 3 4 5 6 8 9 10 11 12 13 14]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "Maternal ur"; "Male ur"]; % DW-XI-9-9 to 9-11

elseif contains(Expt,'9-10')
    TRAPcids = [1506 1078 1762 2011]; % 9-10
    age = 8;
    ExptNo = 0910000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 100 300 1300];
        TRAPcids = [506 78 762];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 250 500 1100];
        TRAPcids = [11];
    end

    GOVA = "GO";

    N = 14;
    OO = [7 1 2 3 4 5 6 8 9 10 11 12 13 14]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "Maternal ur"; "Male ur"]; % DW-XI-9-9 to 9-11

elseif contains(Expt,'9-11')
    TRAPcids = [1067 1509 1688 2365 2046]; % 9-11
    age = 9;
    ExptNo = 0911000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 48 156 1250];
        TRAPcids = [67 509 688];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 170 544 1459];
        TRAPcids = [365 46];
    end

    GOVA = "GO";

    N = 14;
    OO = [7 1 2 3 4 5 6 8 9 10 11 12 13 14]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Et ti"; "2-hex";  "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "Maternal ur"; "Male ur"]; % DW-XI-9-9 to 9-11

elseif contains(Expt,'10-2')
    TRAPcids = [1011 1043 1042]; % 10-2
    age = 6;
    ExptNo = 1002000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 250 550 1250];
        TRAPcids = [11 43 42];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 1 2 3];
    end

    GOVA = "VA";

    N = 14;
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'10-3') && contains(Expt,'APC')
    TRAPcids = [14 26]; % 10-3-APC
    age = 7;
    ExptNo = 1003000;

    depthcutoffs = [0 154 389 1191];

    GOVA = "VA";

    N = 14;
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'11-6')
    TRAPcids = [2291 2311 2070 2073 2517 2095 2112 2343 2044]; % 11-6
    age = 8;
    ExptNo = 1106000;

    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 71 215 1062];
    elseif contains(Expt,'APC')
        TRAPcids = [291 311 70 73 517 95 112 343 44];
        depthcutoffs = [0 45 120 811];
    end

    GOVA = "VA";

    N = 14;
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'12-1')
    TRAPcids = [10032 10129 20020]; % 12-1
    age = 7;
    ExptNo = 1201000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 154 464 1629];
        TRAPcids = [32 129];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 141 409 909];
        TRAPcids = [20];
    end

    GOVA = "GO";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'12-2')
    TRAPcids = [10592 10861 20011 20853]; % 12-2
    age = 8;
    ExptNo = 1202000;

    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 147 504 1543];
        TRAPcids = [592 861];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 156 394 1286];
        TRAPcids = [11 853];
    end

    GOVA = "GO";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'12-3')
    TRAPcids = [10653 10098 10102 10501 20374 20009 20403 20102 20223]; % 12-3 PPC
    age = 9;
    ExptNo = 1203000;

    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 121 389 2331];
        TRAPcids = [653 98 102 501];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 191 527 1268];
        TRAPcids = [374 9 403 102 223];
    end

    GOVA = "GO";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'13-1')
    TRAPcids = [103 115 139 191 193 841 537 538 542 846 545 249]; % 13-1
    age = 40;
    ExptNo = 1301000;

    depthcutoffs = [0 256 586 1248];

    GOVA = "GO";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'13-2')
    TRAPcids = [10089 20082]; % 13-2
    age = 40;
    ExptNo = 1302000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 252 681 1789];
        TRAPcids = [89];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 176 411 1331];
        TRAPcids = [82];
    end

    GOVA = "GO";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'13-3')
    TRAPcids = []; % 13-3
    age = 40;
    ExptNo = 1303000;

    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 276 635 1382];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 213 497 1155];
    end

    GOVA = "VA";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

elseif contains(Expt,'13-4')
    TRAPcids = [10465 10109 20070 20106]; % 13-4
    age = 40;
    ExptNo = 1304000;
    if contains(Expt,'PPC') % if it's NPXspikes from a single probe (ie. when getting Z scores and not using synced probes)
        depthcutoffs = [0 218 488 2484];
        TRAPcids = [465 109];
    elseif contains(Expt,'APC')
        depthcutoffs = [0 237 618 1773];
        TRAPcids = [70 106];
    end

    GOVA = "VA";

    N = 14; % this is for maternal/sigma only
    OO = [7 1 2 3 4 8 9 10 11 12 13 14 5 6]; %Odor order (OO) vector; write elements in desired order
    OID = ["Iso ac"; "Et ac"; "Hex"; "Et bu"; "Maternal ur"; "Male ur"; "MO"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"];

end

% Get TRAPidx from TRAPcids
TRAPidx = [];
if nargin > 1 % if there's an NPXSpikes (if not, I'm only trying to get Expt info, not TRAP idx)
    for i = 1:length(TRAPcids) % find TRAP idx
        TRAPidx = [TRAPidx find(NPXSpikes.cids==TRAPcids(i))];
    end
end

% Getting OOforA based on desiredOID (the OID for DW-XI-12) based on valve numbers from respective expts (OO)
% Determine desiredOID based on GO or VA, so familiar mother odors appear first
if GOVA == "GO"
    desiredOID = ["MO"; "Iso ac"; "Et ac"; "Hex"; "Et bu"; "GO ventrum"; "VA ventrum"; "Male ventrum"; "GO milk"; "VA milk"; "GO AF"; "VA AF"; "Maternal ur"; "Male ur"];
elseif GOVA == "VA"
    desiredOID = ["MO"; "Iso ac"; "Et ac"; "Hex"; "Et bu"; "VA ventrum"; "GO ventrum"; "Male ventrum"; "VA milk"; "GO milk"; "VA AF"; "GO AF"; "Maternal ur"; "Male ur"];
end

% Generate an OOforA, which is the OO according to DW-XI-12 desired odors,
% but based on OIDs from each prior expt; 99 = odor isn't present
OOforA = [];
for i = 1:length(desiredOID) % loop through the desired 14 odors in DW-XI-12

    if ismember(desiredOID(i),OID) % if the i-th desired odor is in OID (the current experiment's odor string array)
        OOforA = [OOforA find(OID == desiredOID(i))]; % find the index of in the OID that matches the i-th desired odor
    
    else % if the i-th desired odor isn't in OID
        OOforA = [OOforA 0]; % use 0 as a "missing odor" index
    end
end

