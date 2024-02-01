Rat_ID = input('What is the Rat ID (1-8) value? ');
Rat = D.animals(Rat_ID).code;

No_of_neu = length(D.animals(Rat_ID).neurons(:));

% Define the path
path = '~/Dropbox/Saptarshi_Postdoc_Projects/Project_MSR/Codes_4/data_extract/';

% Define the file name with RatCode
fileName = sprintf('Rat_%s_data_extracted.json', Rat);


% Preallocate memory for the structure array
dataStruct = repmat(struct('RAT_ID', [], 'NeuID', [], 'N_DID', [],'N_GID', [], 'LOC', [], ...
    'STAT', [], 'rSCC', [], 'bScore', [], 'hdScore', [], 'u_spiketime', [], ...
    'task', [], 'U_GID', [], 'time_range', [], 'duration', [], 'REC_Date', [], ...
    'X', [], 'Y', [], 't', [], 'HD', [], 'MDirection', [], 'Speed', [], ...
    'pmap', []), No_of_neu, 1);
% Inside the loop where you are filling in the data for each row
for i = 1:No_of_neu
    fprintf('%d  ', i);

    % Fill in the data for each row directly
    dataStruct(i).RAT_ID = Rat;
    dataStruct(i).NeuID = D.animals(Rat_ID).neurons(i).userData.Id;
    dataStruct(i).N_DID = D.animals(Rat_ID).neurons(i).id;
    dataStruct(i).N_GID = D.animals(Rat_ID).neurons(i).generalId;
    dataStruct(i).LOC = strsplit(D.animals(Rat_ID).neurons(i).userData.anatomicalLocation,' ');;
    dataStruct(i).STAT = getIfExists(D.animals(Rat_ID).neurons(i).userData, 'stats', 'NA');
    dataStruct(i).rSCC = getIfExists(D.animals(Rat_ID).neurons(i).userData, 'rSCC', 'NA');
    dataStruct(i).bScore = getIfExists(D.animals(Rat_ID).neurons(i).userData, 'bScore', 'NA');
    dataStruct(i).hdScore = getIfExists(D.animals(Rat_ID).neurons(i).userData, 'hdScore', 'NA');

    % Units of the Neurons
    units = D.animals(Rat_ID).neurons(i).members();
    n_units = length(units);
    dataStruct(i).u_spiketime = cell(1, n_units);
    dataStruct(i).task = cell(1, n_units);
    dataStruct(i).U_GID = cell(1, n_units);
    dataStruct(i).time_range = cell(1, n_units);
    dataStruct(i).duration = cell(1, n_units);
    dataStruct(i).REC_Date = cell(1, n_units);

    % Tracker Properties
    dataStruct(i).X = cell(1, n_units);
    dataStruct(i).Y = cell(1, n_units);
    dataStruct(i).t = cell(1, n_units);
    dataStruct(i).HD = cell(1, n_units);
    dataStruct(i).MDirection = cell(1, n_units);
    dataStruct(i).Speed = cell(1, n_units);

    % Plots and Inforates Properties
    dataStruct(i).pmap = cell(1, n_units);

    for j = 1:n_units
        unit = units(j);

        % Units of the Neurons
        dataStruct(i).u_spiketime{j} = unit.t;
        dataStruct(i).task{j} = unit.getParents{4,:}.sessions.type;
        dataStruct(i).U_GID{j} = unit.generalId;

        % Recording Properties
        dataStruct(i).time_range{j} = unit.getParents{4,:}.timeRange;
        dataStruct(i).duration{j} = unit.getParents{4,:}.duration';
        dataStruct(i).REC_Date{j} = unit.getParents{4,:}.date';

        % Tracker Properties
        trackers = unit.getParents{4,:}.trackers(1);
        dataStruct(i).X{j} = trackers.x;
        dataStruct(i).Y{j} = trackers.y';
        dataStruct(i).t{j} = trackers.t';
        dataStruct(i).HD{j} = trackers.hd';
        dataStruct(i).MDirection{j} =  trackers.movingDirection;
        dataStruct(i).Speed{j} = trackers.speed;

        % Plots and Inforates Properties
        dataStruct(i).pmap{j} = dManDe.plots.rateMap(unit, trackers, 'doPlot', false);
    end
end
fprintf("\nSaving the File..");
% Save the structure as a JSON file using jsonencode
jsonStr = jsonencode(dataStruct);

fullFilePath=fullfile(path, fileName);
% Open the file for writing
fid = fopen(fullFilePath, 'w');
fwrite(fid, jsonStr);
fclose(fid);

disp(['Data saved as ' fileName]);


% Helper function to get the field value if it exists, otherwise, return a default value
function value = getIfExists(structure, field, defaultValue)
    if isfield(structure, field)
        value = structure.(field);
    else
        value = defaultValue;
    end
end
