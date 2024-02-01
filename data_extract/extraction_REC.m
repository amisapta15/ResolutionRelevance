Rat_ID = input('What is the Rat ID (1-8) value? ');
Rat = D.animals(Rat_ID).code;

No_of_rec = length(D.animals(Rat_ID).recordings(:));

% Define the path
path = '~/Dropbox/Saptarshi_Postdoc_Projects/Project_MSR/Codes_4/data_extract/';

% Define the file name with RatCode
fileName = sprintf('Rat_%s_RECdata_extracted.json', Rat);


% Preallocate memory for the structure array
dataStruct = repmat(struct('RAT_ID', [], 'REC_GID', [], 'REC_duration', [],'REC_timerange', [], 'REC_Date', [], ... 
    'N_DIDs', [],'Neu_IDs', [], 'N_GIDs', [],'N_LOCs', [], 'task', [], 'U_GIDs', []), No_of_rec, 1);
% Inside the loop where you are filling in the data for each row
for i = 1:No_of_rec
    fprintf('%d  ', i);

    % Fill in the data for each row directly
    dataStruct(i).RAT_ID = Rat;
    dataStruct(i).REC_GID = D.animals(Rat_ID).recordings(i).generalId;
    dataStruct(i).REC_duration = D.animals(Rat_ID).recordings(i).duration;
    dataStruct(i).REC_timerange = D.animals(Rat_ID).recordings(i).timeRange;
    dataStruct(i).REC_Date = D.animals(Rat_ID).recordings(i).date;
    dataStruct(i).task = D.animals(Rat_ID).recordings(i).sessions(1).name;

    n_units = length(D.animals(Rat_ID).recordings(i).units);
    
    dataStruct(i).N_DIDs = cell(1, n_units);
    dataStruct(i).Neu_IDs = cell(1, n_units);
    dataStruct(i).N_GIDs = cell(1, n_units);
    dataStruct(i).N_LOCs = cell(1, n_units);
    dataStruct(i).U_GIDs = cell(1, n_units);
    for j = 1:n_units
        dataStruct(i).N_DIDs{j} = D.animals(Rat_ID).recordings(i).sessions(1).userData.containedNeurons.ids(j);
        
        neus=D.animals(Rat_ID).recordings(i).units(j).getParents{7,:};
        dataStruct(i).Neu_IDs{j} = neus.userData.Id;
        dataStruct(i).N_GIDs{j} = neus.generalId;

        llocs=D.animals(Rat_ID).recordings(i).sessions(1).userData.containedNeurons.locs(j);
        dataStruct(i).N_LOCs{j} = strsplit(llocs{:},' ');
        dataStruct(i).U_GIDs{j} =  D.animals(Rat_ID).recordings(i).units(j).generalId;
    end

end
fprintf("\nSaving the File:");
% Save the structure as a JSON file using jsonencode
jsonStr = jsonencode(dataStruct);

fullFilePath=fullfile(path, fileName);
% Open the file for writing
fid = fopen(fullFilePath, 'w');
fwrite(fid, jsonStr);
fclose(fid);

disp(['Data saved as ' fileName]);


% % Helper function to get the field value if it exists, otherwise, return a default value
% function value = getIfExists(structure, field, defaultValue)
%     if isfield(structure, field)
%         value = structure.(field);
%     else
%         value = defaultValue;
%     end
% end
