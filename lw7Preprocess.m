%% Code for preprocessing on Letswave7
%bandpass, downsample, notch, interpolate, remove ocular artifact,
%rereference, epoch, exclude outlier trials, average epochs

%There should be an "Interpolation.xlsx" file in the base folder where
%three columns are subject id, channel to interpolate, channels to
%interpolate with
ids = [17]; % input ids
ppFilePath = 'C:\Users\eeglab1\Desktop\Yasemin\FaceTempOrd\Preprocessed'; % input baseFilePath
desiredEvents = {'S101','S102','S103','S104','S201','S202','S203','S204'};% input event codes

cd(ppFilePath)

%make folder for fft_mats
fftPath = fullfile(ppFilePath, 'fft');
mkdir(fftPath);

%make folder for fft_mats
epochPath = fullfile(ppFilePath, 'epoch');
mkdir(epochPath);

%make folder for fft_mats
avgPath = fullfile(ppFilePath, 'avg');
mkdir(avgPath);



%Initialize Letswave7
LW_init();

%% Preprocessing Steps
for id = ids

    filePath = [ppFilePath, '\RAW\P', num2str(id), '\'];   %each participants .eeg file should be in a root folder called RAW/P{id}
    files = dir(fullfile(filePath, '*.eeg'));
    for i = 1:length(files)
        fileName = files(i).name;
        FLW_import_data.get_lwdata('filename',fileName,'pathname',filePath,'is_save',1);
        fileName = changeName(fileName, id, i);
        
        %bandpass
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        option=struct('filter_type','bandpass','high_cutoff',100,'low_cutoff',0.5,'filter_order',4,'suffix','butt','is_save',1);
        lwdata= FLW_butterworth_filter.get_lwdata(lwdata,option);
        fileName = lwdata.header.name;
        
        %downsample
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        option=struct('x_dsratio',4,'suffix','ds','is_save',1);
        lwdata= FLW_downsample.get_lwdata(lwdata,option);
        fileName = lwdata.header.name;
        
        %notch
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        option=struct('filter_type','notch','high_cutoff',51,'low_cutoff',49,'filter_order',2,'suffix','notch','is_save',1);
        lwdata= FLW_butterworth_filter.get_lwdata(lwdata,option);
        fileName = lwdata.header.name;

        %interpolate
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        fileName = interpolate(lwdata, id, "Interpolation.xlsx");

        %remove ocular artifact
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        option=struct('ocular_channel',{{'VEOG'}},'suffix','oc_rm','is_save',1);
        lwdata= FLW_ocular_remove.get_lwdata(lwdata,option);
        fileName = lwdata.header.name;
        
        %rereferencing
        chans = {lwdata.header.chanlocs.labels};
        chans(strcmp(chans, 'VEOG')) = [];
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        option=struct('reference_list',{chans},'apply_list',{chans},'suffix','reref','is_save',1);
        lwdata= FLW_rereference.get_lwdata(lwdata,option);
        fileName = lwdata.header.name;
        
        %epoching
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        foundEvents = findEvents(lwdata, desiredEvents);
        option=struct('event_labels',{foundEvents},'x_start',0.5,'x_end',2.5,'x_duration',2,'suffix','','is_save',1);
        lwdataset= FLW_segmentation_separate.get_lwdataset(lwdata,option);
        
        
        for j = 1:length(foundEvents)
            eventCode = foundEvents{j};
            EpochedFileName = [eventCode, ' ', fileName];
            option=struct('filename',EpochedFileName);
            lwdata= FLW_load.get_lwdata(option);
            option=struct('output','power','half_spectrum',1,'suffix','fft','is_save',1);
            lwdata= FLW_FFT.get_lwdata(lwdata,option);
            exportMATfft(lwdata, fftPath);
        end
    end

    %excludeOutliers
    outliers = findOutliers(fftPath, id);
    removeOutliers(outliers, id);
end

%% Average Epochs
files = dir(fullfile(ppFilePath, '*.lw6'));
pattern = '^S\d+\s+reref\s+oc_rm\s+notch\s+ds\s+butt\s+Face_P\d+_block\d+\.lw6$';
matchingFiles = [];

% Loop through the files to find matches
for k = 1:length(files)
    fileName = files(k).name;
    if ~isempty(regexp(fileName, pattern, 'once')) % Use regexp for matching
        matchingFiles = [matchingFiles; files(k)]; % Add matching file to the list
    end
end

% Display matching files
for k = 1:length(matchingFiles)
    fileName = matchingFiles(k).name;
    option=struct('filename',fileName);
    lwdata= FLW_load.get_lwdata(option);
    option=struct('operation','average','suffix','','is_save',1);
    lwdata= FLW_average_epochs.get_lwdata(lwdata,option);
    exportMATepoched(lwdata, epochPath)
end

%% Average Blocks
averageBlocks(ids, desiredEvents, epochPath, avgPath)

%% functions
function [fileName] = interpolate(lwdata, id, excelFile)
    % Read the Excel file
    [num, txt, raw] = xlsread(excelFile); % Assuming Excel has column names

    % Extract subject-specific rows (skip header row, data starts at row 2)
    subj_column = raw(2:end, 1);  % Get the subject column, skipping headers

    % Compare numeric id with numeric values in the subject column
    subjectRows = cellfun(@(x) isequal(x, id), subj_column);
    
    if any(subjectRows)
        % Iterate over all rows for the current subject
        for i = find(subjectRows)'
            toIntp = raw{i+1, 2}; % Channel to interpolate (toIntp)
            forIntp = strsplit(raw{i+1, 3}, ','); % Channels for interpolation (forIntp), assuming comma-separated
            
            % Print the interpolation message
            fprintf('Channel %s is being interpolated using channels: %s\n', toIntp, strjoin(forIntp, ', '));
            
            % Define options for interpolation
            option = struct('channel_to_interpolate', toIntp, ...
                            'channels_for_interpolation_list', {forIntp}, ...
                            'suffix', '', ...
                            'is_save', 1);
            
            % Perform interpolation
            lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
        end
        
        % Return the modified file name after all interpolations
        fileName = lwdata.header.name;
    else
        fprintf('No rows found for subject %d in the Excel file.', id);
        fileName = lwdata.header.name;
    end
end


function [newFileName] = changeName(rawFileName, id, block)
    baseName = erase(rawFileName, '.eeg');
    
    % Define the expected file names
    lw6File = [baseName, '.lw6'];
    matFile = [baseName, '.mat'];
    
    % Define the new names
    newLw6Name = ['Face_P', num2str(id), '_block', num2str(block), '.lw6'];
    newMatName = ['Face_P', num2str(id), '_block', num2str(block), '.mat'];
    
    % Rename the .lw6 and .mat files
    movefile(lw6File, newLw6Name);
    movefile(matFile, newMatName);
    
    newFileName = newMatName;
end

function presentEvents = findEvents(lwdata, desiredEvents)
    % Extract event codes from lwdata header
    allEventCodes = {lwdata.header.events.code};  % Get all event codes from header
    presentEvents = {};
    
    % Loop over the desired event markers and check if they are in the file
    for i = 1:length(desiredEvents)
        if ismember(desiredEvents{i}, allEventCodes)
            presentEvents{end+1} = desiredEvents{i};  % Add to presentEvents if found
        end
    end
end

function exportMATfft(lwdata,fftPath)
    %export to mat (6D to 3D)
    dat = lwdata.data;
    fileName = lwdata.header.name;
    parts = strsplit(fileName);
    newFileName = sprintf('%s_%s_%s', parts{1}, parts{2}, parts{end});
    newFileName = strrep(newFileName, ' ', '_');
    data = squeeze(dat);
    fullFilePath = fullfile(fftPath, newFileName);
    save(fullFilePath, 'data');

end

function exportMATepoched(lwdata,folderPath)
    %export to mat (6D to 3D)
    dat = lwdata.data;
    fileName = lwdata.header.name;
    parts = strsplit(fileName);
    newFileName = sprintf('%s_%s', parts{1}, parts{end});
    newFileName = strrep(newFileName, ' ', '_');
    data = squeeze(dat);
    fullFilePath = fullfile(folderPath, newFileName);
    save(fullFilePath, 'data');

end

function [outliers] = findOutliers(folderPath, participantNum)
    % Function to process FFT data from .mat files, calculate 6 Hz and 7.5 Hz FFT values,
    % perform z-scoring, and identify outliers.
    fprintf('Outliers for participant %d :', participantNum)
    freq6hz_idx = 13;
    freq7_5hz_idx = 16;
    allData = [];
    files = dir(fullfile(folderPath, '*.mat'));
    pattern = sprintf('Face_P%d_block', participantNum);
    matchingFiles = files(contains({files.name}, pattern));
    
    for k = 1:length(matchingFiles)
        fileName = matchingFiles(k).name;
        pat = 'fft_S(?<trigger>\d+)_Face_P\d+_block(?<block>\d+)\.mat';
        matches = regexp(fileName, pat, 'names');
        if ~isempty(matches)
            triggerNumber = str2double(matches.trigger);
            blockNumber = str2double(matches.block);
            dat = load(fullfile(folderPath, fileName));
            data =dat.data; %CHANGE
            chan = size(data,2);
            ep = size(data,1);
            total_channels = chan;
            chans_rm_veog = total_channels - 1;
            freq_6 = data(:, 1:chans_rm_veog, freq6hz_idx);
            freq_6 = mean(freq_6,2);
            freq7_5 = data(:, 1:chans_rm_veog, freq7_5hz_idx);
            freq7_5 = mean(freq7_5,2); 
            epochs = [1:ep].';
            blocks = repelem(blockNumber, ep).';
            triggers = repelem(triggerNumber, ep).';
            comp = [blocks, triggers, epochs, freq_6, freq7_5];
            allData=[allData; comp];
        end
    end
   

    Z_6 = zscore(allData(:,4));
    Z_75 = zscore(allData(:,5));
    
    % Identify outliers (z-scores >= 3 or <= -3)
    outlier_indices = (abs(Z_6) >= 3) | (abs(Z_75) >= 3);
    bep = allData(:,1:3);
    outliers = bep(outlier_indices, :);
    disp(outliers)
    disp(['Outlier detection completed for Participant ', num2str(participantNum)]);
end

function removeOutliers(outliers, id)
    outliers = flipud(outliers);
    for n=1:size(outliers,1)
        block = outliers(n,1);
        trig = outliers(n,2);
        epoch = outliers(n,3);
        identifier = sprintf('S%d reref oc_rm notch ds butt Face_P%d_block%d.lw6', trig, id, block);
        fileName = [identifier];
        option=struct('filename',fileName);
        lwdata= FLW_load.get_lwdata(option);
        eps = size(lwdata.data,1);
        allEp = arrayfun(@num2str, 1:eps, 'UniformOutput', false);
        filteredEp = allEp(~ismember(1:eps, epoch));
        option=struct('type','epoch','items',{filteredEp},'suffix','','is_save',1);
        lwdata= FLW_selection.get_lwdata(lwdata,option);

    end
end

function averageBlocks(ids, desiredEvents, epochPath, avgPath)
    % Loop over participant IDs
    for id = ids
        % Loop over desired events
        for k = 1:length(desiredEvents)
            currentEvent = desiredEvents{k};
            matrices = {};
            % Construct the file pattern for the current event
            filePattern = fullfile(epochPath, sprintf('%s_Face_P%d_block*.mat', currentEvent, id));
            files = dir(filePattern);
            for f = 1:length(files)
                % Load the current .mat file
                fileName = fullfile(epochPath, files(f).name);
                data = load(fileName);
                matrices{end+1} = data.data;
            end
            if ~isempty(matrices)
                % Compute the average across all matrices
                avg_matrix = mean(cat(3, matrices{:}), 3);
                % Save the average matrix to a new file
                fileName = sprintf('%s_P%d.mat', currentEvent, id);
                fullFilePath = fullfile(epochPath, fileName);
                save(fullFilePath, 'avg_matrix');
                
                cd(avgPath);
                option.dimension_descriptors={'channels','X','Y','epochs'};
                option.unit='amplitude';
                option.xunit='time';
                option.yunit='frequency';
                option.xstart=0.5;
                option.xstep=0.004;
                option.ystart=0;
                option.ystep=1;
                option.is_save=1;
                
                fileName = sprintf('%s_P%d', currentEvent, id);
                option.filename=fileName;
                FLW_import_mat.get_lwdata(avg_matrix,option);


                disp(['Average matrix for participant ', num2str(id), ', trigger ', currentEvent, ' saved to ', fileName]);
            end
        end
    end
end
