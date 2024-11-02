%% Write data to CSV
%Parameters
ids = 1:25;
fileName = "TempOrd.csv";
address = "D:\YASEMIN\Dynamic Face Perception\DFPV1\fft\";
writeSNStoCSV(ids, fileName, address);


function writeSNStoCSV(ids, filename, address)
    % Define constants
    freqOI = 0.5:0.5:30;
    triggers = [101:104, 201:204];
    channels = [1:19, 21:25, 27:63]; %removing VEOGs, if no VEOG, electrodes next to ear

    % Initialize empty arrays for table columns
    subj = [];
    cond = [];
    chan = [];
    freq = [];
    sns = [];

    % Loop through each id and trigger to gather data
    for id = ids
        for trigger = triggers
            % Construct the file path
            filepath = address + "sns fft S" + trigger + "_P" + id + "_MAT.mat"; %change for pattern
            
            % Load the data
            snsF = load(filepath);
            sns_data = snsF.data;

            % Loop through each channel and frequency
            for ch = channels
                for fr = freqOI

                    % Calculate index based on frequency
                    ind = fr * 2 + 1;

                    sns_value = sns_data(ch, ind);
                    
                    % Append data to columns
                    subj = [subj; id];
                    cond = [cond; trigger];
                    chan = [chan; ch];
                    freq = [freq; fr];
                    sns = [sns; sns_value];
                end
            end
        end
    end

    T = table(subj, cond, chan, freq, sns);
    
    % Write table to the specified CSV file
    writetable(T, filename);
end
