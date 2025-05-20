function pObsData=readPOBStext(fileSpec,debugEnable)
% This is the common function used to read all data types

try
    filename = fileSpec.loc;
    fsize = dir(filename); fsize = fsize.bytes;
    outpath = ['../processed_POBS_Y2/',fileSpec.inst];
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    if ~exist([outpath,filesep,'parse'],'dir')
        mkdir([outpath,filesep,'parse']);
    end
    if(strcmpi(get(0,'Diary'),'off'))
        diary([outpath,filesep,'parse',filesep,datestr(now,'_yyyyddmm_HHMMSS'),'.log']);
    end
    % some housekeeping
    if(numel(filename) == 0)
        warning('backtrace','off')
        warning('readPobsText:noTaxFiles','readPobsText: No %s files found in %s',fileSpec.type,inpath);
        warning('backtrace','on')
        return;
    end

    % The sensor id number is also the file extension
    SensorID = filename(end-2:end);

    % Extract file number from name i.e. for TAX00004.005 the number is 00004
    fnum = str2double(filename(end-8:end-4));

    % Remove any files with unrecognized file name format from the list
    if isnan(fnum)
        warning('backtrace','off')
        warning('readPobsText:BadFileName','readPobsText: File name format not Recognized %s',filename(isnan(fnum)).name);
        warning('backtrace','on')
        return;
    end

    % Initialize variables
    delimiter = ',';
    eol = '\r\n';
    startRow = 1;
    endRow = inf;

    %--- Start Processing
    lastLine = []; % clear line data from last batch
    pObsData = []; % clear data from last batch

    fprintf('readPobsText: %s\n',filename(end-11:end-4))

    fileID = fopen(filename,'r');
    % Create cleanup object in case of crash
    filecleaner = onCleanup(@() fclose(fileID));
    
    firstLine = fgetl(fileID);
    if ~isempty(firstLine)
        firstLineData = textscan(firstLine,fileSpec.formatSpec,'Delimiter',delimiter,'TextType', 'string','EndOfLine', eol, 'ReturnOnError', false);
    else
        firstLineData = {[]};
    end
    % if the first line is not complete discard the data and
    % generate a warning
    if(any([cellfun(@isempty,firstLineData)]))
        warning('backtrace','off')
        warning('readPobsText:IncompleteFirstLine','readPobsText: Incomplete first line in the first file - Discarding...');
        warning('backtrace','on')
        firstLineData = [];
    end
        
    %% Read the rest of the file using textscan
    % textscan is used with ReturnOnError = true so it will return last
    % good field. The data are returned in cell arrays and any
    % incomplete frames will appear as empty elment in the cell array.
    % The last frame read is almost always an incomplete resulting in
    % array columns having different lengths.
    dataArray = textscan(fileID, fileSpec.formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', true, 'EndOfLine', eol);

    % If there are any data left in the file to read, it means that
    % textscan couldn't read the file and there is some issue with the
    % data. This section will try to read the remaining data and parse
    % into the format until the entire field is read.

    remainingbytes = fsize - ftell(fileID);
    if remainingbytes>0

        % try to scan this line into data later as well
        framesLost = 0;     % keeps count of total frames lost
        numlinesRecovered = 0;  % keeps count of lines recovered in this loop
        totalNumLinesRecovered = 0; % keeps count of total lines recovered

        % read the next line from the file to eliminate any incomplete
        % frame so that textscan can be called again on the rest of the
        % file. The recovered string will parsed into the data format
        % if possible later.
        recoveredString = fgetl(fileID);
        while(fsize - ftell(fileID)>0)
            warning('backtrace','off');
            warning('readPobsText:TextReadFail',['readPobsText: Failed to scan some parts of file %s\n',...
                'at %d bytes of %d\nAttempting to recover remaining data..'],filename(end-11:end-4),ftell(fileID),fsize);
            warning('backtrace','on')
            warning('off','readPobsText:TextReadFail')

            % try to scan the rest of the file again into the data
            % format. this data if any good will be added to the data
            % array from the first scan
            try
                recoveredDataArray = textscan(fileID, fileSpec.formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', true, 'EndOfLine', eol);
                numlinesRecovered = numel(recoveredDataArray{1});

                % If the last entry in the Array is incomplete (most likely)
                % the textscan failed originally because of this frame
                % so we'll discard it
                lastreadLen = numel(dataArray{1});
                for dd = 1:fileSpec.numFields
                    if numel(dataArray{dd}) ~= lastreadLen
                        framesLost = framesLost+1;
                        for de = 1:dd-1
                            dataArray{de}(lastreadLen) = [];
                        end
                        break;
                    end
                end
                % try to recover the line after the textscan fail. This
                % could be a complete frame. this was read earlier with
                % fgetl()
                try
                    recoveredLineData = textscan(recoveredString, fileSpec.formatSpec,'Delimiter', delimiter,'TextType', 'string','EndOfLine', eol, 'ReturnOnError', false);

                    % if the recovered line is valid add do the end of
                    % the data Array
                    if ~any(cellfun(@isempty,recoveredLineData))
                        for dd = 1:fileSpec.numFields
                            dataArray{dd}(end+1) = recoveredLineData{dd}(1);
                        end
                        totalNumLinesRecovered = totalNumLinesRecovered + 1;
                    else
                        framesLost = framesLost + 1;
                    end
                catch
                    framesLost = framesLost + 1;
                end
                recoveredString = [];

                % now let's tack on recoveredDataArray, if any, to the
                % data array
                if numlinesRecovered > 0
                    totalNumLinesRecovered = numlinesRecovered + 1;
                    for dd = 1:fileSpec.numFields
                        recoveredDataLen = numel(recoveredDataArray{dd});    % this can be different for each column
                        dataArray{dd}(end+1:end+recoveredDataLen) = recoveredDataArray{dd};
                    end
                end
            catch
                warning('readPobsText:RecoveryFailed','readPobsText: Failed...\n')
            end

            remainingbytes = fsize - ftell(fileID);

            if(remainingbytes>0)
                % read the next line to eliminate any incomplete frame - will
                % try to scan this line into data later as well
                recoveredString = fgetl(fileID);
            else
                recoveredString = [];
            end
        end
        if(totalNumLinesRecovered>0)
            fprintf('readPobsText: Success!! Recovered %d frames. At least %d frame(s) were lost\n',totalNumLinesRecovered,framesLost)
        end
        warning('on','readPobsText:TextReadFail')
        % This is just cheating to work around a case where the file
        % ends in a '.'
        if numel(recoveredString) ==1
            if recoveredString == '.'
                recoveredString = ',0.';
            end
        end
    else
        recoveredString = [];
    end
    
    lastLine = [];
    lastInd = numel(dataArray{1});
    for dd = 1:fileSpec.numFields
        if numel(dataArray{dd}) ~= lastInd
            for de = 1:dd-1
                dataArray{de}(lastInd) = [];
            end
            break;
        end
    end
    
    %% Close the text file.
    clear filecleaner;

    %% Convert the contents of columns containing numeric text to numbers.
    % Replace non-numeric text with NaN.
    % raw = repmat({''},length(dataArray{1}),length(dataArray)-1);

    lastInd = length(dataArray{1});
    %% Parse and fill up current batch table with data to save.
    % currFileTable stores data from the first (reconstructed) line and
    % the current file.
    currFileTable = table;
    % firstLineData can be two rows if both first and last lines are
    % complete.
    firstLineTable = table;
    try
        % Read data from current file into a table.
        if ~any(cellfun(@isempty,dataArray))
            currFileTable.SampleTime2K = anybaseString2num(82,dataArray{:, 1},45)*2^32 + anybaseString2num(82,dataArray{:, 2},45);
            for dd = fileSpec.nonDataFields+1:fileSpec.numFields
                currFileTable.(fileSpec.fieldNames{dd-fileSpec.nonDataFields}) = (dataArray{:, dd});
                if debugEnable
                    if all(isnan(dataArray{:,dd}))
                        warning('backtrace','off')
                        warning('readPobsText:allNanTexts','One or more data columns are NaN\nMake sure ReadDataWithTimeStamp flag is appropriately set')
                        warning('backtrace','off')
                        warning('off','readPobsText:allNanTexts')
                    end
                end
            end
        end
        if(~isempty(firstLineData))
            firstLineTable.SampleTime2K = anybaseString2num(82,firstLineData{:,1},45)*2^32 + anybaseString2num(82,firstLineData{:,2},45);
            for dd = fileSpec.nonDataFields+1:fileSpec.numFields
                firstLineTable.(fileSpec.fieldNames{dd-fileSpec.nonDataFields}) = firstLineData{:,dd};
            end
            % Insert data from current file and first line into data
            % table for the batch
            pObsData = [pObsData;firstLineTable;currFileTable];
        else
            % Insert data from current file into data table for the
            % batch
            pObsData = [pObsData;currFileTable];
        end
    catch err
        keyboard;
    end
    %% Save the data
    % End of current batch - save it to mat file
    firstFileNameOnly = filename(end-11:end-4);
    if exist([outpath,filesep,firstFileNameOnly,'.mat'],'file')
        warning('backtrace','off');
        warning('readPobsText:outputFilesExist','readPobsText: Found existing mat files. Any existing files will be overwritten')
        warning('off','readPobsText:outputFilesExist');
        warning('backtrace','on');
    end
    % save([outpath,filesep,firstFileNameOnly,'.mat'],'pObsData','SensorID');
    disp('Done!')
catch errMsg
    if debugEnable
        fprintf('\n\n!!readPobsText: Encountered an error!!\n');
        fprintf('readPobsText: Make sure ReadDataWithTimeStamp flag is appropriately set\n');
        fprintf('Error:\n%s\n',errMsg.message);
        fprintf('%s\n',errMsg.identifier);
        for jj = 1:numel(errMsg.stack)
            fprintf('In Function %s, Line %d\n',errMsg.stack(jj).name,errMsg.stack(jj).line)
        end
        fprintf('readPobsText: Halting for debugging\n');
        keyboard
    else
        fprintf('Set the debugEnable flag to stay in the function workspace after a crash\n')
        rethrow(errMsg);
    end

end