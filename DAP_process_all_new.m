clear;

% change working directory to this expt's
cd('/Documents/MATLAB/PROJ');

%% Read in Participant Info

% read in participant info sheet
participant_info = readtable('Participants_PROJ.csv');

% rename vars
participant_info.Properties.VariableNames = {'Subject_No','Include'};

%% Set up some info for use later

% Some Info:
screen_size = [1920,1080]; % first number is x direction,
    % starting in top left corner going right, second numer is in y
    % direction, starting in top left corner going down
sample_rate = 1000; %sample rate of EyeTracker (Hz)
len_wb_scrs = 60; %60 seconds
center_pixels = 400; %number of pixels to consider "center"

recall_time_ms = 3000; %time in milliseconds of pause from the end of 
    %passage/participant's final key press until they are prompted to
    %recall the passage

% may want later for naming and saving things, etc.
folder = '/Users/nicoleayasse 1/Documents/MATLAB/PROJ';

        
%% What Sample Messages mean (for reference)

% 1 = 'MODE RECORD CR 1000'
% 2 = 'Whitescreen60'
% 3 = 'Blackscreen60'
% 4 = 'BaselineStart'
% 5 = 'BaselineEnd'
% 6 = 'Instructions for main displayed'
% 7 = 'PassageInitiated'
% 8 = 'DispCross_PassBegin'
% 9 = 'Main clause sound file begins playing'
% 10 = 'Key Pressed for clause'
% 11 = 'DispDot_PassEnd'
% 12 = 'Start Recalling (pupil end)'
% 13 = 'TRACKER_TIME'


        
%% Initialize some matrices/arrays to add to later

% due to the difficulties of calculating the proper number of rows, or
% dealing with extra rows of NaN's later, these have been made empty and
% vertcat is used at the end
key_lat_allpart = [];
mat_final_allpart = [];
prerecall_mat_allpart = [];
baselineavg_perpass_allpart = [];
baselines_chunks_allpart = [];
baselines_passend_allpart = [];
pupil_wbinfo_allpart = [];
mat_subj_bins_pass_start_allpart = [];
mat_subj_bins_pass_end_allpart = [];
mat_subj_bins_chunk_start_allpart = [];
mat_subj_peak_passend_allpart = [];
mat_subj_peak_chunkend_allpart = [];



%% Loop Through Participants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for participant = 1:size(participant_info,1)
    
if participant_info{participant,2} ~= 0

    subject_no = participant_info{participant,1}
    if isnumeric(subject_no) == 1
        subject_str = int2str(subject_no); %change subject number to string
    elseif isnumeric(subject_no) == 0
        subject_str = subject_no{1};
    end
    
    % if subject is 'testtwo', make it zero to make processing easier
    if strcmp(subject_str, 'testtwo') == 1
        subject_no = 0;
    else
        subject_no = str2double(subject_no{1});
    end
    
    
% Load Sample Report data
sdir = strcat(folder,'/','PROJ_onesubjprocess_',subject_str,'.mat'); %put the name all together
dir_files = dir(sdir); %look for the file that fits that description
    %(in folder, any .mat file with the subject number in it)
file_name = dir_files(1).name; %find the file name of the (hopefully) only file that fits
general_name_samprep = file_name(1:end-4); %char of file name minus the .mat, for later
    %general_name possibly not needed??
load(file_name); %load the data into workspace (table name is tbl_subj_final)


% number of trials
no_passages = nanmax(tbl_subj_final{:,'Passage'});

% calc total number of chunks/clauses for whole expt/session
no_total_chunks = size(tbl_resfile,1);

%%% Create a new column of the chunk # within the entire passage (this is
%%% OverallNo for RO 1, but it continues counting through running orders)
tbl_subj_final{:,'ChunkNo_Overall'} = tbl_subj_final{:,'OverallNo'} - ...
    ((tbl_subj_final{1,'RO'}-1)*no_total_chunks);


%% Create map of variable names, and convert to matrix 
% for faster/more efficient working when needed

% replace SorC
toreplace_strs = {'S','C'};
forreplace_nums = 1:length(toreplace_strs); 
tbl_subj_final{:,'SorCNum'} = NaN;
for i = 1:length(toreplace_strs)
    tf = contains(tbl_subj_final{:,'SorC'},toreplace_strs{i});
        %a logical array of indices
    tbl_subj_final{tf,'SorCNum'} = forreplace_nums(i);
    %replace column SorCNum with 1 or 2 based on value in SorC column
end

% replace HLExp
toreplace_strs = {'L','H'};
forreplace_nums = 1:length(toreplace_strs); 
tbl_subj_final{:,'HLExpNum'} = NaN;
for i = 1:length(toreplace_strs)
    tf = contains(tbl_subj_final{:,'HLExp'},toreplace_strs{i});
        %a logical array of indices
    tbl_subj_final{tf,'HLExpNum'} = forreplace_nums(i);
    %replace column HLExpNum with 1 or 2 based on value in HLExp column
end

% replace EorN
toreplace_strs = {'E','N'};
forreplace_nums = 1:length(toreplace_strs); 
tbl_subj_final{:,'EorNNum'} = NaN;
for i = 1:length(toreplace_strs)
    tf = contains(tbl_subj_final{:,'EorN'},toreplace_strs{i});
        %a logical array of indices
    tbl_subj_final{tf,'EorNNum'} = forreplace_nums(i);
    %replace column EorNNum with 1 or 2 based on value in EorN column
end


%%% So that the columns to be added to the matrix have names in the
%%% mapping:

% pre-alloc new columns in original table 
tbl_subj_final{:,'Time_Rel_PassStart'} = NaN;
tbl_subj_final{:,'Time_Rel_Chunk'} = NaN;
tbl_subj_final{:,'Time_Rel_PassEnd'} = NaN;
% pre-alloc new columns in original table
tbl_subj_final{:,'Pupil_Size_ScBa'} = NaN;
% pre-alloc new col that indicates yes/no baseline period
tbl_subj_final{:,'Baseline'} = 0;
% pre-alloc new col that indicates yes/no prerecall period
tbl_subj_final{:,'Pre_Recall'} = 0;


% Table having the numeric variables of tbl_subj_final:
S = vartype('numeric');
num_subj_final = tbl_subj_final(:,S);

% Given your variable names, and assuming they map in order from columns 1 
% through N, you can create a mapping like so:
    % varNames={'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};
    % col = containers.Map(varNames, 1:numel(varNames));
varNames = num_subj_final.Properties.VariableNames;
col = containers.Map(varNames, 1:numel(varNames));
%%% ALSO: save varNames for use later, to re-convert to table
% And now you can use the map to access columns of your data by 
% variable name. For example, if you want to fetch the columns for 
% variables A and C (i.e. the first and third) from a matrix data, 
% you would do this:
% subData = data(:, [col('A') col('C')]);
% https://stackoverflow.com/questions/44679592/matlab-table-dataset-...
% type-optimization

% make a matrix version of the numeric-only table
mat_subj_final = table2array(num_subj_final);




%% Separate Out White and Black Screens

size_wb_scr = (len_wb_scrs * sample_rate); %for pre-alloc
white_scr = NaN(size_wb_scr,6); %1:SubjNo; 2:RO; 3:SampleIndex (NEW); 
        %4:X; 5:Y; 6: PupSize
black_scr = NaN(size_wb_scr,6);

for iteration = 1:size(mat_subj_final,1) %loop through Data
    
    if mat_subj_final(iteration,col('Trial_Index'))==1 %Trial Index == 1 
                %(wb scr section)
        if mat_subj_final(iteration,col('Sample_Message')) == 2 %White Scr 
                    %Starts message
            white_scr(:,1) = tbl_subj_final{1,'Subject_No'}; 
                %put SubjNo in first col
            white_scr(:,2) = tbl_subj_final{1,'RO'}; %put RO in second col
            white_scr(:,3) = 1:size_wb_scr; %creating an index
            white_scr(:,4:6) = mat_subj_final(iteration:(iteration+...
                size_wb_scr-1),[col('Gaze_X') col('Gaze_Y') ...
                col('Pupil_Size')]); %put X, Y, & PupSize in cols 4-6
            if mat_subj_final(iteration+size_wb_scr-1,...
                col('Sample_Message'))== 3 %if that last index has 
                    %reached Black Scr message
                white_scr(end,4:6) = NaN;
            end
        end %end if White Scr Start message
        if mat_subj_final(iteration,col('Sample_Message')) == 3 %Black Scr 
                    %Starts message
            black_scr(:,1) = tbl_subj_final{1,'Subject_No'}; 
                %put SubjNo in first col
            black_scr(:,2) = tbl_subj_final{1,'RO'}; %put RO in second col
            black_scr(:,3) = 1:size_wb_scr; %creating an index
            black_scr(:,4:6) = mat_subj_final(iteration:(iteration+...
                size_wb_scr-1),[col('Gaze_X') col('Gaze_Y') ...
                col('Pupil_Size')]); %put X, Y, & PupSize in cols 4-6
            if mat_subj_final(iteration+size_wb_scr-1,...
                col('Sample_Message'))== 4 %if that last index has 
                    %reached next message
                black_scr(end,4:6) = NaN;
            end
        end %end if Black Scr Start message
    end %end if trial index == 1
    
end %end for loop through mat_subj_final


% calculate Pupil Dynamic Range from Pupil Max Black and Pupil Min White
pupil_maxblack = nanmean(black_scr((len_wb_scrs*.75)*sample_rate:end,6));
pupil_minwhite = nanmean(white_scr((len_wb_scrs*.75)*sample_rate:end,6));
pupil_dr = pupil_maxblack - pupil_minwhite;

pupil_wbinfo = [tbl_subj_final{1,'Subject_No'}, pupil_minwhite, ...
    pupil_maxblack, pupil_dr];


%% Key Press

% create an output matrix
key_lat = NaN(no_total_chunks,13);
    % 1 subj no
    % 2 RO
    % 3 trial index
    % 4 overall number
    % 5 passage
    % 6 chunk
    % 7 SorC
    % 8 HLExp
    % 9 Prob
    % 10 EorN
    % 11 Processing
    % 12 File Length
    % 13 key press latency

% calculate key press latency from sample report instead of results file
for chunk = 1:no_total_chunks
    
    % initialize to NaN
    mclause_start = NaN;
    keypress = NaN;
    
    % subset that table to have only current chunk
    current_data = tbl_subj_final(tbl_subj_final.ChunkNo_Overall==chunk,:);

    for sample = 1:size(current_data,1) %loop through Data
    
                
            if current_data{sample,'Sample_Message'}==9 %Main clause plays
                mclause_start = current_data{sample,'Time_Stamp'};
                    % save this as time stamp of main clause starting
            end
            if current_data{sample,'Sample_Message'}==10 %key press
                keypress = current_data{sample,'Time_Stamp'};
                    % save this as time stamp of key press after clause
            end

            if isnan(keypress) == 0 % if keypress has been assigned

                % convert filelength from sec to msec
                filelength_ms = current_data{sample,'FileLength'} *1000;
                % calculate time from END of clause playing to key press
                keypress_latency = keypress - (mclause_start+filelength_ms);

                key_lat(chunk,13) = keypress_latency; %put that in matrix to save

                % put other info into the matrix to save
                key_lat(chunk,1) = subject_no;
                key_lat(chunk,2) = current_data{sample,'RO'};
                key_lat(chunk,3) = current_data{sample,'Trial_Index'};
                key_lat(chunk,4) = current_data{sample,'ChunkNo_Overall'};
                key_lat(chunk,5) = current_data{sample,'Passage'};
                key_lat(chunk,6) = current_data{sample,'Chunk'};
                if strcmpi(current_data{sample,'SorC'},'S')
                    key_lat(chunk,7) = 1;
                elseif strcmpi(current_data{sample,'SorC'},'C')
                    key_lat(chunk,7) = 2;
                end
                if strcmpi(current_data{sample,'HLExp'},'L')
                    key_lat(chunk,8) = 1;
                elseif strcmpi(current_data{sample,'HLExp'},'H')
                    key_lat(chunk,8) = 2;
                end
                key_lat(chunk,9) = current_data{sample,'Prob'};
                if strcmpi(current_data{sample,'EorN'},'E')
                    key_lat(chunk,10) = 1;
                elseif strcmpi(current_data{sample,'EorN'},'N')
                    key_lat(chunk,10) = 2;
                end
                key_lat(chunk,11) = current_data{sample,'Processing'};
                key_lat(chunk,12) = current_data{sample,'FileLength'};

                break %break out of sample loop

            end % if keypress has been assigned
        
    end % end sample loop
    
end % end chunk loop


%% Make Relative Times

% Sample Messages (for reference):
% 1 = 'MODE RECORD CR 1000'
% 2 = 'Whitescreen60'
% 3 = 'Blackscreen60'
% 4 = 'BaselineStart'
% 5 = 'BaselineEnd'
% 6 = 'Instructions for main displayed'
% 7 = 'PassageInitiated'
% 8 = 'DispCross_PassBegin'
% 9 = 'Main clause sound file begins playing'
% 10 = 'Key Pressed for clause' (next clause)
% 11 = 'DispDot_PassEnd'
% 12 = 'Start Recalling (pupil end)'
% 13 = 'TRACKER_TIME'

%%% to subset columns wanted from matrix:
% subData = data(:, [col('A') col('C')]);

% % new columns have been pre-alloc in original table so matrix has the
% % mapped names:
% tbl_subj_final{:,'Time_Rel_PassStart'} = NaN;
% tbl_subj_final{:,'Time_Rel_Chunk'} = NaN;
% tbl_subj_final{:,'Time_Rel_PassEnd'} = NaN;



for pass = 1:no_passages
    
    % time relative to passage start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==8,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 8 which is display cross
        %passage begins (7 = passage initiated)
    pass_start_time = mat_subj_final(start_iter_pass,col('Time_Stamp'));
        %get the time stamp for when the passage started
        
    % time relative to passage end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==11,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 11 which is display dot
        %passage end (7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 12=start recalling)
    pass_end_time = mat_subj_final(end_iter_pass,col('Time_Stamp'));
        %get the time stamp for when the passage started
        
    % find time of end of pre-recall period
    end_pup_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12=start recalling 
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot passage
        %end)
    
    for iteration = start_iter_pass:end_pup_pass
        
        % make time relative to passage start for that passage trial
        mat_subj_final(iteration,col('Time_Rel_PassStart')) = ... 
            mat_subj_final(iteration,col('Time_Stamp')) - pass_start_time;
        
        % make time relative to passage end for that passage trial
        mat_subj_final(iteration,col('Time_Rel_PassEnd')) = ...
            mat_subj_final(iteration,col('Time_Stamp')) - pass_end_time;
        
    end %end iteration (sample) index for time rel to pass start/end
    
end %end passage index



for chunk = 1:no_total_chunks
    
    % time relative to chunk start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==9,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 9 which is 'Main clause sound 
        %file begins playing' (8 = passage begins, 7 = pass initiated)
    chunk_start_time = mat_subj_final(start_iter_chunk,col('Time_Stamp'));
        %get the time stamp for when the chunk started
    
    % time of chunk end (key press)
    end_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==10,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 10 which is 'Key Pressed for 
        %clause' (next clause) (8 = passage begins, 7 = pass initiated, 
        %9 = main clause sound file begins playing)
       
    
    for iteration = start_iter_chunk:end_iter_chunk
        
        % make time relative to chunk start for that chunk
        mat_subj_final(iteration,col('Time_Rel_Chunk')) = ...
            mat_subj_final(iteration,col('Time_Stamp')) - chunk_start_time;
        
    end %end iteration (sample) index for time rel to pass start
    
end %end passage index





%% Pupil Pre-Processing


% Sample Messages (for reference):
% 1 = 'MODE RECORD CR 1000'
% 2 = 'Whitescreen60'
% 3 = 'Blackscreen60'
% 4 = 'BaselineStart'
% 5 = 'BaselineEnd'
% 6 = 'Instructions for main displayed'
% 7 = 'PassageInitiated'
% 8 = 'DispCross_PassBegin'
% 9 = 'Main clause sound file begins playing'
% 10 = 'Key Pressed for clause' (next clause)
% 11 = 'DispDot_PassEnd'
% 12 = 'Start Recalling (pupil end)'
% 13 = 'TRACKER_TIME'

trial_pupmean = NaN; %set trial/passage pupil mean to NaN

for pass = 1:no_passages
    
    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==7,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 7 which is passage initiated 
        %(8 = display cross passage begins)
        %NOTE: this is starting in a different spot from time rel pass
        %start
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12=start recalling
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot
        %passage end)
    
    temp_pup_holder = mat_subj_final(start_iter_pass:end_iter_pass,...
        col('Pupil_Size'));
    trial_pupmean = nanmean(temp_pup_holder); %find the mean pupil diam
    trial_pupstd = std(temp_pup_holder, 'omitnan'); %standard dev pupil 

    count_blinklength = 1;
    
    
    for iteration = start_iter_pass:size(mat_subj_final,1)
        
        if isnan(mat_subj_final(iteration,col('Pupil_Size'))) == 0 
            %if pupil size is not NaN

            if mat_subj_final(iteration,col('Pupil_Size'))<(trial_pupmean-...
                    (3*trial_pupstd)) %blink identifier

                mat_subj_final(iteration,col('Pupil_Size')) = NaN; %should 
                %re-assign the value of the current pupil size

                if iteration > 1 %if it is currently first iteration, 
                    %can't look at one before
                    if mat_subj_final(iteration-1,col('Pupil_Size')) < ...
                            (trial_pupmean - (3*trial_pupstd)) 
                        %if previous pupil size is part of blink
                        count_blinklength = count_blinklength + 1;
                    else
                        count_blinklength = 1; %if previous pupil size is 
                        %not part of blink, then counter should be one
                    end
                end

            else

                if count_blinklength >= 2 %if equal to or greater than 2 
                    %(1ms per sample)
                    % would only still be greater than one if it is
                    % just after a blink (because it resets to 1 right
                    % after this block of code)
                    % applying this to any length blink now (Winn 2016)

                    if iteration >= 80+count_blinklength %if we have had 
                        %enough samples prior to this blink

                        for iterate_blink = -(count_blinklength+80):(1+160) 
                            %start at beginning of blink and go to end, 
                            % but starting 80ms before and going
                            % to 160ms after the blink (see
                            % Zekveld, Rudner, et al., 2014)

                            vector = [mat_subj_final(iteration+...
                                iterate_blink-1,col('Pupil_Size')); ...
                                mat_subj_final(iteration,...
                                col('Pupil_Size'))];
                                %vector is made up of the pupil size just
                                %before the blink (or the just-previously-
                                %assigned newly calculated pupil size and 
                                %the current pupil size which is just after
                                %the blink
                            mat_subj_final(iteration+iterate_blink,...
                                col('Pupil_Size')) = nanmean(vector);
                                %assign current pupil size to be the mean
                                %of those two values

                            vector = [NaN;NaN];

                        end %for iterate_blink

                    end %if we had enough samples prior to this blink

                end %if count_blinklength >= 2

                count_blinklength = 1;

            end %blink identifier (if in blink)

        end %if pupil size is not NaN
        
        if (iteration<size(mat_subj_final,1)) && ...
                (mat_subj_final(iteration+1,col('Passage'))~=pass) 
            %if next iter will no longer be in the same passage
            %(but we are not at the end)
            break %break out of for loop to move on to next passage
        end
        
    end %end iteration
    
end %end passage


for pass = 1:no_passages

    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==7,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 7 which is passage initiated 
        %(8 = display cross passage begins)
        %NOTE: this is starting in a different spot from time rel pass
        %start
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12=start recalling
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot
        %passage end)
    
    for iteration = start_iter_pass:size(mat_subj_final,1)

        if mat_subj_final(iteration,col('Pupil_Size')) == 0 
            %if pupil size is 0 
            %(even after previous processing step)

            mat_subj_final(iteration,col('Pupil_Size')) = NaN;
            
        end
        
        if (iteration<size(mat_subj_final,1)) && ...
                (mat_subj_final(iteration+1,col('Passage'))~=pass) 
            %if next iter will no longer be in the same passage
            %(but we are not at the end)
            break %break out of for loop to move on to next passage
        end
    end %for iteration
end %for passage



for pass = 1:no_passages
    
    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==7,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 7 which is passage initiated 
        %(8 = display cross passage begins)
        %NOTE: this is starting in a different spot from time rel pass
        %start
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12=start recalling
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot
        %passage end)
        
    mat_subj_final(start_iter_pass:end_iter_pass,col('Pupil_Size')) = ...
        movmedian(mat_subj_final(start_iter_pass:end_iter_pass,...
        col('Pupil_Size')),250,'omitnan');
    %creates a moving MEDIAN (to help with outliers) of vector, length 250, 
    %centered, omitting nans

end


%% Scale Pupil Sizes

% % new columns have been pre-alloc in original table (so matrix has names
% mapped on already)
% tbl_subj_final{:,'Baseline'} = 0;
% tbl_subj_final{:,'Pupil_Size_ScBa'} = NaN;
% tbl_subj_final{:,'Pre_Recall'} = 0;

% pre-alloc matrix for baselines
baselineavg_perpass = NaN(no_passages,3);
    %1 = subj no
    %2 = passage number
    %3 = baseline average



for pass = 1:no_passages
    
    pass_is = pass; %this seems to sometimes be necessary for some reason

    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==8,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 8 which is display cross
        %passage begins (7 = passage initiated)
        
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12=start recalling 
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot
        %passage end)
        
    % get baseline start and end indices for current passage
    base_start_ind = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==4,1,'first');
        % 4 = BaselineStart
    base_end_ind = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==5,1,'first');
        % 5 = BaselineEnd
        
    % get pre-recall start and end indices for current passage
    prerecall_start_ind = find(mat_subj_final(:,col('Passage'))==pass &...
        mat_subj_final(:,col('Sample_Message'))==11,1,'first'); 
        %should be scalar
        % 11 = DispDot_PassEnd
    prerecall_end_ind = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first');
        % 12 = Start Recalling (pupil end)
        
        
    % find average for baseline period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get just baseline pupil info for current passage
    cur_baseline_pup = mat_subj_final(base_start_ind:base_end_ind,...
        col('Pupil_Size'));
    % get the mean of that
    cur_baseline_avg = nanmean(cur_baseline_pup); %find the mean pupil diam
    % put in avg base pupil
    baselineavg_perpass(pass_is,3) = cur_baseline_avg; %put avg base pupil
    % put in the passage number
    baselineavg_perpass(pass_is,2) = pass_is; %indicate the passage number
    % put in the subject number
    baselineavg_perpass(pass_is,1) = subject_no; %indicate the subject number
        
    
    % indicate yes(1)/no(0) in new col whether in a baseline period %%%%%%%
    
    for iteration = base_start_ind:base_end_ind
        
        %indicates yes/no baseline period
        mat_subj_final(iteration,col('Baseline')) = 1;
        

    end % end iterate through baseline
    
    
    % scale pupil size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iteration = start_iter_pass:end_iter_pass
    
        mat_subj_final(iteration,col('Pupil_Size_ScBa')) = ...
            (mat_subj_final(iteration,col('Pupil_Size')) -...
            baselineavg_perpass(pass_is,3)) / pupil_dr;
        
    end %end iterate through passage
        
        
    
end %passage index



%% Pupil in Bins through Time

% Sample Messages (for reference):
% 1 = 'MODE RECORD CR 1000'
% 2 = 'Whitescreen60'
% 3 = 'Blackscreen60'
% 4 = 'BaselineStart'
% 5 = 'BaselineEnd'
% 6 = 'Instructions for main displayed'
% 7 = 'PassageInitiated'
% 8 = 'DispCross_PassBegin'
% 9 = 'Main clause sound file begins playing'
% 10 = 'Key Pressed for clause' (next clause)
% 11 = 'DispDot_PassEnd'
% 12 = 'Start Recalling (pupil end)'
% 13 = 'TRACKER_TIME'

%%% to subset columns wanted from matrix:
% subData = data(:, [col('A') col('C')]);

% create new tables and matrices, since doing bins will re-size the length
mat_subj_bins_pass_start = NaN((10000*no_passages),13); 
%go from passage start until recall
    %1: subj no
    %2: trial index (order played? minus 1 for b/w scrs)
    %3: passage
    %4: HLExpNum
    %5: Prob
    %6: EorNNum
    %7: Processing
    %8: Pupil dr
    %9: Pupil base
    %10: Bin number
    %11: Time stamp at start of bin
    %12: Pupil size (raw)
    %13: Pupil size scaled/baselined
var_names_bins_pass_start = {'Subject_No','Trial_Index','Passage',...
    'HLExpNum','Prob','EorNNum','Processing','Pupil_DR','Pupil_Base',...
    'Bin_No','Time_Stamp_BinStart','Pupil_Size','Pupil_Size_ScBa'};
col_bins_pass_start = containers.Map(var_names_bins_pass_start, ...
    1:numel(var_names_bins_pass_start));

mat_subj_bins_pass_end = NaN((13*1000)*no_passages,13); 
%from 10sec before to 3sec after
%(recall_time_ms = 3000ms)
    %1: subj no
    %2: trial index (order played? minus 1 for b/w scrs)
    %3: passage
    %4: HLExpNum
    %5: Prob
    %6: EorNNum
    %7: Processing
    %8: Pupil dr
    %9: Pupil base
    %10: Bin number
    %11: Time stamp at start of bin
    %12: Pupil size (raw)
    %13: Pupil size scaled/baselined
    %14: Pupil size scaled/baselined to end of passage
var_names_bins_pass_end = {'Subject_No','Trial_Index','Passage',...
    'HLExpNum','Prob','EorNNum','Processing','Pupil_DR','Pupil_Base',...
    'Bin_No','Time_Stamp_BinStart','Pupil_Size','Pupil_Size_ScBa',...
    'Pupil_Size_ScBaE'};
col_bins_pass_end = containers.Map(var_names_bins_pass_end, ...
    1:numel(var_names_bins_pass_end));

mat_subj_bins_chunk_start = NaN(1000*no_total_chunks,17); 
%go from chunk start until key press
    %1: subj no
    %2: trial index (order played? minus 1 for b/w scrs)
    %3: passage
    %4: HLExpNum
    %5: Prob
    %6: EorNNum
    %7: Processing
    %8: ChunkNo_Overall (chunk # in overall expt, not just within passage)
    %9: chunk
    %10: SorCNum
    %11: Pupil dr
    %12: Pupil base
    %13: Bin number
    %14: Time stamp at start of bin
    %15: Pupil size (raw)
    %16: Pupil size scaled/baselined
    %17: Pupil size scaled/baselined to chunk start
var_names_bins_chunk_start = {'Subject_No','Trial_Index','Passage',...
    'HLExpNum','Prob','EorNNum','Processing','ChunkNo_Overall','Chunk',...
    'SorCNum','Pupil_DR','Pupil_Base',...
    'Bin_No','Time_Stamp_BinStart','Pupil_Size','Pupil_Size_ScBa',...
    'Pupil_Size_ScBaC'};
col_bins_chunk_start = containers.Map(var_names_bins_chunk_start, ...
    1:numel(var_names_bins_chunk_start));

% % columns for bin numbers have been pre-alloc:
% tbl_subj_final{:,'Bins_PassStart'} = NaN;
% tbl_subj_final{:,'Bins_ChunkStart'} = NaN;
% tbl_subj_final{:,'Bins_PassEnd'} = NaN;

bin_size = 50; %50ms or 50 sample bins


%%% For Bins Relative to Passage Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count_bin = 0; %initialize a counter for which bin we are on,
    %OUTSIDE of loop, since this is long form
    %(this way each passage does not re-start at the top of the matrix)

for pass = 1:no_passages
    
    % iteration at passage start 
    start_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==8,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 8 which is display cross
        %passage begins (7 = passage initiated)
        
    % iteration at passage end 
    end_iter_pass = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12 which is start recalling
        %(pupil ends) (7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot passage
        %ends
        
    count_bin_inside = 0; %initialize a counter for which bin we are on
        %INSIDE the loop, to bin count within a passage 
    
    time_stamp_binstart = NaN; %initialize time stamp at bin start, in case
    
    for bin_iter = (start_iter_pass+50):50:end_iter_pass
            %note: any remainder samples that do not add up to 50 at the
            %end will not be included (e.g. 5:5:16 = [5 10 15] with the
            %remaining 1 not included
            % Starting at start_iter_pass + 50 since we then back up 50
            % samples with for sample = -49:0
            
        temp_raw_pup_bin = NaN(bin_size,1);
        temp_sc_pup_bin = NaN(bin_size,1);
        
        count_bin = count_bin + 1; %keep track of bin number
        count_bin_inside = count_bin_inside +1; %bin number within passage
        
        count_sample = 0; %initialize counter for which sample 
            %within bin we are on
            
        time_stamp_binstart = mat_subj_final(bin_iter-49,...
            col('Time_Stamp')); %get the time stamp at the start of the 
            %current bin (49 samples ago --> see for sample below)
        
        for sample = -49:0
            
            if isnan(mat_subj_final(bin_iter+sample,col('Pupil_Size')))...
                    ==0 && mat_subj_final(bin_iter+sample,...
                    col('Pupil_Size')) ~= 0
                %if pupil size is NOT NaN and NOT Zero (should have been
                %de-blinked already though)
                
                count_sample = count_sample + 1; %keep track of sample #
                temp_raw_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size'));
                temp_sc_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size_ScBa'));
                
            end %end if pupil size exists
            
        end %end for sample -49:0
        
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Subject_No')) = subject_no; %put subj no in col 1
        for i = col_bins_pass_start('Trial_Index'):...
                col_bins_pass_start('Processing')
            %want to loop through from trial index to processing
            mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
                (var_names_bins_pass_start{i})) = ...
                mat_subj_final(bin_iter+sample,...
                col(var_names_bins_pass_start{i}));
        end
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Pupil_DR')) = pupil_dr;
            %put in pupillary dynamic range for this participant
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Pupil_Base')) = baselineavg_perpass(pass,3);
            %put in baseline average (col 3) for this passage
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Bin_No')) = count_bin_inside; %put in bin number 
            %INSIDE the passage bin count
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Time_Stamp_BinStart')) = time_stamp_binstart; 
            %time stamp at bin start
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Pupil_Size')) = nanmean(temp_raw_pup_bin);
        mat_subj_bins_pass_start(count_bin,col_bins_pass_start...
            ('Pupil_Size_ScBa')) = nanmean(temp_sc_pup_bin);
        
    end %end for bin (jumping by 50's)
    
end %end for pass


%%% For Bins Relative to Passage End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start 10 sec (10,000 ms) before passage ends, End 3 sec (3,000 ms) later
% when Recall period starts

count_bin = 0; %initialize a counter for which bin we are on,
    %OUTSIDE of loop, since this is long form
    %(this way each passage does not re-start at the top of the matrix)
    
baselines_passend = NaN(no_passages,3);

for pass = 1:no_passages
    
    % iteration at passage 10 sec before passage end
    pass_snd_end_iter = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==11,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 11: display dot passage end 
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 12=start recalling
        %(pupil ends) )
    pass_snd_10msbefend_iter = pass_snd_end_iter - 10000; %10sec before
        
    % iteration at start of recall period 
    recall_starts_iter = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12 which is start recalling
        %(pupil ends) (7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot passage
        %ends
        
    count_bin_inside = 0; %initialize a counter for which bin we are on 
        %for INSIDE the passage
        
    % brief baseline at passage end
    baseline_passend_cur = mat_subj_final(pass_snd_end_iter,col('Pupil_Size'));
    baselines_passend(pass,3) = baseline_passend_cur;
    baselines_passend(pass,1) = subject_no;
    baselines_passend(pass,2) = pass;
    
    time_stamp_binstart = NaN; %initialize time stamp at bin start, in case
    
    for bin_iter = pass_snd_10msbefend_iter:50:recall_starts_iter
            %note: any remainder samples that do not add up to 50 at the
            %end will not be included (e.g. 5:5:16 = [5 10 15] with the
            %remaining 1 not included
            
        temp_raw_pup_bin = NaN(bin_size,1);
        temp_sc_pup_bin = NaN(bin_size,1);
        
        count_bin = count_bin + 1; %keep track of bin number
        count_bin_inside = count_bin_inside +1; %bin count INSIDE passage
        
        count_sample = 0; %initialize counter for which sample 
            %within bin we are on
            
        time_stamp_binstart = mat_subj_final(bin_iter-49,...
            col('Time_Stamp')); %get the time stamp at the start of the 
            %current bin (49 samples ago --> see for sample below)
        
        for sample = -49:0
            
            if isnan(mat_subj_final(bin_iter+sample,col('Pupil_Size')))...
                    ==0 && mat_subj_final(bin_iter+sample,...
                    col('Pupil_Size')) ~= 0
                %if pupil size is NOT NaN and NOT Zero (should have been
                %de-blinked already though)
                
                count_sample = count_sample + 1; %keep track of sample #
                temp_raw_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size'));
                temp_sc_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size_ScBa'));
                
            end %end if pupil size exists
            
        end %end for sample -49:0
        
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Subject_No')) = subject_no; %put subj no in col 1
        for i = col_bins_pass_end('Trial_Index'):...
                col_bins_pass_end('Processing')
            %want to loop through from trial index to processing
            mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
                (var_names_bins_pass_end{i})) = ...
                mat_subj_final(bin_iter+sample,...
                col(var_names_bins_pass_end{i}));
        end
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Pupil_DR')) = pupil_dr;
            %put in pupillary dynamic range for this participant
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Pupil_Base')) = baselineavg_perpass(pass,3);
            %put in baseline average (col 3) for this passage
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Bin_No')) = count_bin_inside; %put in bin number
            %INSIDE the passage bin counter
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Time_Stamp_BinStart')) = time_stamp_binstart; 
            %time stamp at bin start
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Pupil_Size')) = nanmean(temp_raw_pup_bin);
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Pupil_Size_ScBa')) = nanmean(temp_sc_pup_bin);
        mat_subj_bins_pass_end(count_bin,col_bins_pass_end...
            ('Pupil_Size_ScBaE')) = (nanmean(temp_raw_pup_bin)-...
            baseline_passend_cur)/pupil_dr;
        
    end %end for bin (jumping by 50's)
    
end %end for pass


%%% For Bins Relative to Chunk Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count_bin = 0; %initialize a counter for which bin we are on,
    %OUTSIDE of loop, since this is long form
    %(this way each passage does not re-start at the top of the matrix)
    
baselines_chunks = NaN(no_total_chunks,3);

for chunk = 1:no_total_chunks
    
    % iteration at chunk start 
    start_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==9,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 9 which is 'Main clause sound 
        %file begins playing' (8 = passage begins, 7 = pass initiated)
        
    % iteration at chunk end (key pressed) 
    end_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==10,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 10 which is 
        %'Key Pressed for clause' (next clause)
        
    % brief baseline at chunk start
    baseline_chunkstart_arr = mat_subj_final((start_iter_chunk-24):...
        start_iter_chunk,col('Pupil_Size'));
    baseline_chunkstart = nanmean(baseline_chunkstart_arr);
    baselines_chunks(chunk,3) = baseline_chunkstart;
    baselines_chunks(chunk,1) = subject_no;
    baselines_chunks(chunk,2) = chunk;
        
    count_bin_inside = 0; %initialize a counter for which bin we are on
        %INSIDE the chunk
    
    time_stamp_binstart = NaN; %initialize time stamp at bin start, in case
    
    for bin_iter = start_iter_chunk:50:end_iter_chunk
            %note: any remainder samples that do not add up to 50 at the
            %end will not be included (e.g. 5:5:16 = [5 10 15] with the
            %remaining 1 not included
            
        temp_raw_pup_bin = NaN(bin_size,1);
        temp_sc_pup_bin = NaN(bin_size,1);
        
        count_bin = count_bin + 1; %keep track of bin number
        count_bin_inside = count_bin_inside +1; 
            %keep track of INSIDE of chunk bin number
        
        count_sample = 0; %initialize counter for which sample 
            %within bin we are on
            
        time_stamp_binstart = mat_subj_final(bin_iter-49,...
            col('Time_Stamp')); %get the time stamp at the start of the 
            %current bin (49 samples ago --> see for sample below)
        
        for sample = -49:0
            
            if isnan(mat_subj_final(bin_iter+sample,col('Pupil_Size')))...
                    ==0 && mat_subj_final(bin_iter+sample,...
                    col('Pupil_Size')) ~= 0
                %if pupil size is NOT NaN and NOT Zero (should have been
                %de-blinked already though)
                
                count_sample = count_sample + 1; %keep track of sample #
                temp_raw_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size'));
                temp_sc_pup_bin(count_sample,1) = mat_subj_final...
                    (bin_iter+sample,col('Pupil_Size_ScBa'));
                
            end %end if pupil size exists
            
        end %end for sample -49:0

        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Subject_No')) = subject_no; %put subj no in col 1
        for i = col_bins_chunk_start('Trial_Index'):...
                col_bins_chunk_start('SorCNum')
            %want to loop through from trial index to SorCNum
            mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
                (var_names_bins_chunk_start{i})) = ...
                mat_subj_final(bin_iter+sample,...
                col(var_names_bins_chunk_start{i}));
        end
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Pupil_DR')) = pupil_dr;
            %put in pupillary dynamic range for this participant
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Pupil_Base')) = baselineavg_perpass(pass,3);
            %put in baseline average (col 3) for this passage
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Bin_No')) = count_bin_inside; %put in bin number
            %INSIDE the chunk
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Time_Stamp_BinStart')) = time_stamp_binstart; 
            %time stamp at bin start
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Pupil_Size')) = nanmean(temp_raw_pup_bin);
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Pupil_Size_ScBa')) = nanmean(temp_sc_pup_bin);
        mat_subj_bins_chunk_start(count_bin,col_bins_chunk_start...
            ('Pupil_Size_ScBaC')) = (nanmean(temp_raw_pup_bin)-...
            baseline_chunkstart)/pupil_dr;
        
    end %end for bin (jumping by 50's)
    
end %end for chunk



%% Peak Pupils

% Sample Messages (for reference):
% 1 = 'MODE RECORD CR 1000'
% 2 = 'Whitescreen60'
% 3 = 'Blackscreen60'
% 4 = 'BaselineStart'
% 5 = 'BaselineEnd'
% 6 = 'Instructions for main displayed'
% 7 = 'PassageInitiated'
% 8 = 'DispCross_PassBegin'
% 9 = 'Main clause sound file begins playing'
% 10 = 'Key Pressed for clause' (next clause)
% 11 = 'DispDot_PassEnd'
% 12 = 'Start Recalling (pupil end)'
% 13 = 'TRACKER_TIME'

%%% to subset columns wanted from matrix:
% subData = data(:, [col('A') col('C')]);

% number of samples we want the peak to be to be stable
no_samples_peak = 10;


%%% For Peak at Passage End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create new tables and matrices
mat_subj_peak_passend = NaN(no_passages,13); 
%go from passage end until recall
    %1: subj no
    %2: trial index (order played? minus 1 for b/w scrs)
    %3: passage
    %4: HLExpNum
    %5: Prob
    %6: EorNNum
    %7: Processing
    %8: Pupil dr
    %9: Pupil base
    %10: Pupil base passage end
    %11: Time stamp at end of pass
    %12: Pupil size (raw)
    %13: Pupil size scaled/baselined to passage beginning
    %14: Pupil size scaled/baselined to passage end
var_names_peak_passend = {'Subject_No','Trial_Index','Passage',...
    'HLExpNum','Prob','EorNNum','Processing','Pupil_DR','Pupil_Base',...
    'Pupil_Base_PassEnd','Time_Stamp_PassEnd','Pupil_Size',...
    'Pupil_Size_ScBa','Pupil_Size_ScBaE'};
col_peak_passend = containers.Map(var_names_peak_passend, ...
    1:numel(var_names_peak_passend));
    
for pass = 1:no_passages
     
    % iteration at passage end
    pass_snd_end_iter = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==11,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 11: display dot passage end 
        %(7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 12=start recalling
        %(pupil ends) )
    
    % iteration at start of recall period 
    recall_starts_iter = find(mat_subj_final(:,col('Passage'))==pass & ...
        mat_subj_final(:,col('Sample_Message'))==12,1,'first'); 
        %finds the first 1 index that equals trial index in the passage
        %column and where sample message == 12 which is start recalling
        %(pupil ends) (7=passage initiated, 8=passage begins, 9=sound file 
        %begins playing, 10=key pressed for clause, 11=display dot passage
        %ends

    
    time_stamp_passend = mat_subj_final(pass_snd_end_iter,...
        col('Time_Stamp')); 
    
    % Get a temporary array of the pupil sizes within that window
    % (pre-recall period)
    temp_pupil_raw = mat_subj_final(pass_snd_end_iter:...
        (recall_starts_iter-10),col('Pupil_Size'));
    % grab the maximum pupil size and its index (note that already smoothed and de-blinked)
    [max_pup,index_max_pup] = max(temp_pupil_raw,[],'omitnan');
    % make a small array of 5 samples before to 5 samples after that max,
    % to create a more stable peak:
        % NOTE: modified so that it can account for if max is last sample
        % NOTE 2: it automatically takes first index if multiple maxes
    size_temppup = size(temp_pupil_raw,1); %size of temp pup array
    diff = size_temppup - index_max_pup;
    if index_max_pup <= no_samples_peak %if peak is within first 10 samples
        peak_pupil_temparr = temp_pupil_raw(1:no_samples_peak);
            % if peak is within 10 samples of the beginning,
            % just take the 10 samples at the beginning
    elseif diff > no_samples_peak %if the difference is greater than 
            %number of samples for peak
        peak_pupil_temparr = temp_pupil_raw((index_max_pup-5):...
            (index_max_pup+5));
    else
        peak_pupil_temparr = temp_pupil_raw(end-no_samples_peak:end);
            % if peak is within 10 samples of the end, just take the 10
            % samples at the end
    end
    % calculate the mean of that 10-sample peak
    peak_pupil = mean(peak_pupil_temparr,'omitnan');
        
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Subject_No')) = subject_no; %put subj no in col 1
    for i = col_peak_passend('Trial_Index'):...
            col_peak_passend('Processing')
        %want to loop through from trial index to processing
        mat_subj_peak_passend(pass,col_peak_passend...
            (var_names_peak_passend{i})) = ...
            mat_subj_final(pass_snd_end_iter,...
            col(var_names_peak_passend{i}));
    end
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_DR')) = pupil_dr;
        %put in pupillary dynamic range for this participant
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_Base')) = baselineavg_perpass(pass,3);
        %put in baseline average (col 3) for this passage
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_Base_PassEnd')) = baselines_passend(pass,3);
        %put in baseline average for end of this passage
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Time_Stamp_PassEnd')) = time_stamp_passend; 
        %time stamp at passage end
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_Size')) = nanmean(temp_raw_pup_bin);
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_Size_ScBa')) = (peak_pupil-baselineavg_perpass(pass,3))/...
        pupil_dr;
    mat_subj_peak_passend(pass,col_peak_passend...
        ('Pupil_Size_ScBaE')) = (peak_pupil-...
        baselines_passend(pass,3))/pupil_dr;
        
    
end %end for pass


%%% For Peak at Chunk End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create new tables and matrices
mat_subj_peak_chunkend = NaN(no_total_chunks,17); 
%go from chunk start until key press
    %1: subj no
    %2: trial index (order played? minus 1 for b/w scrs)
    %3: passage
    %4: HLExpNum
    %5: Prob
    %6: EorNNum
    %7: Processing
    %8: ChunkNo_Overall (chunk # in overall expt, not just within passage)
    %9: chunk
    %10: SorCNum
    %11: Pupil dr
    %12: Pupil base
    %13: Pupil base at chunk start
    %14: Time stamp at start of chunk
    %15: Pupil size (raw)
    %16: Pupil size scaled/baselined
    %17: Pupil size scaled/baselined to chunk start
var_names_peak_chunkend = {'Subject_No','Trial_Index','Passage',...
    'HLExpNum','Prob','EorNNum','Processing','ChunkNo_Overall','Chunk',...
    'SorCNum','Pupil_DR','Pupil_Base','Pupil_Base_ChunkStart',...
    'Time_Stamp_ChunkStart','Pupil_Size','Pupil_Size_ScBa',...
    'Pupil_Size_ScBaC'};
col_peak_chunkend = containers.Map(var_names_peak_chunkend, ...
    1:numel(var_names_peak_chunkend));
    
for chunk = 1:no_total_chunks
     
    % iteration at chunk start 
    start_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==9,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 9 which is 'Main clause sound 
        %file begins playing' (8 = passage begins, 7 = pass initiated)
        
    % iteration at chunk end (key pressed) 
    end_iter_chunk = find(mat_subj_final(:,col('ChunkNo_Overall'))==...
        chunk & mat_subj_final(:,col('Sample_Message'))==10,1,'first'); 
        %finds the first 1 index that equals chunk index in the overallno
        %column and where sample message == 10 which is 
        %'Key Pressed for clause' (next clause)
    
    time_stamp_chunkstart = mat_subj_final(start_iter_chunk,...
        col('Time_Stamp')); 
    
    % Get a temporary array of the pupil sizes within that window
    % (chunk until 10 samples before press):
    temp_pupil_raw = mat_subj_final(start_iter_chunk:...
        (end_iter_chunk-10),col('Pupil_Size'));
    % grab the maximum pupil size and its index:
    [max_pup,index_max_pup] = max(temp_pupil_raw,[],'omitnan');
    % make a small array of 5 samples before to 5 samples after that max,
    % to create a more stable peak:
        % NOTE: modified so that it can account for if max is last sample
        % NOTE 2: it automatically takes first index if multiple maxes
    size_temppup = size(temp_pupil_raw,1); %size of temp pup array
    diff = size_temppup - index_max_pup;
    if index_max_pup <= no_samples_peak %if peak is within first 10 samples
        peak_pupil_temparr = temp_pupil_raw(1:no_samples_peak);
            % if peak is within 10 samples of the beginning,
            % just take the 10 samples at the beginning
    elseif diff > no_samples_peak %if the difference is greater than 
            %number of samples for peak
        peak_pupil_temparr = temp_pupil_raw((index_max_pup-5):...
            (index_max_pup+5));
    else
        peak_pupil_temparr = temp_pupil_raw(end-no_samples_peak:end);
            % if peak is within 10 samples of the end, just take the 10
            % samples at the end
    end
    % calculate the mean of that 10-sample peak:
    peak_pupil = mean(peak_pupil_temparr,'omitnan');
        
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Subject_No')) = subject_no; %put subj no in col 1
    for i = col_peak_chunkend('Trial_Index'):...
            col_peak_chunkend('SorCNum')
        %want to loop through from trial index to processing
        mat_subj_peak_chunkend(chunk,col_peak_chunkend...
            (var_names_peak_chunkend{i})) = ...
            mat_subj_final(start_iter_chunk,...
            col(var_names_peak_chunkend{i}));
    end
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_DR')) = pupil_dr;
        %put in pupillary dynamic range for this participant
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_Base')) = baselineavg_perpass(pass,3);
        %put in baseline average (col 3) for this passage
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_Base_ChunkStart')) = baselines_chunks(chunk,3);
        %put in baseline average for end of this passage
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Time_Stamp_ChunkStart')) = time_stamp_chunkstart; 
        %time stamp at passage end
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_Size')) = nanmean(temp_raw_pup_bin);
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_Size_ScBa')) = (peak_pupil-baselineavg_perpass(pass,3))/...
        pupil_dr;
    mat_subj_peak_chunkend(chunk,col_peak_chunkend...
        ('Pupil_Size_ScBaC')) = (peak_pupil-...
        baselines_chunks(chunk,3))/pupil_dr;
        
    
end %end for pass




%% Vertically cocatenate

key_lat_allpart = vertcat(key_lat_allpart,key_lat);
mat_final_allpart = vertcat(mat_final_allpart,mat_subj_final);
prerecall_mat_allpart = vertcat(prerecall_mat_allpart,prerecall_mat);
baselineavg_perpass_allpart = vertcat(baselineavg_perpass_allpart,...
    baselineavg_perpass);
baselines_chunks_allpart = vertcat(baselines_chunks_allpart,...
    baselines_chunks);
baselines_passend_allpart = vertcat(baselines_passend_allpart,...
    baselines_passend);
pupil_wbinfo_allpart = vertcat(pupil_wbinfo_allpart,pupil_wbinfo);
mat_subj_bins_pass_start_allpart = vertcat(...
    mat_subj_bins_pass_start_allpart,mat_subj_bins_pass_start);
mat_subj_bins_pass_end_allpart = vertcat(mat_subj_bins_pass_end_allpart,...
    mat_subj_bins_pass_end);
mat_subj_bins_chunk_start_allpart = vertcat(...
    mat_subj_bins_chunk_start_allpart,mat_subj_bins_chunk_start);
mat_subj_peak_passend_allpart = vertcat(mat_subj_peak_passend_allpart,...
    mat_subj_peak_passend);
mat_subj_peak_chunkend_allpart = vertcat(mat_subj_peak_chunkend_allpart,...
    mat_subj_peak_chunkend);



end % end if participant is Include~=0 (first column participant sheet)

end % end participant for loop

% convert back to table so we have the variable names
tbl_final_allpart = array2table(mat_final_allpart,'VariableNames',...
    varNames);
prerecall_tbl_allpart = array2table(prerecall_mat_allpart,'VariableNames',...
    varNames);
tbl_subj_bins_pass_start_allpart = array2table(...
    mat_subj_bins_pass_start_allpart,'VariableNames',...
    var_names_bins_pass_start);
tbl_subj_bins_pass_end_allpart = array2table(...
    mat_subj_bins_pass_end_allpart,'VariableNames',...
    var_names_bins_pass_end);
tbl_subj_bins_chunk_start_allpart = array2table(...
    mat_subj_bins_chunk_start_allpart,'VariableNames',...
    var_names_bins_chunk_start);
tbl_subj_peak_passend_allpart = array2table(...
    mat_subj_peak_passend_allpart,'VariableNames',var_names_peak_passend);
tbl_subj_peak_chunkend_allpart = array2table(...
    mat_subj_peak_chunkend_allpart,'VariableNames',var_names_peak_chunkend);



%% Write info to csvs etc.

csvwrite('PROJ_key_lat_allpart.csv',key_lat_allpart);

csvwrite('PROJ_baselineavg_perpass_allpart.csv',...
    baselineavg_perpass_allpart);
csvwrite('PROJ_baselines_chunks_allpart.csv',...
    baselines_chunks_allpart);
csvwrite('PROJ_baselines_passend_allpart.csv',...
    baselines_passend_allpart);

csvwrite('PROJ_pupil_wbinfo_allpart.csv',...
    pupil_wbinfo_allpart);

writetable(tbl_final_allpart,'PROJ_tbl_sample_final_allpart.csv',...
    'Delimiter',',');

writetable(prerecall_tbl_allpart,'PROJ_prerecall_tbl_allpart.csv',...
    'Delimiter',',');

writetable(tbl_subj_bins_pass_start_allpart,...
    'PROJ_tbl_subj_bins_pass_start_allpart.csv','Delimiter',',');
writetable(tbl_subj_bins_pass_end_allpart,...
    'PROJ_tbl_subj_bins_pass_end_allpart.csv','Delimiter',',');
writetable(tbl_subj_bins_chunk_start_allpart,...
    'PROJ_tbl_subj_bins_chunk_start_allpart.csv','Delimiter',',');

writetable(tbl_subj_peak_passend_allpart,...
    'PROJ_tbl_subj_peak_passend_allpart.csv','Delimiter',',');
writetable(tbl_subj_peak_chunkend_allpart,...
    'PROJ_tbl_subj_peak_chunkend_allpart.csv','Delimiter',',');




% change working directory back
cd('/Documents/MATLAB');



