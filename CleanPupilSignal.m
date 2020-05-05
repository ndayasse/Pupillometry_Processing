function [output_vector] = CleanPupilSignal(vector,sample_rate,...
    method,moving_window_fntn,remove_prop)
%CleanPupilSignal pre-processes and cleans the time-course pupillary signal
%   Takes a numeric vector of pupil samples
%   Removes blinks and performs spline or linear interpolation
%   Applies a smoothing function
%   Note: sample_rate should be given in Hz
%   Method must be 'spline' for spline interpolation or 'linear' for linear
%       interpolation for imputing the missing values
%   moving_window_fntn is  'median' for using movmedian function, or
%       smoothing data over a 250ms window using median function that is
%       less affected by outliers, can also set to 'mean' to use movmean
%       instead to us mean function
%   remove_prop is optional input to remove trials (make all NaN) that have
%       more than the given proportion missing


%% Initial checks:
% check vector is appropriately set up:
size_vec = size(vector);
if size_vec(1)==1 && size_vec(2)==1
    error('vector must be a one-dimensional vector longer than 1')
elseif size_vec(1) > 1 && size_vec(2) > 1
    error('vector must be a one-dimensional vector')
elseif length(size_vec) > 2
    error('vector must be a one-dimensional vector')
end

% check number of input arguments:
switch nargin
    
    case 5 %all inputs used including remove_prop
        
        %% Calculate pupil mean and standard deviation (will need for identifying
        % blinks):
        pupil_mean = nanmean(vector);
        pupil_std = std(vector,'omitnan');

        % Go through and identify blinks and remove them
        count_blink_samples = 1;
        for iteration = 1:size(vector,1)

            if ~isnan(vector(iteration)) %if pupil size is not NaN

                if vector(iteration) < (pupil_mean-(3*pupil_std)) 
                        %blink identifier

                    vector(iteration) = NaN; 
                        %should re-assign the value of the current pupil size
                    count_blink_samples = count_blink_samples + 1;
                        %count a blink sample

                elseif vector(iteration) > (pupil_mean+(3*pupil_std)) 
                    %suggests an odd event occurred such as head moving away

                    vector(iteration) = NaN; 
                        %should re-assign the value of the current pupil size
                    count_blink_samples = count_blink_samples + 1;
                        %count a blink sample

                end

            end

        end %end going through vector looking for blinks to remove

        %% Count percentage of vector that is missing/blink:
        prop_missing = count_blink_samples / length(vector);
        if prop_missing > 0.2 %if more than 20% is missing
            warning('>20% of your pupil data is missing, consider removing trial')
        end
        if prop_missing > remove_prop
            output_vector = NaN(size(vector));
        else

            %% Use interpolation to replace missing values:
            if method=="spline"
                vector_splinefillmiss = fillmissing(vector,'spline');
            elseif method=="linear"
                vector_splinefillmiss = fillmissing(vector,'linear');
            else
                error('method must be given as spline or linear')
            end

            %% Smooth using moving window:

            % calculate how many samples ~250ms to smooth over is based on sampling
            % rate (sampling rate is samples per sec):
            ms_per_sample = 1000 / sample_rate; %calc how many ms per sample
                %(since 1000ms = 1 sec and sample_rate is samples per sec)
            samples_per_250ms = ceil(250 / ms_per_sample); %calc how many samples 
                %make up 250ms and round up
            % perform smoothing function:
            if moving_window_fntn=="median"
                output_vector = movmedian(vector_splinefillmiss,...
                    samples_per_250ms,'omitnan');
            elseif moving_window_fntn=="mean"
                output_vector = movmean(vector_splinefillmiss,...
                    samples_per_250ms,'omitnan');
            end
        
        end %end if proportion missing is less than the threshold 
        
    case 4 %remove_prop missing/not provided
        
        %% Calculate pupil mean and standard deviation (will need for identifying
        % blinks):
        pupil_mean = nanmean(vector);
        pupil_std = std(vector,'omitnan');

        % Go through and identify blinks and remove them
        count_blink_samples = 1;
        for iteration = 1:size(vector,1)

            if ~isnan(vector(iteration)) %if pupil size is not NaN

                if vector(iteration) < (pupil_mean-(3*pupil_std)) 
                        %blink identifier

                    vector(iteration) = NaN; 
                        %should re-assign the value of the current pupil size
                    count_blink_samples = count_blink_samples + 1;
                        %count a blink sample

                elseif vector(iteration) > (pupil_mean+(3*pupil_std)) 
                    %suggests an odd event occurred such as head moving away

                    vector(iteration) = NaN; 
                        %should re-assign the value of the current pupil size
                    count_blink_samples = count_blink_samples + 1;
                        %count a blink sample

                end

            end

        end %end going through vector looking for blinks to remove

        %% Count percentage of vector that is missing/blink:
        prop_missing = count_blink_samples / length(vector);
        if prop_missing > 0.2 %if more than 20% is missing
            warning('>20% of your pupil data is missing, consider removing trial')
        end

        %% Use interpolation to replace missing values:
        if method=="spline"
            vector_splinefillmiss = fillmissing(vector,'spline');
        elseif method=="linear"
            vector_splinefillmiss = fillmissing(vector,'linear');
        else
            error('method must be given as spline or linear')
        end

        %% Smooth using moving window:

        % calculate how many samples ~250ms to smooth over is based on sampling
        % rate (sampling rate is samples per sec):
        ms_per_sample = 1000 / sample_rate; %calc how many ms per sample
            %(since 1000ms = 1 sec and sample_rate is samples per sec)
        samples_per_250ms = ceil(250 / ms_per_sample); %calc how many samples 
            %make up 250ms and round up
        % perform smoothing function:
        if moving_window_fntn=="median"
            output_vector = movmedian(vector_splinefillmiss,...
                samples_per_250ms,'omitnan');
        elseif moving_window_fntn=="mean"
            output_vector = movmean(vector_splinefillmiss,...
                samples_per_250ms,'omitnan');
        end
        
end %end switch based on nargin


end % end function

