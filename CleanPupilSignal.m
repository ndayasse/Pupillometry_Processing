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
        
        %% De-blink:
        
        % Calculate pupil mean and standard deviation (will need for 
        % identifying blinks):
        pupil_mean = nanmean(vector);
        pupil_std = std(vector,'omitnan');

        % Go through and identify blinks and remove them (make NaN)
        indices = find(vector(:,1)<(pupil_mean-(3*pupil_std)));
            %save the index numbers of where it is a blink
        blink_samples = sum(vector(:,1)<(pupil_mean-(3*pupil_std)));
            %count number of samples that are blinks
        vector(vector(:,1)<(pupil_mean-(3*pupil_std)),:) = NaN;
            %turn the blinks into NaN's
        D = (diff([0,diff(indices')==1,0]))'; 
            %find differences in indices to find consecutive blinked
            %samples --> because need to take out 80ms before to 160ms
            %after a blink (see Zekveld, Rudner, et al., 2014)
        first = indices(D>0); %first indices of consecutive blinks
        last = indices(D<0); %last indices of consecutive blinks
        % Deal with those blinks that occur within the first 80ms or the
        % last 160ms:
        ind_less_80 = indices(:,1)<=80; %logical array of any indices <=80
        if any(ind_less_80) %if any indices are within the first 80
            vector(1:80,1) = NaN; %then make all the first 80 NaN
        end
        ind_last_160 = indices(:,1)>=(size(vector,1)-160); 
            %logical array of any indices within last 160
        if any(ind_last_160) %if any indices are within last 160
            vector((size(vector,1)-160):size(vector,1)) = NaN;
                %then make all the last 160 NaN
        end
        % now take out the 80ms before and 160ms after a blink:
        for i = 1:length(first)
            if first(i) > 80
                vector((first(i)-80):first(i),1) = NaN;
            end
        end
        for i = 1:length(last)
            if last(i) < 160
                vector(last(i):(last(i)+160),1) = NaN;
            end
        end

        %% Count percentage of vector that is missing/blink:
        prop_missing = blink_samples / length(vector);
        if prop_missing > 0.2 %if more than 20% is missing
            warning('>20% of your pupil data is missing, consider removing trial')
        end
        if prop_missing > remove_prop
            output_vector = NaN(size(vector));
        else

            if method=="spline"
                vector_splinefillmiss = fillmissing(vector,'spline',1,...
                    'EndValues','nearest'); %will use spline interp but replace 
                        %the EndValues with the nearest values (and operates
                        %along dimension 1, which is what the 1 indicates)
            elseif method=="linear"
                vector_splinefillmiss = fillmissing(vector,'linear',1,...
                    'EndValues','nearest'); %will use linear interp but replace 
                        %the EndValues with the nearest values (and operates
                        %along dimension 1, which is what the 1 indicates)
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
        
        %% De-blink:
        
        % Calculate pupil mean and standard deviation (will need for 
        % identifying blinks):
        pupil_mean = nanmean(vector);
        pupil_std = std(vector,'omitnan');

        % Go through and identify blinks and remove them (make NaN)
        indices = find(vector(:,1)<(pupil_mean-(3*pupil_std)));
            %save the index numbers of where it is a blink
        blink_samples = sum(vector(:,1)<(pupil_mean-(3*pupil_std)));
            %count number of samples that are blinks
        vector(vector(:,1)<(pupil_mean-(3*pupil_std)),:) = NaN;
            %turn the blinks into NaN's
        D = (diff([0,diff(indices')==1,0]))'; 
            %find differences in indices to find consecutive blinked
            %samples --> because need to take out 80ms before to 160ms
            %after a blink (see Zekveld, Rudner, et al., 2014)
        first = indices(D>0); %first indices of consecutive blinks
        last = indices(D<0); %last indices of consecutive blinks
        % Deal with those blinks that occur within the first 80ms or the
        % last 160ms:
        ind_less_80 = indices(:,1)<=80; %logical array of any indices <=80
        if any(ind_less_80) %if any indices are within the first 80
            vector(1:80,1) = NaN; %then make all the first 80 NaN
        end
        ind_last_160 = indices(:,1)>=(size(vector,1)-160); 
            %logical array of any indices within last 160
        if any(ind_last_160) %if any indices are within last 160
            vector((size(vector,1)-160):size(vector,1)) = NaN;
                %then make all the last 160 NaN
        end
        % now take out the 80ms before and 160ms after a blink:
        for i = 1:length(first)
            if first(i) > 80
                vector((first(i)-80):first(i),1) = NaN;
            end
        end
        for i = 1:length(last)
            if last(i) < 160
                vector(last(i):(last(i)+160),1) = NaN;
            end
        end

        %% Count percentage of vector that is missing/blink:
        prop_missing = blink_samples / length(vector);
        if prop_missing > 0.2 %if more than 20% is missing
            warning('>20% of your pupil data is missing, consider removing trial')
        end

        %% Use interpolation to replace missing values:
        if method=="spline"
            vector_splinefillmiss = fillmissing(vector,'spline',1,...
                'EndValues','nearest'); %will use spline interp but replace 
                    %the EndValues with the nearest values (and operates
                    %along dimension 1, which is what the 1 indicates)
        elseif method=="linear"
            vector_splinefillmiss = fillmissing(vector,'linear',1,...
                'EndValues','nearest'); %will use linear interp but replace 
                    %the EndValues with the nearest values (and operates
                    %along dimension 1, which is what the 1 indicates)
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

    %% Make sure we don't have any negative values left after that 
    % (since pupil size can't be negative) and if we do then replace it
    % with the nearest value

    output_vector(output_vector<0) = NaN; 
        %replace any negative values with NaN
    output_vector = fillmissing(output_vector, 'nearest', 1);
        %replace any NaN's (which used to be negative) with the nearest
        %value (will default to the right if otherwise equal)


end % end function

