function [pnn20, pnn50, stdnn, rmssd] = time_var_measures(rr_int)
% ----------------------------------------------------------------------- %
% calculate time domain measures of RR interval sequence variability
% input rr_int = input sequence of RR intervals
% output pnn20 = percentage of successive RR interval differences greater than 20 ms
% output pnn50 = percentage of successive RR interval differences greater than 50 ms
% output stdnn = standard deviation of RR intervals
% output rmssd = root mean square of differences of successive RR intervals
% ----------------------------------------------------------------------- %

rr_int_diffs = diff(rr_int);
pnn20 = sum(rr_int_diffs > 20) / length(rr_int_diffs);
pnn50 = sum(rr_int_diffs > 50) / length(rr_int_diffs);
stdnn = std(rr_int);

% first, calculate time differences between successive RR intervals. Square
% each difference value. Take the mean of all squared difference values,
% then take the square root.
rmssd = sqrt(mean(diff(rr_int).^2));

end