function [avg, err] = jh_mean_err(data,err_type)
% author: Junghan Shin
% contact: cautious@kaist.ac.kr
% last modified: 2019.12.07.
% description
%       return mean and error from samples
% inputs
%       data: samples in 'cell' type
%       err_type: 1 - standard error(SE) / 2 - 95% CI(1.96 * SE)
% outputs
%       mean and errors
% example
%       jh_mean_err(1,{[1 2 3 4 5],[1 3 5],[6 8 10],[5 8 9 8]})

% argument parsing
switch nargin
    case 1
        err_type = 1;
    case 2
    otherwise
        error('invalid input');
end

% error handling
if ~iscell(data)
    error('data should be cell type');
end

if err_type ~=1 && err_type ~=2
    error('invalid error type');
end

% calculate average and error
avg = cellfun(@nanmean,data); 
err = cellfun(@(x) nanstd(x)/sqrt(length(x)), data);
if err_type == 2
    err = err.* 1.96;
end    