function progressBar(curr_iter_num, tot_iter_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A progress indicator expressed in % of the total number of iterations.
%
% ARGUMENTS:
% curr_iter_num - integer indicating current iteration.
% tot_iter_num - integer indicating total number of iterations.
%
% OUTPUT: 
% None. It prints % to Matlab's workspace.
%
% REQUIRES: 
% None
%
% USAGE:
%{
	tot_iter_num = 1000;
    for this_it = 1:tot_iter_num,
        progressBar(this_it, tot_iter_num)
    end 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if curr_iter_num == 1,
        fprintf(1, '%3d%%', 0);
    else
        percentage_gap = 2;
        period = floor(percentage_gap*tot_iter_num/100);
        if mod(curr_iter_num, period) == 0
            fprintf('\b\b\b\b')
            fprintf(1, '%3d%%', percentage_gap*curr_iter_num/period);
        end
    end
end