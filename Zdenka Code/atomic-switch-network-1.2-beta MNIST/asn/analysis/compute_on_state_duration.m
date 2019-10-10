function compute_on_state_duration(filename)

if isempty(filename)

filename = 'script_cluster_run_dcwait_duration_simulation';

end

dc_amplitude_on  = 2^-10:2^-6:2^1;     % Volt (step 2^-6: 128 simulations)
dc_bias_residual = 2^-10:2^-10:2^-4; 

duration_map = zeros(length(dc_amplitude_on), length(dc_bias_residual));

for ii = 1:length(dc_amplitude_on)
	for kk=1:length(dc_bias_residual)

		load([datestr(now, 29) '-' filename '_on-' num2str(ii, '%04d') '_off-' num2str(kk, '%04d') '.mat'], 'Output', 'SimulationOptions', 'Stimulus')

        
	    c = 1./ Output.networkResistance;
	    c_off = c(floor(SimulationOptions.NumberOfIterations/2):end);

	    idx = find(c_off < 1e-4);
	    idx = idx + floor(SimulationOptions.NumberOfIterations/2) - 1;

	    if isempty(idx)
	    	duration_map(ii, kk) = SimulationOptions.T - Stimulus.OffTime;
	    else
	    	duration_map(ii, kk) = SimulationOptions.T - SimulationOptions.TimeVector(idx(1));

	    end


end
end


save('2017-03-26_nw-00100-max_av_duration_map_off_1e-4', 'duration_map')

end
