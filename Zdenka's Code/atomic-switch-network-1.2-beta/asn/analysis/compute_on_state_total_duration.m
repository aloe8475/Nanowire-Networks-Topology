function compute_on_state_total_duration(filename)

if isempty(filename)

filename = 'script_single_switch_dcwait_duration';

end

dc_amplitude_on  = 2^-9:2^-9:2^-1;       % Volt (step 2^-6: 256 simulations)
dc_bias_residual = 2^-10:2^-10:2^-6; 

duration_map = zeros(length(dc_amplitude_on), length(dc_bias_residual));

for ii = 1:length(dc_amplitude_on)
	for kk=1:length(dc_bias_residual)
        
        try

            load(['2017-03-29' '-' filename '_on-' num2str(ii, '%04d') '_off-' num2str(kk, '%04d') '.mat'], 'Output', 'SimulationOptions', 'Stimulus')
        catch
            disp('not a file')
        end
	    c = 1./ Output.networkResistance;

	    idx = find(c > 2e-7);

	    if isempty(idx)
	    	duration_map(ii, kk) = 0;
	    else
	    	duration_map(ii, kk) = SimulationOptions.TimeVector(idx(end)) - SimulationOptions.TimeVector(idx(1));

	    end


end
end


save('2017-03-29_single_switch_total_duration_map_off_2e-7', 'duration_map')

end
