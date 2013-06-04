function [garbage_units] = get_sane(dounits, spiketimes, bestchannel, ...
    wavedir, sampling_rate, ispen)
    %pp thresholds are for the maximum allowable peak-to-peak voltages and
    %the minimum allowable peak to peak voltages
    VPP_TH = 350;   %Vpp Threshold High
    VPP_TL = 60;    %Vpp Treshold Low
    %Absolute thresholds are for the maximum allowable voltage and the
    %minimum allowable voltages
    AV_TH = 350;    %Absolute Voltage High
    AV_TL = -350;   %Absolute Voltage Low
    %Coefficient of variation thresholds are used to prune units that are
    %noisy around the minimum, lower values discard more units
    COFV_MIN_TH = 0.22;  %Coefficient of Variation Threshold High
    COFV_MAX_TH = 0.25;  %Currently UNUSED 
    %Quasi-derivative test, used to eliminate excessively "flat" spikes
    DIFF_MAX_TL = 25;
    
    %Set the minimum number of spikes allowed based on whether
    %get_final_units has finished
    if ispen==1
        WC_TL = 0;
    else
        WC_TL = 200;    %Wave Count Threshold Low
    end
    
    %Preallocate arrays, NaN prevents confusion e.g. if a min() is used and
    %there are zeros on the matrix
    garbage_indicies = [];
    cofv_min = NaN(1,length(dounits));
    cofv_max = NaN(1,length(dounits));
    vpp = NaN(1,length(dounits));
    vmax = NaN(1,length(dounits));
    vmin = NaN(1,length(dounits));
    diff_max = NaN(1,length(dounits));
    diff_mean = NaN(1,length(dounits));
        
    parfor i=1:length(dounits)    
        %Load each waveform and compute the values of each of the
        %parameters
        unit = dounits(i);
        if length(spiketimes{unit}) < WC_TL
            garbage_indicies = [garbage_indicies i];
            continue;
        end
        current = load([wavedir 'waveforms_i' num2str(1) ...
            '_cl' num2str(unit) '.mat']);
        channel = bestchannel(unit);
        best_channel_waves = current.waveforms{cell2mat(channel)};
        [n_waves, ~] = size(best_channel_waves);
        if(n_waves > 1)
            mean_waveform = mean(best_channel_waves);
        else
            mean_waveform = best_channel_waves;
        end
        vmax(i) = max(mean_waveform);
        vmin(i) = min(mean_waveform);
        vpp(i) = max(mean_waveform) - min(mean_waveform);
        
        min_index = find(mean_waveform == vmin(i));
        all_mins = best_channel_waves(:,min_index);
        cofv_min(i) = abs(std(all_mins)/mean(all_mins));
        
        max_index = find(mean_waveform == vmax(i));
        all_max = best_channel_waves(:,max_index);
        cofv_max(i) = abs(std(all_max)/mean(all_max));
        diff_max(i) = max(abs(diff(mean_waveform)));    
    end
    %Concatenate all the indicies that fail the thresholds, then removed
    %redunant indicies
    garbage_indicies = [garbage_indicies find(vpp > VPP_TH)];
    garbage_indicies = [garbage_indicies find(vpp < VPP_TL)];
    garbage_indicies = [garbage_indicies find(vmax > AV_TH)];
    garbage_indicies = [garbage_indicies find(vmin < AV_TL)];
    garbage_indicies = [garbage_indicies find(diff_max < DIFF_MAX_TL)];
    cofv_min_fail = find(cofv_min > COFV_MIN_TH);
    cofv_max_fail = find(cofv_max > COFV_MAX_TH);
    %garbage_indicies = [garbage_indicies cofv_max_fail];
    garbage_indicies = [garbage_indicies cofv_min_fail];
    garbage_units = dounits(garbage_indicies); 
    garbage_units = sort(unique(garbage_units));
end

