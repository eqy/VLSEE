function [garbage_units] = get_sane(dounits, spiketimes, bestchannel, wavedir, sampling_rate, ispen)
    VPP_TH = 350;   %Vpp Threshold High
    VPP_TL = 60;    %Vpp Treshold Low
    AV_TH = 350;    %Absolute Voltage High
    AV_TL = -350;   %Absolute Voltage Low
    
    COFV_TH = 0.25;  %Coefficient of Variation Threshold Low
    if ispen==1
        WC_TL = 0;
    else
        WC_TL = 200;    %Wave Count Threshold Low
    end
    garbage_indicies = [];
    cofv = NaN(1,length(dounits));
    vpp = NaN(1,length(dounits));
    vmax = NaN(1,length(dounits));
    vmin = NaN(1,length(dounits));
    parfor i=1:length(dounits)
        unit = dounits(i);
        if length(spiketimes{unit}) < WC_TL
            garbage_indicies = [garbage_indicies i];
            continue;
        end
        current = load([wavedir 'waveforms_i' num2str(1) '_cl' num2str(unit) '.mat']);
        channel = bestchannel(unit);
        best_channel_waves = current.waveforms{cell2mat(channel)};
        mean_waveform =  mean(best_channel_waves);
        vmax(i) = max(mean_waveform);
        vmin(i) = min(mean_waveform);
        vpp(i) = max(mean_waveform) - min(mean_waveform);
        
        min_index = find(mean_waveform == min(mean_waveform));
        all_mins = best_channel_waves(:,min_index);
        cofv(i) = abs(std(all_mins)/mean(all_mins));
    end
   
    garbage_indicies = [garbage_indicies find(vpp > VPP_TH)];
    garbage_indicies = [garbage_indicies find(vpp < VPP_TL)];
    garbage_indicies = [garbage_indicies find(vmax > AV_TH)];
    garbage_indicies = [garbage_indicies find(vmin < AV_TL)];
    garbage_indicies = [garbage_indicies find(cofv > COFV_TH)];
    
    garbage_units = dounits(garbage_indicies); 
    garbage_units = [garbage_units simple_snr(dounits, bestchannel, wavedir, 4)]; %Add simple snr pruning to get_sane
   
    garbage_units = sort(unique(garbage_units));
end