
function [unitquality] = fully_auto_snr(timesdir, savedir, wavedir,...
    final_leftpoints,final_rightpoints)
    %set_plot_parameters
    final_spiketimes = load([timesdir 'finalspiketimes.mat']);
    load([savedir 'final_params.mat']);  %loads parameters file.
    n_units = length(final_spiketimes.spiketimes);
    unitquality = zeros(n_units, 1);
    bestchannel=parameters.bestchannels;
    percentiles = zeros(n_units,1);
    medians = zeros(n_units,1);
    n_ones = 0;
    n_twos = 0;
    n_threes = 0;
    
    
    for i = 1:n_units
        spikes = load([wavedir 'waveforms_i' num2str(1) '_cl' num2str(i) '.mat']);
        mean_waveform = mean(spikes.waveforms{bestchannel{i}});
        [n_waveforms, ~] = size(spikes.waveforms{bestchannel{i}});
        total_noise  = spikes.waveforms{bestchannel{i}} - repmat(mean_waveform,n_waveforms,1);
        beg_index = round(final_leftpoints/2);
        end_index = beg_index + round(final_rightpoints/2);
        noise = total_noise(:, beg_index:end_index);
        signal = mean_waveform(beg_index:end_index);
        
        noise_power = sum(noise.^2, 2);
        signal_power = sum(signal.^2);
        snr = signal_power./noise_power;
        percentiles(i) = prctile(snr,2);
        medians(i) = median(snr(snr > percentiles(i)));
        disp(['unit ' num2str(i) ' 2nd prctle = ' num2str(percentiles(i))]);
        disp(['unit ' num2str(i) ' median      = ' num2str(medians(i))]);
        %pause
    end
    
    ones = intersect(find(percentiles > 2), find(medians > 8));
    n_ones = length(ones);
    unitquality(ones) = 1;
    threes = setdiff(1:n_units, ones);
    n_threes = length(threes);
    unitquality(threes) = 3;
    disp(['Ones ' num2str(n_ones)]);
    disp(['Twos ' num2str(n_twos)]);
    disp(['Threes ' num2str(n_units - n_ones - n_twos)]);
    
end
