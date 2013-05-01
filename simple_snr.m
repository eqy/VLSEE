%%%Simple implementation of signal to noise pruning
function [garbage_units] = simple_snr(dounits, snr_bestchannel, snr_wavedir)
    TL = 3;
    dounits_len = length(dounits);
    snr = NaN(dounits_len,1);
    
   
    for i=1:dounits_len       
        unit_i_label = dounits(i);
        wave_i_data = load([snr_wavedir 'waveforms_i' num2str(1) '_cl' num2str(unit_i_label) '.mat']);
        waves_i = wave_i_data.waveforms{cell2mat(snr_bestchannel(unit_i_label))};
        template_i = mean(waves_i);
        [nrows, ncols] = size(waves_i);
        sig_ampl =  (1/ncols)*sum(template_i.^2);
        noise = waves_i - repmat(template_i, nrows, 1);  
        noise_ampl = (1/ncols)*sum(noise.^2,2);
        snr(i) = mean(sig_ampl./noise_ampl);
    end
    garbage_units_indicies = snr < TL;
    garbage_units = dounits(garbage_units_indicies);
end