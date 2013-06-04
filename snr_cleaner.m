%Currently unused, potentially use to "clean up" units later on

function [waveforms] = snr_cleaner(timesdir, wavedir, savedir, i, final_leftpoints, final_rightpoints))
 final_spiketimes = load([timesdir 'finalspiketimes.mat']);
 load([savedir 'final_params.mat']);  %loads parameters file.

end