%Created May 10 2011, updated January 2013.  This program is for analysis of electrophysiological, behavioral, and stimulus data.  

probetype='probe3D_2x64E_1';                    %probe type.  see get_probegeometry for definitions.  
badchannels=[73 93 124];              %ok to leave empty. specifies the faulty channels on the probe.
backgroundchans1=['all'];            %default=['all']. can leave empty or write numeric list. The channels in the current set are not used in backgroundchans.  badchannels are removed from this list.
laser_artifact_removal='n';

trainingtrials=[20:-1:1];           %default=[20:-1:1]. can go backwards. used by make_seed_templates. setting the upper limit to more trials than exist is ok.

minamplitude=45;                    %default=50; can also leave empty and will calculate in make_seed_templates, but tends to be overestimate in high-spiking areas.   spike detection threshold, if use spikedetectionmethod=1.                                                                                                      

samplingrate=25000;                 %default=25000 Hz. Used to use 1/44.8e-6 (~22 kHz) pre-2012.  
f_low=600;                          %default=600.  bandpass filter, low frequency. 
f_high=6500;                        %default=6500. bandpass filter, high frequency.
n=3;                                %default=3.    number of poles butterworth filter

clusterstdev=12;                     %default=8. used in make_seed_templates.  Higher values produce fewer clusters.  Note that spikes larger/smaller than minamplitude are assigned bigger/smaller value of clusterstdev.
mergeclusterstdev=3;                %default=3. used in prune_templates & prune_penult_times
allcluststdev=[15];                 %default=[10]. used in run_template_matching. can use multiple values e.g. [6 8 10]. can be combined with subtract_templates='y'. note: values don't proportionally alter number of 'good' units.
max_tempupdatefraction=0.2;        %default=0.01. used in run_template_matching. maximum fraction to update template for next trial by the mean of matched spikes in the current trial.

discardSDfactor=4;                %default=4. used in get_final_units. minimum ratio of Vpp to s.d. to qualify as a unit. higher values disqualifies more units.
origmergeSDfactor=2;                %default=1.5. used in get_final_units. higher values merge more units.
basevoltagerange=200;               %default=200. used in a merging correction factor in get_final_units. used to prevent factor from growing too large if range(wavei)=range(wavej)
goodratebonus=1;                  %default=2. used in a merging correction factor get_final_units.
maxisidiff=0.005;
minburstfraction=0.05; 
maxburstisi=0.08;                   %default=0.08 sec. used in get_final_units.
do_pca_merge = 'n';                 %use old merge routine or pca_merge?

runAutoUnitQuality='s';             %default is 's' (semi-automatic). 'a'=automatic. 'n'=no. 'd'=demo mode that shows results for each unit.
final_minspikesperunit=200;         %default=200; used in get_final_units & auto_unit_quality. minimum number of spikes required for a unit to be scored >3.
autoVoltageCutoff=65;               %default=65 uV; used in auto_unit_quality.  minimum voltage required for a unit to be scored >3.

%note: remember there are additional plotting variables in make_savedir, set_plot_parameters, and get_unit_quality.

% *******************Finished setting parameters***********************

make_savedir;                       %asks user for raw data directory & makes savedir in corresponding \data analysis\ directory.

if startfromscratch=='y'

% share_data                        %asks user whether to share all files containing filename between specified network drives.
%rename_crashfiles                  %renames files if ephys program crashed during data acquisition.

split_channelsets                   %asks whether to split template making and matching tasks on multiple network drives. each drive will get a group of channel sets to work with.

% get_stimuli                         %extracts stimulus event (cue, solenoid & lick,..) times and saves in stimuli folder in data analysis. note: this program may be run separately from spikesort_muxi.

get_laserartifacts                  %optional: collects mean laser-induced artifact for every good channel, and subtracts this from raw data before doing any spike collection. note: artifact removal is not perfect.

make_seed_templates                 %creates templates from 'seed' sets of channels (usually nearest neighbors). stores as sortspikes_seti.mat
                                                                                                                        
make_orig_templates                 %final selection of unique templates. 
                                                                                          
run_template_matching               %matches data to templates on each set of channels. can use multiple values for allcluststdev and can do template subtraction from data.
                                                                     
get_penultimate_units               %does some light merging of units in their sets, and collects waveforms across all channels on a shaft.

combine_channelsets                 %if template making and matching tasks were split, will upload penultimate results to master drive, rename files and continue with get_final_units on the master drive.

get_final_units                     %discards units if they fail SD test, and merges similar units based on SD of waveforms across all channels. bad units are called 'badunits', and merged units are in 'mergeclusts'
  
plot_summary_muxi                   %plots waveforms, isi, psth for each unit and stores figures in \units\. 

get_unitquality                     %assigns score of 1, 2, or 3 to each unit. 

get_positions                       %assigns the position of the probe tips and determines unit coordinates relative to bregma.

%***Below this section: Post-Sorting Analysis***

% get_raw_LFP                         %extracts the unfiltered, downsampled signal for every channel. Used in LFP analysis.  

get_triggedspectra                  %gets event-triggered power spectral density for each trial and channel. 

get_triggedLFPpower                 %gets event-triggered LFP power in the specified LFP frequency band.

event_triggered                     %aligns unit activity to selected events & makes plots.

get_unitproperties 

plot_LFPstack                       %plots depth stack of event-triggered LFP for selected trials. contains additional plot settings inside subroutine script.

get_FTA                             %gets field-triggered average on all electrodes in specified frequency range.

plot_FTA                            %plots FTA data (raster of LFP peaks, FTA spike rate).

spike_LFP                           %obtains average LFP signal triggered on peaks, and plots the triggered LFP and PSTH vs depth. triggered on a reference channel. 

get_STA                             %gets spike-triggered average of LFP on all electrodes for each unit.

plot_STA                            %plots STA data (average STA signal, Fourier transform, power spectrum. 

get_multiunit                       %multiunit spike detection.

triggered_multiunit                 %event-triggered multiunit depth stacks. need stimuli file.
% cue_lick_multiephys               %cues (olfactory, visual, etc), lickometer & multiunit ephys analysis.

pairwise_correlation                %cross-correlation between all pairs of units.
 
% source_target_correlations        %average cross-correlation between groups of units in target and source locations.
% 
% spike_coherence_analysis          %calculates pairwise correlation (coincidence index) vs time, and coherence vs time & interneuronal distance.

% get_bar_parameters                 %visual stimulus analysis: obtains start and end times and angle of moving bar stimulus.
% bar_visual_response                %visual stimulus analysis: plots the psth for each unique moving bar angle.
% stimtrig_spike_LFP                 %visual stimulus analysis: visual grating stimuli. stimulus-triggered LFP & PSTH profiles. 
% movie_muxi                         %creates movie of raw data with spike audio. need to use virtualdub to compress and combine audio with video. does not require spike sorting.                                                                       
% spikes_position                   %use to evaluate any place-selective firing cells. need to have tracking.mat file.
% extracellular_analysis            %plots spike amplitude vs electrode distance from best channel.

end