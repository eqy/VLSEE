disp(['get_final_units'])
disp('Custom version that uses get_sane and pca_merge')
%Dec 29 2011.  A limitation of the current implementation is that when
%trying to merge units, if there are sufficient similar but slightly
%different unique units they will all be combined if origmergeSDfactor is too
%high. To avoid this "unit sink" problem, can try to implement a more
%gradual annealing process in which waveform means are updated after every
%merged pair of units, rather than lump everything together at the end. 
%Can also try to use ISI information.

set_plot_parameters
maxtrial=parameters.maxtrial;
trialduration=parameters.trialduration;
load([timesdir 'penultimate_params.mat']);  %loads parameters file.
numberoftrials=maxtrial-1;
maxtime=(numberoftrials*trialduration)/samplingrate;   %added 1/3/13 for new merge correction factor

isiplottime=0:isibinsize:isirange;

minamplitude=parameters.minamplitude;
lengthperchan=parameters.lengthperchan{1};   %this variable should be same for all units.

load([timesdir 'penult_spiketimes.mat'])   %loads spiketimes created in collect_spiketimes;
load([timesdir 'penult_baretimes.mat'])   %loads spiketimes created in collect_spiketimes;
load([timesdir 'penult_jittertimes.mat'])   %loads spiketimes created in collect_spiketimes; 
load([penultwavedir 'bestchannel.mat']);  %loads parameters file.

parameters.bestchannel=bestchannel;


close all
scrsz=get(0,'ScreenSize');
timestarting=datenum(clock)*60*24;   %starting time in minutes.

dounits=1:length(spiketimes);
origdounits=dounits;

maxjitter=min([(leftpoints-origleftpoints) (rightpoints-origrightpoints)]);
if maxjitter>=upsamplingfactor;
jittersamples=-upsamplingfactor:2:upsamplingfactor;
jittersamples=setdiff(jittersamples,0);
jittersamples=[0 jittersamples];
a=find(abs(jittersamples)>rejecttime);
jittersamples(a)=[];
else jittersamples=0;
end

disp(['final merging step: checking ' num2str(length(dounits)) ' candidate units.'])
disp(['checking amplitude to SD of waveforms.'])
meanwaves=[]; sdwaves=[]; badunits=[]; 
for unitind=1:length(dounits);
    unit=dounits(unitind);   
%     lengthperchan=parameters.lengthperchan{unit};
    bestchan=bestchannel{unit};
    
    if length(bestchan)==0    %fixes very rare bug
        badunits=[badunits unit];
        continue
    end
    
    currentshaft=s.shaft(bestchan);    %added Mar 10 2012.
    parameters.shaft{unit}=currentshaft;
    timesuniti=spiketimes{unit};
    if length(timesuniti)<2
        continue
    end
    
    stimesi=spiketimes{unit};
    difftimes=diff(stimesi);
    fractionbursts=length(find(abs(difftimes)<maxburstisi))/length(difftimes);
      
    load([penultwavedir 'waveforms_i' num2str(1) '_cl' num2str(unit) '.mat'])
    if length(waveforms{bestchan})==0
        badunits=[badunits unit];
        continue
    end
    
  meanwavesik=[]; sdwavesik=[];  
  for k=1:length(jittersamples);
  jitterk=jittersamples(k);
    
    meanwaveuniti=[]; sdwaveuniti=[]; previousmaxmean=0;
    for j=1:length(waveforms);   %for each channel.
        
        if length(waveforms{j})==0    %added June 22 2012.
            continue
        end
        
        if s.shaft(j)~=currentshaft    %added Mar 10 2012.
            continue
        end

        t0=leftpoints-origleftpoints-jitterk;
        tf=t0+origleftpoints+origrightpoints;
        
        waveschanj=waveforms{j}(:,t0:tf);
        
        if size(waveforms{j},1)>1
        meanwaveuniti=[meanwaveuniti, mean(waveschanj)];
        else
        meanwaveuniti=[meanwaveuniti, waveschanj];
        end
        sdwaveuniti=[sdwaveuniti, std(waveschanj)];
        
        isbadchannel=length(find(badchannels==j));
        newdiscardSDfactor=discardSDfactor*(1-fractionbursts/2);
           
        %note: I found this metric to be the most reliable in identifying bad units. 28/12/11
        ratiomeanstd=range(mean(waveschanj))/max(std(waveschanj));  %calculates the ratio of mean range to the largest s.d.
       
        if j==bestchan & isbadchannel==0 & k==1 & ratiomeanstd<newdiscardSDfactor     
            badunits=[badunits unit];   
            continue
        end    
        
        if j==bestchan & isbadchannel==0 & k==1 & abs(mean(waveschanj(:,extraleft)))<minamplitude  %discard unit if peak position is < minamplitude
            badunits=[badunits unit];       
        end

    end
    meanwavesik{k}=meanwaveuniti;
    sdwavesik{k}=sdwaveuniti;
    
  end
    meanwaves{unit}=meanwavesik;
    sdwaves{unit}=sdwavesik;

end



badunits=unique(badunits);
if do_get_sane == 'y'
    display('using get_sane')
    tic
    badunits2=get_sane(dounits,spiketimes,bestchannel,penultwavedir,25000,1);
    toc 
    display('during get_sane')
    old_len = length(badunits);
    badunits=[badunits badunits2];
    badunits=unique(badunits);
    new_len = length(badunits);
    disp([num2str(new_len-old_len) ' units were uniquely identified as bad by get_sane'])
else
    display('not using get_sane')
end

dounits=setdiff(dounits,badunits);
disp(['discarded ' num2str(length(badunits)) ' units which failed quality checks, ' num2str(length(dounits)) ' units remaining.'])

sdrate=[];  %added 1/3/13
for i=1:length(dounits);
    uniti=dounits(i);
    ratei=histc(spiketimes{uniti},0:100:maxtime)/100;
    sdrate{uniti}=std(ratei);
end

disp(['checking for units to merge.'])
samemeans=[];

if do_pca_merge == 'n'
    disp('not using pca_merge')
    for i=1:length(dounits);
        uniti=dounits(i);

        stimesi=spiketimes{uniti};

        meanwavei=meanwaves{uniti}{1};

        isii=diff(stimesi);
        isii=isii(find(isii<=isirange));
        if length(isii)>0
        histisii=100*histc(isii,isiplottime)/length(isii);
        histisii=smooth(histisii,20);
        peakisitimei=isiplottime(find(histisii==max(histisii)));
        peakisitimei=peakisitimei(1);
        else peakisitimej=100;
        end

        for j=2:length(dounits);
          unitj=dounits(j);
          if i>=j | parameters.shaft{uniti}~=parameters.shaft{unitj}
             continue
          end

          stimesj=spiketimes{unitj}; 

          isij=diff(stimesj);
          isij=isij(find(isij<isirange));
          if length(isij)>0
          histisij=100*histc(isij,isiplottime)/length(isij);
          histisij=smooth(histisij,20);
          peakisitimej=isiplottime(find(histisij==max(histisij)));
          peakisitimej=peakisitimej(1);
          else peakisitimej=100;
          end

          combinedrate=histc(sort([stimesi stimesj]),0:100:maxtime)/100; %added 1/3/13 for new merge correction factor
          sdcombinedrate=std(combinedrate);

          timeratio=length(stimesi)/length(stimesj);
          if timeratio>=1  %equalize # of times.
          newstimesi=stimesi(1:length(stimesj));
          newstimesj=stimesj;
          else
          newstimesi=stimesi;
          newstimesj=stimesj(1:length(stimesi));
          end
          difftimes=abs(newstimesj-newstimesi);                                        

          for k=1:length(jittersamples);
          jitterk=jittersamples(k);

            meanwavej=meanwaves{unitj}{k};

            diffwaves=abs(meanwavei-meanwavej);     %THIS IS THE SLOWEST LINE IN LOOP. use only ~1 to 3 points instead of entire waveform; implemented on 7/12/11 to deal with Dec6b data.

            minsd=min([sdwaves{uniti}{1}; sdwaves{unitj}{k}]);   %default minsd=min([sdwaves{uniti}{1}; sdwaves{unitj}{k}]);

            mindiff=minsd;       

    %       mergeSDfactor=origmergeSDfactor; %no fudge factor (27/2/13)

            if bestchannel{uniti}==bestchannel{unitj}     %added 28/2/13
            fudgefactor=(range(meanwavei)+range(meanwavej))/(abs(range(meanwavei)-range(meanwavej))+basevoltagerange); 
                if sdcombinedrate<sqrt((sdrate{uniti})^2+(sdrate{uniti})^2)  %multiplies a bonus factor if the S.D. of the combined firing rate is < than that of the individual rates.
                    fudgefactor=goodratebonus*fudgefactor;
                end

    %             if abs(peakisitimei-peakisitimej)<maxisidiff & peakisitimei<2*isirange & peakisitimej<2*isirange;
    %                 fudgefactor=(range(meanwavei)+range(meanwavej))/(abs(range(meanwavei)-range(meanwavej))+basevoltagerange); 
    %             end

            else fudgefactor=-origmergeSDfactor+0.2;
            end
            mergeSDfactor=origmergeSDfactor+fudgefactor;

            nomatches=length(find(diffwaves>mergeSDfactor*mindiff));

            if sum(nomatches)>0    %if it appears there is no match between uniti and unitj, then try a second test to see if they match
            test2=length(find(diffwaves>(2*mergeSDfactor)*mindiff));  %second test relaxes merge criteria by a factor of two, but requires that the units essentially be identical.
            diffcombinedtimes=diff(sort([stimesi stimesj]));
            fractionreject=length(find(diffcombinedtimes<=rejecttime/samplingrate))/length(diffcombinedtimes);
                if sum(test2)<=0 & fractionreject>0.32   %default=0.68
                nomatches=0; 
                end
            end

            if sum(nomatches)<=0     
                samemeans=[samemeans; uniti unitj];   %the left and right columns specify indices to be merged. 
                continue   %continue is for the jittersamples loop         
           end

          end

        end

    end
else
    disp(['using pca_merge'])
    tic
    [samemeans, debug] = pca_merge(dounits,bestchannel,penultwavedir);
    toc
    disp(['during pca_merge for ' num2str(length(dounits)) ' units'])
    disp(['pca_merge decided on ' num2str(length(samemeans)) ' merges'])
end
  

newtimes=[]; newbaretimes=[]; newjittertimes=[]; newshaft=[];
unitcounter=1;  mergeunits=[]; usedunits=[];

if length(samemeans)>0
 
a=unique(samemeans);

    for i=1:length(a);
        uniti=a(i);
        if length(intersect(usedunits,uniti))>0
            continue
        end
        
        tempmerge=uniti;
        for j=1:100;
            b=ismember(samemeans, tempmerge);
            [c,d]=find(b==1);         
            tempmerge=unique([samemeans(c,1) samemeans(c,2)]);
        end
                      
       mergeunits{i}=tempmerge;     
    
       usedunits=[usedunits mergeunits{i}'];
    end
     
    newmergeunits=[]; mergecounter=1;
    for j=1:length(mergeunits)
        if length(mergeunits{j})>0
            newmergeunits{mergecounter}=mergeunits{j};
            mergecounter=mergecounter+1;
        end
    end
    mergeunits=newmergeunits;
   
    for ind=1:length(mergeunits);  
           
        allinds=unique(mergeunits{ind});  %all indices to merge.
        firstind=allinds(1); 
        mergetimes=[]; mergebaretimes=[]; mergejittertimes=[];        
        for i=1:length(allinds);
            indi=allinds(i);       
            
            mergetimes=[mergetimes spiketimes{indi}];   
            mergebaretimes=[mergebaretimes baretimes{indi}];
            mergejittertimes=[mergejittertimes jittertimes{indi}];
        end
            
        [mergetimes, sortindex]=sort(mergetimes);   
        mergebaretimes=mergebaretimes(sortindex);
        mergejittertimes=mergejittertimes(sortindex);
        
        difftimes=diff(mergetimes);
        overlaps=find(difftimes<(rejecttime+minoverlap)/samplingrate);          %default=rejecttime (modified Jan 13 2012 to reduce ISI violations). indices of spike which has isi<rejecttime with another spike.
        mergetimes(overlaps)=[];  %removes duplicate events found on multiple channels.
        mergebaretimes(overlaps)=[];
        mergejittertimes(overlaps)=[];

        newtimes{unitcounter}=mergetimes;
        newbaretimes{unitcounter}=mergebaretimes;
        newjittertimes{unitcounter}=mergejittertimes; 
        newshaft{unitcounter}=parameters.shaft{firstind};
        
        unitcounter=unitcounter+1; 
 
    end
    
end


allinds=setdiff(dounits,usedunits);  %now do the unmerged units, rejecting the bad units.
 
    for ind=1:length(allinds);    
        unit=allinds(ind);
        
        stimesi=spiketimes{unit};
        baretimesi=baretimes{unit};
        jittertimesi=jittertimes{unit};
        
        [stimesi, sortindex]=sort(stimesi);   
        baretimesi=baretimesi(sortindex);
        jittertimesi=jittertimesi(sortindex);
        
        difftimes=diff(stimesi);
        overlaps=find(difftimes<(rejecttime+minoverlap)/samplingrate);          %default=rejecttime (modified Jan 13 2012 to reduce ISI violations). indices of spike which has isi<rejecttime with another spike.  
        stimesi(overlaps)=[];  %removes duplicate events found on multiple channels.
        baretimesi(overlaps)=[];
        jittertimesi(overlaps)=[]; 
        
        if length(stimesi)<final_minspikesperunit
            badunits=sort([badunits unit]);
            continue
        end
        
        newtimes{unitcounter}=stimesi;
        newbaretimes{unitcounter}=baretimesi;
        newjittertimes{unitcounter}=jittertimesi;
        newshaft{unitcounter}=parameters.shaft{unit};
        
        unitcounter=unitcounter+1; 
 
    end

spiketimes=newtimes;
baretimes=newbaretimes;
jittertimes=newjittertimes;

disp(['discarded another ' num2str(length(intersect(allinds,badunits))) ' units because they contain < ' num2str(final_minspikesperunit) ' spikes.'])

disp(['pruned final number of units from ' num2str(length(origdounits)) ' to ' num2str(length(spiketimes)) '.'])

save([timesdir 'finalspiketimes.mat'],'spiketimes', '-mat')
save([timesdir 'finalbaretimes.mat'],'baretimes','-mat')
save([timesdir 'finaljittertimes.mat'],'jittertimes','-mat')

save_finalparameters

final_times_to_waves

disp(['finding final spike amps & best channel per unit.'])
bestchannel=[]; spikeamps=[];
for unitind=1:length(dounits);
    unit=dounits(unitind); 
    
    load([finalwavedir 'waveforms_i' num2str(1) '_cl' num2str(unit) '.mat'])
    
    maxamp=0; bestchaninunit=[];
    shaftinuse=parameters.shaft{unit};   %finds the shaft containing these channels; should only be single shaft.
    channelsonshaft=find(s.shaft==shaftinuse);    %finds all the channels on that shaft.    
    for j=1:length(channelsonshaft);
        chj=channelsonshaft(j);
        spikeamps{unit}{chj}=range(waveforms{chj}');
        if range(mean(waveforms{chj}))>maxamp
            maxamp=range(mean(waveforms{chj}));
            bestchaninunit=chj;
        end
    end
    bestchannel{unit}=bestchaninunit;
end
save([finalwavedir 'spikeamps.mat'],'spikeamps','-mat')
save([finalwavedir 'bestchannel.mat'],'bestchannel','-mat')
clear spikeamps
timeending=datenum(clock)*60*24; %final time in minutes
telapsed=timeending-timestarting;
disp([num2str(telapsed) ' minutes elapsed during get_final_units.'])
