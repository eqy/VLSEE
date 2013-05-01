close all;
%Pick two units
i=1;
j=2;

for i=1:length(spiketimes)
for j=2:length(spiketimes)
    if j <= i
        continue
    end
n_bins=40;

load([finalwavedir 'waveforms_i' num2str(1) '_cl' num2str(i) '.mat'])    
bestchan_i=bestchannel{i};
waves_i = waveforms{bestchan_i};
load([finalwavedir 'waveforms_i' num2str(1) '_cl' num2str(j) '.mat'])    
bestchan_j=bestchannel{j};
waves_j = waveforms{bestchan_j};

[n_waves_i, cols] = size(waves_i);
[n_waves_j, cols] = size(waves_j);

%Gives n: we are operating in R^n (currently unused)
n_dimension = length(waves_i(1,:));

%This idea of fake noise is difficult to justify, especially given the
%current unclarity of understanding of noise whitening. Don't have a solid
%understanding of how it's done in Rutishauser at the moment, just
%foolishly copying equations at face value. Unfortunately, this
%noise-whitening step seems to be essential in order to properly implement
%the projection-test. 
%Currently generating synthetic Gaussian noise with a standard deviaction
%of 8.9 uV.
noise = normrnd(0, 8.9, 1000000,47);
C = cov(noise);
R = chol(C);
waves_i = waves_i*(C\R);
waves_j = waves_j*(C\R);
temp=cov(waves_i);

%Assumes templates are just the means of the waveforms
u_i = mean(waves_i);
u_j = mean(waves_j);

%Unit vector defining the u_ij axis
n_u_ij = (u_j-u_i)/norm((u_j-u_i));

%Will store the projections of the [the difference between each event 
%and unit i] and [the u_ij axis
projections = zeros(n_waves_i + n_waves_j, 1);

%Projection for all the events in unit_i
for k=1:n_waves_i
    diff_i_k = waves_i(k,:) - u_i;
    projections(k) = dot(diff_i_k, n_u_ij);
end
%Projection for all the events in unit_j
for k=1:n_waves_j
    diff_i_k = waves_j(k,:) - u_i;
    projections(k+n_waves_i) = dot(diff_i_k, n_u_ij);
end
[nelements,xcenters] = hist(projections,n_bins);

gauss_dist_i = fitdist(projections(1:n_waves_i), 'Normal');
gauss_dist_j = fitdist(projections(1+n_waves_i:length(projections)), 'Normal');
gauss_i = pdf(gauss_dist_i, xcenters);
gauss_j = pdf(gauss_dist_j,xcenters);


bar(xcenters,nelements/trapz(xcenters,nelements));
hold on;
plot(xcenters,gauss_i, 'color', [1 0 0], 'linewidth', 2);
plot(xcenters,gauss_j, 'color', [0 1 0], 'linewidth', 2);
pause;
distance = abs(gauss_dist_i.Params(1) - gauss_dist_j.Params(1));
disp(['Projection test distance of unit ' num2str(i) ' and unit ' num2str(j) ' is ' num2str(distance)])
disp(['distance/std is '  num2str(distance/std(projections))])
hold off;
end
end
