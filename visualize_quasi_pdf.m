%%%2013/02/25 Some tests to check for "split units"
%%Requires merging to be complete in get_final_units
spiketimes_len = length(spiketimes);
for i=1:spiketimes_len
    %Generate colors to be used in the plot
    rand('state',23)
    cmap = rand(10,3);

        %% 2013/02/27 New step to check ISI distribution for short intervals 
    isi_i = diff(spiketimes{i}) .* 1000;
    rounded_isi_i = round(isi_i);
    if(mode(rounded_isi_i) ~= 0)
        max_isi = mode(rounded_isi_i)*5;
    else
        max_isi = mean(isi_i)*5;
    end   
    max_isi = min(1000, max_isi);
    %max_isi=floor(median(isi_i)*2);
    isi_i(isi_i>max_isi)=[];
    if(length(isi_i) <= 10)
        continue;
    end
    
    w = input('Continue', 's');
    close all;
    %isi_i = round(isi_i);
        isi_dist_i = fitdist(isi_i', 'gamma');
        isi_dist_i2 = fitdist(isi_i', 'lognormal');
        isi_dist_i3 = fitdist(isi_i', 'rayleigh');
        isi_dist_i4 = fitdist(isi_i', 'weibull');
        isi_dist_i5 = fitdist(isi_i', 'nakagami');
        isi_dist_i6 = fitdist(isi_i', 'rician');
        isi_dist_i7 = fitdist(isi_i', 'gev');
        isi_dist_i8 = fitdist(isi_i', 'inversegaussian');
        
       
        set(0,'defaultlinelinewidth',3)
        %disp(num2str(isi_dist_i.Params));
        scrsz = get(0,'ScreenSize');
        figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
        ylabel('value');
        xlabel('ISI (ms)');
        title(['Unit ' num2str(i)]);
        hold on;
        

        
        h = hist(isi_i, 1:max_isi);
        h = h/sum(h);
        hist_plot = plot(h,'.','color',cmap(1,:)); w = input('Continue', 's');
        
        gamma_plot = plot(1:max_isi, gampdf( 1:max_isi, isi_dist_i.Params(1)...
            , isi_dist_i.Params(2)),'color',cmap(2,:)); w = input('Continue', 's');
        logn_plot = plot(1:max_isi, lognpdf(1:max_isi, isi_dist_i2.Params(1)...
            , isi_dist_i2.Params(2)),'color',cmap(3,:)); w = input('Continue', 's');
        rayleigh_plot = plot(1:max_isi, raylpdf(1:max_isi, isi_dist_i3.Params(1)...
            ),'color',[0 1 1]); w = input('Continue', 's');
        weibull_plot = plot(1:max_isi, wblpdf(1:max_isi, isi_dist_i4.Params(1)...
           , isi_dist_i4.Params(2)),'color',cmap(4,:)); w = input('Continue', 's');
        nakagami_plot = plot(1:max_isi, pdf('nakagami', 1:max_isi, isi_dist_i5.Params(1)...
            , isi_dist_i5.Params(2)),'color', cmap(5,:)); w = input('Continue', 's');
        rician_plot = plot(1:max_isi, pdf('rician', 1:max_isi, isi_dist_i6.Params(1)...
            , isi_dist_i6.Params(2)),'color', cmap(6,:)); w = input('Continue', 's');
        gev_plot = plot(1:max_isi, pdf('gev', 1:max_isi, isi_dist_i7.Params(1)...
            , isi_dist_i7.Params(2)),'color',cmap(7,:)); w = input('Continue', 's');
        invgaussian_plot = plot(1:max_isi, pdf('inversegaussian', 1:max_isi, isi_dist_i8.Params(1)...
            , isi_dist_i8.Params(2)),'color',cmap(8,:)); w = input('Continue', 's');
        
        
        
        
        
        legend([gamma_plot, logn_plot, rayleigh_plot, weibull_plot, nakagami_plot, ...
            rician_plot, gev_plot, invgaussian_plot, hist_plot]  ...
            ,'Gamma quasi-pdf', 'Log-Normal quasi-pdf', 'Rayleigh quasi-pdf'...
            , 'Weibull quasi-pdf', 'Nakagami quasi-pdf', 'Rician quasi-pdf' ...
            , 'Generalized Extreme Value quasi-pdf', 'Inverse Gaussian quasi-pdf'...
            , 'Normalized Histogram');
        hold off;
        disp(num2str(mean(isi_i)));
end   
        %plot(hist(isi_i,1:100));
        %hold on;
        %plot(0:100, poisspdf(0:100,isi_dist_i.Params));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit_i_max_isi = max(diff(spiketimes{i}));
%disp(['Unit ' num2str(i) ' has max ISI of ' num2str(unit_i_max_isi)])
%     if unit_i_max_isi > 100
%         for j=1:spiketimes_len
%             if j <= i  
%             continue
%             end
%             unit_j_max_isi = max(diff(spiketimes{j}));
%             if unit_j_max_isi > 100
%                 tempspiketimes = [spiketimes{i} spiketimes{j}];
%                 tempspiketimes = unique(tempspiketimes);
%                 tempspiketimes = sort(tempspiketimes);
%                 new_max_isi = max(diff(tempspiketimes));
%                 if new_max_isi < 50
%                     disp(['Unit ' num2str(i) ' and unit ' num2str(j) ' are favored to merge'])
%                     disp(['Original Max ISI of unit ' num2str(i) ' : ' num2str(unit_i_max_isi)])
%                     disp(['Original Max ISI of unit ' num2str(j) ' : ' num2str(unit_j_max_isi)])
%                     disp(['Merged Max ISI : ' num2str(new_max_isi)])
%         
%                 end
%             end
%         end 
%     end
%end
