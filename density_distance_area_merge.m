%%%%%Density Distance Area, requires get_final_units to be complete
n_units = length(spiketimes);
max_isi=500;                                                                %micro-seconds
bin_size=5;                                                                %micro-seconds
for i=1:n_units;
    stimes_i = spiketimes{i};
    isi_i = diff(stimes_i) .* 1000;
    isi_i(isi_i>max_isi)=[];
    hist_i = hist(isi_i, max_isi/bin_size);
    
    %A = hist(isi_i, max_isi/bin_size);
    %A = A/sum(A);
    %plot(1:bin_size:max_isi, A, '.');
    
    for j=2:n_units;
        if j <= i
        continue;
        end
        stimes_j = spiketimes{j};
        isi_j = diff(stimes_j) .* 1000;
        isi_j(isi_j>max_isi)=[];
        hist_j = hist(isi_j, max_isi/bin_size);
        
        diff_i_j = abs(hist_i - hist_j);
        %if length(isi_i) < 500 || length(isi_j) < 500
        %    continue;
        %end
        density_distance_area_i_j = sum( (bin_size*diff_i_j) );
        if density_distance_area_i_j > 1000
            continue;
        end
        disp(['Density Distance Area of unit ' num2str(i) ' and unit ' num2str(j) ...
            ' : ' num2str(density_distance_area_i_j)])
    end
    
end
