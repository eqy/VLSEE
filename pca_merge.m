function [merges] = pca_merge(dounits, pca_bestchannel, pca_wavedir)
    PCA_TH = 5;
    PCT_DIFF_TH = 0.30;                   %Vmax and Vmin percentage difference threshold 
    dounits_len = length(dounits);
    distances = NaN(length(dounits));
    vmax_pctdiff = NaN(length(dounits));
    vmin_pctdiff = NaN(length(dounits));
    
    parfor i=1:dounits_len
        unit_i_label = dounits(i);
        wave_i_data = load([pca_wavedir 'waveforms_i' num2str(1) '_cl' num2str(unit_i_label) '.mat']);
        bestchan_i = pca_bestchannel{unit_i_label};
        waves_i = wave_i_data.waveforms{bestchan_i};
        [n_i, m_i] = size(waves_i);
        distances_j = NaN(1,dounits_len);
        vmax_pctdiff_j = NaN(1,dounits_len);
        vmin_pctdiff_j = NaN(1,dounits_len);
        for j=1:dounits_len
            unit_j_label = dounits(j);
            if unit_j_label <= unit_i_label
               continue;
            end       

            wave_j_data = load([pca_wavedir 'waveforms_i' num2str(1) '_cl' num2str(unit_j_label) '.mat']);
            bestchan_j = pca_bestchannel{unit_j_label};
            if bestchan_j ~= bestchan_i
                continue;
            end
            waves_j = wave_j_data.waveforms{bestchan_j};
            %%%Percentage Difference Calculation
                template_i = mean(waves_i);
                template_j = mean(waves_j);
                max_i = max(template_i);
                max_j = max(template_j);
                min_i = min(template_i);
                min_j = min(template_j);
                vmax_pctdiff_j(j) =  abs(max_i - max_j)/(abs(max_i+max_j)/2);
                vmin_pctdiff_j(j) =  abs(min_i - min_j)/(abs(min_i+min_j)/2); 
            %%%
            
            waves_i_j = [waves_i ; waves_j];
            standard_waves_i_j = zscore(waves_i_j);
            [c,s] = pca(standard_waves_i_j);

            [n_j, ~] = size(waves_j);
            
            m_i_x = mean(s(1:n_i,1));
            m_i_y = mean(s(1:n_i,2));
            m_i_z = mean(s(1:n_i,3));

            m_j_x = mean(s(1+n_i:n_i+n_j,1));
            m_j_y = mean(s(1+n_i:n_i+n_j,2));
            m_j_z = mean(s(1+n_i:n_i+n_j,3));
                        
            dist_i_j = sqrt((m_i_x-m_j_x)^2 + (m_i_y-m_j_y)^2 + (m_i_z-m_j_z)^2);

            distances_j(j) = dist_i_j;

        end
        distances(i,:) = distances_j;
        vmax_pctdiff(i,:) = vmax_pctdiff_j;
        vmin_pctdiff(i,:) = vmin_pctdiff_j;
    end
    [pca_pass_i, pca_pass_j] = find(distances < PCA_TH);
    [vmax_i_pass, vmax_j_pass] = find(vmax_pctdiff <= PCT_DIFF_TH);
    [vmin_i_pass, vmin_j_pass] = find(vmin_pctdiff <= PCT_DIFF_TH);
    pca_pass = [pca_pass_i pca_pass_j];
    vmax_pass = [vmax_i_pass vmax_j_pass];
    vmin_pass = [vmin_i_pass vmin_j_pass];
    
    merge_indicies = intersect(vmax_pass, vmin_pass, 'rows');
    merge_indicies = intersect(pca_pass, merge_indicies, 'rows');
    merge_i = dounits(merge_indicies(:,1));
    merge_j = dounits(merge_indicies(:,2));
    merges = [merge_i' merge_j'];
end