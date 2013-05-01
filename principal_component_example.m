%i = 23;
%j = 24;

example=NaN(length(spiketimes));
example2=NaN(length(spiketimes));
for i=1:length(spiketimes)
    
        wave_i_data = load([finalwavedir 'waveforms_i' num2str(1) '_cl' num2str(i) '.mat']);
        bestchan_i=bestchannel{i};
        waves_i = wave_i_data.waveforms{bestchan_i};


        %waves_i_s_1 = median(s(:,1));
        %waves_i_s_2 = median(s(:,2));
        %waves_i_s_3 = median(s(:,3));
        

        [n_i, m_i] = size(waves_i);
        example2_j = NaN(1,length(spiketimes));
        
        tic
    for j=1:length(spiketimes)
        if j <= i
           continue;
        end       
       
        wave_j_data = load([finalwavedir 'waveforms_i' num2str(1) '_cl' num2str(j) '.mat']) ;   
        bestchan_j=bestchannel{j};
        
        if bestchan_j ~= bestchan_i
            continue;
        end
        waves_j = wave_j_data.waveforms{bestchan_i};

        waves_i_j = [waves_i ; waves_j];
        zs_i_j = zscore(waves_i_j);
            
        [c,s] = pca(waves_i_j);
        [zs_i_j,zs] = pca(zs_i_j);
        [n_j, ~] = size(waves_j);
        %scatter(zs(1:n_i,1), zs(1:n_i,2),3,[1 0 0]);
        hold on;
        %scatter(zs(n_i+1:n_i+n_j,1), zs(n_i+1:n_i+n_j,2),3,[0 0 1]);
        
        m_i_x = mean(s(1:n_i,1));
        m_i_y = mean(s(1:n_i,2));
        m_i_z = mean(s(1:n_i,3));
        
        m_j_x = mean(s(1+n_i:n_i+n_j,1));
        m_j_y = mean(s(1+n_i:n_i+n_j,2));
        m_j_z = mean(s(1+n_i:n_i+n_j,3));
        
        zm_i_x = mean(s(1:n_i,1));
        zm_i_y = mean(s(1:n_i,2));
        zm_i_z = mean(s(1:n_i,3));
        
        zm_j_x = mean(s(1+n_i:n_i+n_j,1));
        zm_j_y = mean(s(1+n_i:n_i+n_j,2));
        zm_j_z = mean(s(1+n_i:n_i+n_j,3));
        
        dist_i_j = sqrt((m_i_x-m_j_x)^2 + (m_i_y-m_j_y)^2 + (m_i_z-m_j_z)^2);
        si_dist_i_j = sqrt((zm_i_x-zm_j_x)^2 + (zm_i_y-m_j_y)^2 + (zm_i_z-m_j_z)^2);
        %mahal_dist_i_j = mahal([m_j_x, m_j_y, m_j_z], s(1:n_i,1:3));
        %disp(['Unit ' num2str(i) ' and Unit ' num2str(j)]);
        %disp(['Distance of means ' num2str(dist_i_j)]);
        %disp(['Mahalanobis ' num2str(m_dist_i_j)]);
        %pause;
        example(i,j)= dist_i_j;
        example2_j(j)=si_dist_i_j;
        hold off;
        %
        %waves_j_s_1 = median(s(:,1));
        %waves_j_s_2 = median(s(:,2));
        %waves_j_s_3 = median(s(:,3));

    end
    example2(i,:) = example2_j;
    
end
