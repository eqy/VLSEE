%%Second merge step following get_final_units
n_units = length(spiketimes);
for i=1:n_units;
    stimesi = spiketimes{i};
    isi_i = diff(stimesi) .* 1000;
    isi_i(isi_i>max_isi)=[];
    if(length(isi_i) < min_spikes_for_fit)
        continue
    end
    gamma_i = fitdist(isi_i', 'gamma');
    for j=1:n_units;
       if j <= i
           continue
       end
       stimesj = spiketimes{j};
       isi_j = diff(stimesj) .* 1000;
       isi_j(isi_j>max_isi)=[];
       if(length(isi_j) < min_spikes_for_fit)
           continue
       end
       gamma_j = fitdist(isi_j', 'gamma');
       
       %%Main comparison
       %alpha_diff = abs(gamma_i.Params(1) - gamma_j.Params(1));
       %beta_diff = abs(gamma_i.Params(2) - gamma_j.Params(2));
       
%        if gamma_i.Params(1) > gamma_j.Params(1)
%             alpha_ratio = gamma_i.Params(1)/gamma_j.Params(1);
%        else
%             alpha_ratio = gamma_j.Params(1)/gamma_i.Params(1);
%        end
%        
%        if gamma_i.Params(2) > gamma_j.Params(2)
%             beta_ratio = gamma_i.Params(2)/gamma_j.Params(2);
%        else
%             beta_ratio = gamma_j.Params(2)/gamma_i.Params(2);
%        end
      
       alpha_p_diff = 100*abs(gamma_i.Params(1)-gamma_j.Params(1))/mean([gamma_i.Params(1), gamma_j.Params(1)]);
       beta_p_diff  = 100*abs(gamma_i.Params(2)-gamma_j.Params(2))/mean([gamma_i.Params(2), gamma_j.Params(2)]);
       %if(alpha_diff < max_alpha_diff && beta_diff < max_alpha_diff)
       if(alpha_p_diff < max_alpha_p_diff || beta_p_diff < max_beta_p_diff)
           disp(['alpha percentage difference: ' num2str(alpha_p_diff)])
           disp(['beta percentage difference: ' num2str(beta_p_diff)])
           %favored to merge at this point
       %%Decision already made a this point    
       disp(['unit i: ' num2str(i) ' unit j: ' num2str(j)])
       %disp(['alpha ratio: ' num2str(alpha_ratio) ' beta ratio: ' num2str(beta_ratio)])
       end
       
    end
       
end
