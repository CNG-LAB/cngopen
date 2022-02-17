% Fig 4. Structure-function coupling in different macaque samples. 
% A) Anesthetized macaque (sample: Davis) structure-function coupling 
% and correlation with human structure-function coupling projected to 
% macaque space; B) Awake macaque (sample Newcastle) structure-function 
% coupling and correlation with human structure-function coupling projected 
% to macaque space; C) Anesthetized macaque (sample Oxford) structure-
% function coupling and correlation with human structure-function coupling 
% projected to macaque space; D) Structure-function coupling averaged in 
% cytoarchitectural classes 6 across 182 parcels of the Markov parcellation 
% for HCP, Davis, Newcastle and Oxford samples, boxes show the median and 
% interquartile (25-75%) range, whiskers depict the 1.5*IQR from the quartile. 
% Source data are provided as a Source Data file.


% for preprocessing please have a look here: https://github.com/CNG-LAB/PRIME-DE

% for structure-function associations in macaques we used the same approach
% as for humans: 

% correlate both metrics for Davis / Newcastle and Oxford samples
for i = 1:182
    [r p] = corr(fc_mac(i,:)',fc_mac(i,:)','type','spearman')
    fc_mpc_mac(i) = r;
end

mac_map = macaque_template;
for i = 1:91
    mac_map(:,find(markovparcellation==i+1)) = fc_mpc_mac(i);
end
for i = 1:91
    mac_map(:,find(markovparcellation==i+1001)) = fc_mpc_mac(i+91);
end

f = figure,
BoSurfStatViewData(mac_map,surf.macaque,'')
colormap(flipud(cbrewer('div','RdBu',99)))

% scatter
 f=figure,
    scatter(davis_mcc_fc_mpc,fc_mpc_human2macaque_m,100, mesu_mac,'filled'),lsline
    colormap((mes_color))
   
% spin tests - based on brainspace.readthedocs.de
macaque_sample = newc_mcc_fc_mpc;
    for spin_test = 1
        mc_tdiff = zeros(1,length(maskmaskvo.x1));
        for i = 1:91
            mc_tdiff(:,find(mask_markov==i)) = macaque_sample(i);
        end
        for i = 1:91
            mc_tdiff(:,find(mask_markov==i+1000)) = macaque_sample(i+91);
        end
      
        n_permutations = 1000;
        y_rand = spin_permutations({[mc_tdiff(1:32492)'],[mc_tdiff(32493:end)']}, ...
            {LHmSP, RHmSP}, ...
            n_permutations);
        
        % Merge the rotated data into single vectors
        mc_rotated = squeeze([y_rand{1}(:,1,:);y_rand{2}(:,1,:)]);
        
        [r_original_thick, pval_thick_spin] = corr(fc_mpc_human2macaque_m',macaque_sample, ...
            'rows','pairwise','type','spearman')
        
        for i = 1:91
            mc_rotated_p(i,:) = mean(mc_rotated(find(mask_markov==i),:));
        end
        for i = 1:91
            mc_rotated_p(i+91,:) = mean(mc_rotated(find(mask_markov==i+1000),:));
        end      
        
        r_rand_thick = corr(mc_rotated_p,fc_mpc_human2macaque_m', ...
            'rows','pairwise','type','spearman');
        
        % Compute percentile rank.
        prctile_rank_thick = mean(r_original_thick > r_rand_thick)
        
    end
  
    
% spintest
% Raincloud projection
for macaque_rain = 1
    for yeo_1 = 1
        
        to_brain_m = fc_mpc_mac;
        change_p = [];
        change_p.v = to_brain_m(find(yeo_mcc2==1));
        change_p.sn = to_brain_m(find(yeo_mcc2==2));
        change_p.da = to_brain_m(find(yeo_mcc2==3));
        change_p.va = to_brain_m(find(yeo_mcc2==4));
        change_p.l = to_brain_m(find(yeo_mcc2==5));
        change_p.fp = to_brain_m(find(yeo_mcc2==6));
        change_p.dmn = to_brain_m(find(yeo_mcc2==7));
        
        cl = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164;...
            230 148 34; 205 62 78]./256;
      
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d{1} = change_p.v';
        d{2} = change_p.sn';
        d{3} = change_p.da';
        d{4} = change_p.va';
        d{5} = change_p.l';
        d{6} = change_p.fp';
        d{7} = change_p.dmn';
        
        means = cellfun(@mean, d);
        variances = cellfun(@std, d);
        
        f = figure('Position', fig_position);
        h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 0);
        h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
        h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
        h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
            'box_col_match', 0);
        h5 = raincloud_plot(d{5}, 'box_on', 1, 'color', cl(5,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', .95, 'box_col_match', 0);
        h6 = raincloud_plot(d{6}, 'box_on', 1, 'color', cl(6,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', 1.15, 'dot_dodge_amount', 1.15, 'box_col_match', 0);
        h7 = raincloud_plot(d{7}, 'box_on', 1, 'color', cl(7,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', 1.35, 'dot_dodge_amount', 1.35,...
            'box_col_match', 0);
        
        set(gca,'XLim', [-1 1], 'YLim', [-3 3]);
        box off
        exportfigbo(f,[RPATH 'fc-mpc.yeo.mc.rain.png'],'png', 10)
    end
    
    for mesul_1 = 1
        to_brain_m = fc_mpc_mac;
        change_p.pl  = to_brain_m(find(mesul_m==4));
        change_p.he = to_brain_m(find(mesul_m==3));
        change_p.un = to_brain_m(find(mesul_m==2));
        change_p.pr = to_brain_m(find(mesul_m==1));
        
        mes_color = (cbrewer('div','Spectral',5));
        cl = flipud(mes_color([1:2,4:5],:));
        
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d = [];
        d{1} = change_p.pl';
        d{2} = change_p.he';
        d{3} = change_p.um';
        d{4} = change_p.pr';
        
        means = cellfun(@mean, d);
        variances = cellfun(@std, d);
        
        f = figure('Position', fig_position);
        h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 0);
        h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
        h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
        h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75,...
            'box_col_match', 0);
         
        set(gca,'XLim', [-1 1], 'YLim', [-2.5 3]);
        box off
        exportfigbo(f,[RPATH 'fc-mpc.mesul.mc.rain.png'],'png', 10)
    end
end

    