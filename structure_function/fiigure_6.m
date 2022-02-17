
% Fig 6. Difference in organizational gradients of MPC and rsFC in humans 
% and macaques. A) Left panel: Principal MPC and rsFC gradient; Middle 
% panel: alignment; Right panel: delta_rsFC-MPCG1; B) Principal gradient of 
% MPC and tertiary gradient of rsFC in macaques; C) Difference between 
% principal gradients of MPC and rsFC in humans, mapped to macaque space, 
% and difference between corresponding gradients in macaques (lower panel); 
% right: correlation between human and macaque maps; D) Raincloud plots of 
% organizational differences as a function of cytoarchitectural class, 
% and functional networks in humans (400 parcels, Schaefer parcellation) 
% and macaques (182 parcels, Markov parcellation), boxes show the median and 
% interquartile (25-75%) range, whiskers depict the 1.5*IQR from the quartile. 
% Source data are provided as a Source Data file.

% human
gm = GradientMaps('kernel','na','approach','dm','align','pa');
gm = gm.fit({MPCm,fc400m});

mpc2mpc = gm.aligned{1}(:,1)
fc2mpc  = gm.aligned{2}(:,1)

[c1] = rescale(-mpc2mpc)
[f1] = rescale(-fc2mpc)

heri_ct = zeros(1,20484);
for i = 1:200
    heri_ct(:,find(parcels400==i+1)) = c1(i) - f1(i);
end
for i = 1:200
    heri_ct(:,find(parcels400==i+1001)) = c1(i+200)-f1(i+200);
end

f = figure,
BoSurfStatViewData(heri_ct,SN,'')
colormap(cork)
BoSurfStatColLim([-1 1])

% macaque
gm = GradientMaps('kernel','na','approach','dm','align','pa');
gm = gm.fit({mpc_matrix.cov_mean,FC_Markov});

mc_mpc_a_g = gm.gradients{1};
mc_fc_a_g  = gm.gradients{2};

[mc_rank_mpc] =rescale(-mc_mpc_a_g(:,1))
[mc_rank_fc] = rescale(mc_fc_a_g(:,3))

mc_diff = mc_rank_mpc-mc_rank_fc;

heri_ct = zeros(1,length(maskmaskvo.x1));
for i = 1:91
    heri_ct(:,find(mask_markov==i)) = mc_diff(i);
end
for i = 1:91
    heri_ct(:,find(mask_markov==i+1000)) = mc_diff(i+91);
end

f = figure,
BoSurfStatViewData(heri_ct,surf_m,'')
colormap(cork)
BoSurfStatColLim([-1 1])


% raincloudplots analogue to  Figure 2
% Raincloud projection
for human_rain = 1
    for yeo_1 = 1
        
        to_brain_m = c1-f1;
        change_p.v = to_brain_m(find(yeo_post_hoc==1));
        change_p.sn = to_brain_m(find(yeo_post_hoc==2));
        change_p.da = to_brain_m(find(yeo_post_hoc==3));
        change_p.va = to_brain_m(find(yeo_post_hoc==4));
        change_p.l  = to_brain_m(find(yeo_post_hoc==5));
        change_p.fp = to_brain_m(find(yeo_post_hoc==6));
        change_p.dmn = to_brain_m(find(yeo_post_hoc==7));
        
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
        
        set(gca,'XLim', [-1 1], 'YLim', [-2 2]);
        box off
    end
    
    for mesul_1 = 1
        to_brain_m = c1-f1
        change_p.pl  = to_brain_m(find(mesulam400==4));
        change_p.hm = to_brain_m(find(mesulam400==3));
        change_p.um = to_brain_m(find(mesulam400==2));
        change_p.pr = to_brain_m(find(mesulam400==1));
        
        mes_color = (cbrewer('div','Spectral',5));
        cl = flipud(mes_color([1:2,4:5],:));
        
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d = [];
        d{1} = change_p.pl';
        d{2} = change_p.hm';
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
        
        set(gca,'XLim', [-1 1], 'YLim', [-3.5 5]);
        box off
    end
end
for macaque_rain = 1
    for yeo_1 = 1
        
        to_brain_m = mc_diff;
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
    end
    
    for mesul_1 = 1
        to_brain_m = mc_diff;
        change_p.pl  = to_brain_m(find(mesul_m==4));
        change_p.hm = to_brain_m(find(mesul_m==3));
        change_p.um = to_brain_m(find(mesul_m==2));
        change_p.pr = to_brain_m(find(mesul_m==1));
        
        mes_color = (cbrewer('div','Spectral',5));
        cl = flipud(mes_color([1:2,4:5],:));
        
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d = [];
        d{1} = change_p.pl';
        d{2} = change_p.hm';
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
    end
end
