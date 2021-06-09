% gradient construction
% for more details see also brainspace.readthedocs.io

% gradient in humans and their differences: 
for heritable_gradients_mpc = 1
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({MPCm, hr_mpc});
    
    mpc_g1 = gm.aligned{1}
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =mpc_g1(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = mpc_g1(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(batlow)
  
     mpc_h = rescale(gm.aligned{2}(:,1))
     
     heri_ct = zeros(1,20484);
     for i = 1:200
         heri_ct(:,find(parcels400==i+1)) =mpc_h(i)-mpc_g1(i);
     end
     for i = 1:200
         heri_ct(:,find(parcels400==i+1001)) = mpc_h(i+200);
     end
     
     f = figure,
     BoSurfStatViewData(-heri_ct,SN,'')
     colormap(bilbao)
     
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1),'filled','k'),lsline
    xlim([-0.15 0.21])
    
    % heritability along principal gradient
    bins1   = [];
    bins1   = quantileranks(-gm.aligned{1}(:,1),10);
  
    for j=1:10
        for i = 1:10
            mean_herit_m(j,i) = mean(mean(MPCherit(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_herit_m)
    colormap((cbrewer('seq','Reds',99)))
    colorbar
  
end

for heritable_gradients_fc = 1
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({fc400m, hr_fc});
    
    mpc_g1 = gm.aligned{1}
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =mpc_g1(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = mpc_g1(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(batlow)
  
     mpc_h = rescale(gm.aligned{2}(:,1))
     
     heri_ct = zeros(1,20484);
     for i = 1:200
         heri_ct(:,find(parcels400==i+1)) =mpc_h(i)-mpc_g1(i);
     end
     for i = 1:200
         heri_ct(:,find(parcels400==i+1001)) = mpc_h(i+200);
     end
     
     f = figure,
     BoSurfStatViewData(-heri_ct,SN,'')
     colormap(bilbao)
     
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1),'filled','k'),lsline
    xlim([-0.15 0.21])
    
    % heritability along principal gradient
    bins1   = [];
    bins1   = quantileranks(-gm.aligned{1}(:,1),10);
  
    for j=1:10
        for i = 1:10
            mean_herit_m(j,i) = mean(mean(fc400herit(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_herit_m)
    colormap((cbrewer('seq','Reds',99)))
    colorbar
end


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
        change_p.l = to_brain_m(find(yeo_post_hoc==5));
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
        exportfigbo(f,[RPATH 'fc-mpc.yeo.rain1.png'],'png', 10)
    end
    
    for mesul_1 = 1
        to_brain_m = c1-f1
        change_p.v  = to_brain_m(find(mesulam400==4));
        change_p.sn = to_brain_m(find(mesulam400==3));
        change_p.da = to_brain_m(find(mesulam400==2));
        change_p.va = to_brain_m(find(mesulam400==1));
        
        mes_color = (cbrewer('div','Spectral',5));
        cl = flipud(mes_color([1:2,4:5],:));
        
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d = [];
        d{1} = change_p.v';
        d{2} = change_p.sn';
        d{3} = change_p.da';
        d{4} = change_p.va';
        
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
        exportfigbo(f,[RPATH 'gradients.mesul.rain.png'],'png', 10)
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
        exportfigbo(f,[RPATH 'gradients.yeo.mc.rain.png'],'png', 10)
    end
    
    for mesul_1 = 1
        to_brain_m = mc_diff;
        change_p.v  = to_brain_m(find(mesul_m==4));
        change_p.sn = to_brain_m(find(mesul_m==3));
        change_p.da = to_brain_m(find(mesul_m==2));
        change_p.va = to_brain_m(find(mesul_m==1));
        
        mes_color = (cbrewer('div','Spectral',5));
        cl = flipud(mes_color([1:2,4:5],:));
        
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d = [];
        d{1} = change_p.v';
        d{2} = change_p.sn';
        d{3} = change_p.da';
        d{4} = change_p.va';
        
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
        exportfigbo(f,[RPATH 'gradients.mesul.mc.rain.png'],'png', 10)
    end
end






