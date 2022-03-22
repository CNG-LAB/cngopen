
% Fig 3. Microstructure-function coupling in macaques A) Creating MPC in 
% macaques; B) MPC matrix in macaques, ordered along cytoarchitectural 
% class based on 59 and Markov labels 61; C) rsFC matrix in macaques, 
% ordered along cytoarchitectural class; D) Correspondence between MPC 
% and rsFC in macaques; E) Cytoarchitectural classes; F) Functional 
% communities based on 28; G) Row-wise association of MPC and rsFC; 
% upper panel: human map in macaque space; lower panel: macaque map; 
% right: scatter between human and macaque MPC-rsFC coupling; H) 
% Raincloud plots of coupling in humans and macaques as a function of 
% cytoarchitectural class and functional communities 28 in macaque space 
% (182 parcels of Markov parcellation), boxes show the median and 
% interquartile (25-75%) range, whiskers depict the 1.5*IQR from the 
% quartile. Source data are provided as a Source Data file.

% for preprocessing please have a look here: https://github.com/CNG-LAB/PRIME-DE

% for structure-function associations in macaques we used the same approach
% as for humans: 

% aux data: 
yeo_mcc2 = aux_data.yeo_macaque;
mesul_m  = aux_data.mesul_macaque;

% correlate both metrics
for i = 1:182
    [r p] = corr(fc_mac(i,:)',MPC_mac(i,:)','type','spearman')
    fc_mpc_mac(i) = r;
end

mac_map = macaque_template;
for i = 1:91
    mac_map(:,find(aux_data.markovparcellation==i+1)) = fc_mpc_mac(i);
end
for i = 1:91
    mac_map(:,find(aux_data.markovparcellation==i+1001)) = fc_mpc_mac(i+91);
end

f = figure,
BoSurfStatViewData(mac_map,aux_data.surf_macaque,'')
colormap(flipud(cbrewer('div','RdBu',99)))

% then we projected the human data in macaque space (see the prime-de
% github for code and computed the correlation 

corr(fc_mpc_mac,fc_mpc_human2mac,'type','spearman')

f=figure,
scatter(fc_mpc_mac,fc_mpc_human2mac,100,'filled','k'),lsline

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





