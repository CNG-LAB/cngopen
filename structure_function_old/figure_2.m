
% Fig 2. Genetic basis of structure-function coupling. A) Heritability of 
% MPC and rsFC; B) Row-wise association between mean and heritable seed-wise 
% connectivity – reflecting genetic control over connectivity profiles, 
% Right: distribution of coupling in 400 parcels per cytoarchitectural class, 
% box shows the median and interquartile (25-75%) range, whiskers depict the 
% 1.5*IQR from the quartile; C) Heritability of microstructure-function 
% coupling, Right: distribution of heritability in 400 parcels per 
% cytoarchitectural class, box as in B). Source data are provided 
% as a Source Data file.

% aux data: 
SN = aux_data.surf;
parcels400 = aux_data.parcels400;
mesulam400 = aux_data.mesul;

% similar approach for the heritable components
for i = 1:400
    [r p] = corr(fc400m(i,:)',hr_fc(i,:)','type','spearman')
    fc_fch(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = fc_fch(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = fc_fch(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))


% similar approach for the heritable components
for i = 1:400
    [r p] = corr(MPCm(i,:)',hr_mpc(i,:)','type','spearman')
    mpc_mpch(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = mpc_mpch(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = mpc_mpch(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))


% Plot heritability as computed by solar

heri_ct = zeros(1,20484);
for i = 1:200
    heri_ct(:,find(parcels400==i+1)) =   coupling_her(i,2);
end
for i = 1:200
    heri_ct(:,find(parcels400==i+1001)) = coupling_her(i+200,2);
end


f=figure,
BoSurfStatViewData(heri_ct,SN,'')
colormap((cbrewer('seq','Reds',11)))


% for Rainclooudplots https://github.com/RainCloudPlots/RainCloudPlots

for visuals = 1
    change_p.pl  = coupling_her(find(mesulam400==4),2);
    change_p.hm  = coupling_her(find(mesulam400==3),2);
    change_p.um  = coupling_her(find(mesulam400==2),2);
    change_p.pr  = coupling_her(find(mesulam400==1),2);
    mes_color = flipud(cbrewer('div','Spectral',5));
    cl = mes_color([1:2,4:5],:);
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
    box off
    exportfigbo(f,[RPATH 'F2.heritability.coupling.mesulam.png'],'png', 10)
end

% Computing mean and standard deviations
for i = 1:4
    MN = mean(coupling_her(mesulam400==i,2));
    ST = std(coupling_her(mesulam400==i,2));
    mean_std_Mc(i,1) = MN;
    mean_std_Mc(i,2) = ST;
end
mean_std_Mc

% Difference between networks
for i = 1:4
    for j = 1:4
        [H,P,CI,STATS] = ttest2(coupling_her(mesulam400==i,2),coupling_her(mesulam400==j,2))
        mes_coup_ttest(i,j) = STATS.tstat
        mes_coup_p(i,j) = P;
    end
end










   