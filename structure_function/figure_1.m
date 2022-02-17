% Fig 1. Structure-function coupling and heritability in human cortical 
% regions. A) Microstructural profile covariance (MPC) was chosen to map 
% networks of microstructural similarity for each cortical node, sorted 
% along cytoarchitectural class 22,59; B) Resting-state functional connectivity 
% (rsFC) analysis maps nodal patterns of intrinsic functional connectivity, 
% sorted along cytoarchitectural class; C) Row-wise coupling of MPC and 
% rsFC Middle: raincloud plot of distribution within cytoarchitectural 
% classes of coupling in the 400 Schaefer parcels, box shows the median 
% and interquartile (25-75%) range, whiskers depict the 1.5*interquartile 
% range (IQR) from the quartiles; Right: Reference visualization 
% cytoarchitectural class. Source data are provided as a Source Data file.


% compute the correlation between microstructure profile covariance and
% functional connectivity

% aux data: 
SN = aux_data.surf;
parcels400 = aux_data.parcels400;
mesulam400 = aux_data.mesul;

% make the mean maps of both metrics
fc400m = squeeze(mean(fc400z(keep,:,:),1));
MPCm   = squeeze(mean(MPC(keep,:,:),1));

% correlate both metrics
for i = 1:400
    [r p] = corr(fc400m(i,:)',MPCm(i,:)','type','spearman')
    fc_mpc(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = fc_mpc(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = fc_mpc(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))

% for Rainclooudplots https://github.com/RainCloudPlots/RainCloudPlots

for visuals = 1
    change_p.pl  = fc_mpc(find(mesulam400==4),2);
    change_p.hm  = fc_mpc(find(mesulam400==3),2);
    change_p.um  = fc_mpc(find(mesulam400==2),2);
    change_p.pr  = fc_mpc(find(mesulam400==1),2);
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
end


