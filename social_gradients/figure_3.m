
% First we averaged the values within 400 parcels to optimize bootstrapping

for i = 1:200
    G1_400(i,:) = mean(G1_last(:,parcels400==i+1),2)';
    G2_400(i,:) = mean(G2_last(:,parcels400==i+1),2)';
    G3_400(i,:) = mean(G3_last(:,parcels400==i+1),2)';
    Z_400(i,:)  = mean(Z_last(:,parcels400==i+1),2)';
end
for i = 1:200
    G1_400(i+200,:) = mean(G1_last(:,parcels400==i+1001),2)';
    G2_400(i+200,:) = mean(G2_last(:,parcels400==i+1001),2)';
    G3_400(i+200,:) = mean(G3_last(:,parcels400==i+1001),2)';
    Z_400(i+200,:)  = mean(Z_last(:,parcels400==i+1001),2)';
end



for change_behavioral_tom = 1
    
    val     = (Tom);
    Xn      = [ G1_400', G2_400', G3_400'];
    Xn(isnan(Xn)) = 0;
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,find(strcmp(groupN,'Presence')))
    keep1   = union(keep1,find(strcmp(groupN,'Control1')))
    keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
    keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0)); %Affect tc3 is excluded
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep6   = find(val>-666);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6);
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    gn      = groupN(keep);
    sub     = id(keep);
    A       = term(ak);
    S       = term(sk);
    GN      = term(gn);
    Sub     = term(sub);
    
    M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Perspective - (0.5*(GN.Affect+GN.Presence)))
    pp  = 1 - tcdf(slm.t,slm.df) 
    
     M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Perspective - (GN.Control1))
    pp  = 1 - tcdf(slm.t,slm.df) 
end

for change_behavioral_compassion = 1
    
    val     = (Comp);
    Xn      = [ G1_400', G2_400', G3_400'];
    Xn(isnan(Xn)) = 0;
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,find(strcmp(groupN,'Presence')))
    keep1   = union(keep1,find(strcmp(groupN,'Control1')))
    keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
    keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0)); %Affect tc3 is excluded
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep6   = find(val>-666);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6);
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    gn      = groupN(keep);
    sub     = id(keep);
    A       = term(ak);
    S       = term(sk);
    GN      = term(gn);
    Sub     = term(sub);
    
    M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Affect - (0.5*(GN.Perspective+GN.Presence)))
    pp  = 1 - tcdf(slm.t,slm.df) 
    
     M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Affect - (GN.Control1))
    pp  = 1 - tcdf(slm.t,slm.df) 
end

 for change_behavioral_attention = 1
    
    val     = (Att);
    Xn      = [ G1_400', G2_400', G3_400'];
    Xn(isnan(Xn)) = 0;
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,find(strcmp(groupN,'Presence')))
    keep1   = union(keep1,find(strcmp(groupN,'Control1')))
    keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
    keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0)); %Affect tc3 is excluded
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep7   = find(abs(val)<30);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6,keep7);
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    gn      = groupN(keep);
    sub     = id(keep);
    A       = term(ak);
    S       = term(sk);
    GN      = term(gn);
    Sub     = term(sub);
    
    M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Presence - (0.5*(GN.Perspective+GN.Affect)))
    pp  = 1 - tcdf(slm.t,slm.df) 
    
     M = 1 + A + S + GN + random(Sub) + I;
    slm = SurfStatLinMod(val(keep), M);
    slm = SurfStatT(slm,GN.Presence - (GN.Control1))
    pp  = 1 - tcdf(slm.t,slm.df) 
end


for prediction_tom = 1
    
    val     = (Tom);
    Xn      = [ G1_400', G2_400', G3_400'];
    Xn(isnan(Xn)) = 0;
    keep1   = find(strcmp(groupN,'Perspective'));
    keep2   = intersect(find(abs(mean(Xn,2)) <100), find(sum(Xn,2)~=0));
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6);
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    
    A       = term(ak);
    S       = term(sk);
    
    algoname = 'cv_lassopcr'; % cross-validated penalized regression. Predict pain ratings
    
    M = 1 + A + S;
    slm = SurfStatLinMod(Ck, M);
    feature_matrix = (Ck - slm.X*slm.coef);
    
    dat = dat;
    dat.dat = (feature_matrix)';
    dat.Y =  quantileranks(val(keep),10);
    [cverr, stats, optout] = predict(dat,  'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', 0.5,...
        'bootweights','bootsamples', 5000);
    F3.prediction_stats_perspecitve = stats
    
    f= figure,
    plot_correlation_samefig(stats.Y, stats.yfit);
    xlabel('original outcome (y)');
    ylabel('predicted outcome (y-fit)');
    saveas(f, [RPATH 'tom._p_lasso.5fold.nobin.png'])
    
    to_brain_w = zeros(20484,3)
    for j = 1:3
        for i = 1:200
            to_brain_w(parcels400==i+1,j) = stats.WTS.wZ(i+((j-1)*400));
        end
        for i = 1:200
            to_brain_w(parcels400==i+1001,j) = stats.WTS.wZ((i+200)+((j-1)*400));
        end
    end
    
    f = figure,
    BoSurfStatViewData(mean((to_brain_w),2),SN,'')
    colormap([flipud(cbrewer('div','RdBu',11))])
    BoSurfStatColLim([-2 2])
    saveas(f, [RPATH 'tom._p_lasso.5fold.brainZabs.nobin.png'])
    
    for social_gradient_1 = 1
        
        for i = 1:3
            to_brain_m = (to_brain_w(:,i))
            change_p.a = to_brain_m(find(ntw(1,:)>0));
            change_p.c = to_brain_m(find(ntw(2,:)>0));
            change_p.t = to_brain_m(find(ntw(3,:)>0));
            
            
            
            cl(1, :) = [0.8844    0.7828    0.0195];
            cl(2, :) = [0.9412    0.2314    0.1255];
            cl(3, :) = [0.1922    0.6392    0.3294];
            
            fig_position = [200 200 600 400]; % coordinates for figures
            
            d{1} = change_p.a';
            d{2} = change_p.c';
            d{3} = change_p.t';
            
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
            set(gca,'XLim', [-3 3], 'YLim', [-1 1]);
            box off
            exportfigbo(f,[RPATH 'F3.pred.tom',num2str(i), '.png'],'png', 10)
        end
    end
    
end

for prediction_compassion = 1
    
    
    val     = (Comp);
    
    
    Xn      = [ G1_400', G2_400', G3_400'];
    
    
    Xn(isnan(Xn)) = 0;
    keep1   = find(strcmp(groupN,'Affect'));
    keep2   = intersect(find(abs(mean(Xn,2)) <100), find(sum(Xn,2)~=0));
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep6   = find(abs(val)<30);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6);
    
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    
    A       = term(ak);
    S       = term(sk);
    
    algoname = 'cv_lassopcr'; % cross-validated penalized regression. Predict pain ratings
    
    M = 1 + A + S;
    slm = SurfStatLinMod(Ck, M);
    feature_matrix = (Ck - slm.X*slm.coef);
    
    dat = dat;
    dat.dat = (feature_matrix)';
    dat.Y =  quantileranks(val(keep),10);
    [cverr, stats, optout] = predict(dat,  'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', 0.5,...
        'bootweights','bootsamples', 5000);
    F3.prediction_stats_affect = stats;
    
    f= figure,
    plot_correlation_samefig(stats.Y, stats.yfit);
    xlabel('original outcome (y)');
    ylabel('predicted outcome (y-fit)');
    saveas(f, [RPATH 'com._p_lasso.5fold.nobin.png'])
    
    
    to_brain_w = zeros(20484,3)
    for j = 1:3
        for i = 1:200
            to_brain_w(parcels400==i+1,j) = stats.WTS.wZ(i+((j-1)*400));
        end
        for i = 1:200
            to_brain_w(parcels400==i+1001,j) = stats.WTS.wZ((i+200)+((j-1)*400));
        end
    end
    
    f = figure,
    BoSurfStatViewData(mean((to_brain_w),2),SN,'')
    colormap([flipud(cbrewer('div','RdBu',11))])
    BoSurfStatColLim([-2 2])
    saveas(f, [RPATH 'com._p_lasso.5fold.brainZabs.nobin.png'])
    
    
    for social_gradient_1 = 1
        
        for i = 1:3
            to_brain_m = (to_brain_w(:,i))
            change_p.a = to_brain_m(find(ntw(1,:)>0));
            change_p.c = to_brain_m(find(ntw(2,:)>0));
            change_p.t = to_brain_m(find(ntw(3,:)>0));
            
            
            
            cl(1, :) = [0.8844    0.7828    0.0195];
            cl(2, :) = [0.9412    0.2314    0.1255];
            cl(3, :) = [0.1922    0.6392    0.3294];
            
            fig_position = [200 200 600 400]; % coordinates for figures
            
            d{1} = change_p.a';
            d{2} = change_p.c';
            d{3} = change_p.t';
            
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
            set(gca,'XLim', [-3 3], 'YLim', [-1 1.5]);
            box off
            exportfigbo(f,[RPATH 'F3.pred.com',num2str(i), '.png'],'png', 10)
        end
    end
    
end

for prediction_attention = 1
    
    
    val     = (Att);
    Xn      = [ G1_400', G2_400', G3_400'];
    
    Xn(isnan(Xn)) = 0;
    keep1   = find(strcmp(groupN,'Presence'));
    keep2   = intersect(find(abs(mean(Xn,2)) <100), find(sum(Xn,2)~=0));
    keep4   = find(~(strcmp(group4,'Group3')));
    keep5   = find(tpnum>0);
    keep6   = find(abs(val)<30);
    keep    = mintersect(keep1,keep2,keep4,keep5,keep6);
    
    Ck      = Xn(keep,:);
    ak      = age(keep);
    sk      = cellstr(sex(keep));
    
    A       = term(ak);
    S       = term(sk);
    
    M = 1 + A + S ;
    slm = SurfStatLinMod(Ck, M);
    feature_matrix = (Ck - slm.X*slm.coef);
    
    dat = dat; %dummy fmri_data for input
    dat.dat = (feature_matrix)';
    dat.Y =  quantileranks(val(keep),10);
    [cverr, stats, optout] = predict(dat,  'algorithm_name', 'cv_lassopcr', 'nfolds', 5, 'error_type', 'mse', 'numcomponents', 10, 'Alpha', 0.5,...
        'bootweights','bootsamples', 5000);
    F3.prediction_stats_presence = stats;
    
    f= figure,
    plot_correlation_samefig(stats.Y, stats.yfit);
    xlabel('original outcome (y)');
    ylabel('predicted outcome (y-fit)');
    saveas(f, [RPATH 'att._p_lasso.5fold.10bins.png'])
    
    
    to_brain_w = zeros(20484,3)
    for j = 1:3
        for i = 1:200
            to_brain_w(parcels400==i+1,j) = stats.WTS.wZ(i+((j-1)*400));
        end
        for i = 1:200
            to_brain_w(parcels400==i+1001,j) = stats.WTS.wZ((i+200)+((j-1)*400));
        end
    end
    
    f = figure,
    BoSurfStatViewData(mean((to_brain_w),2),SN,'')
    colormap([flipud(cbrewer('div','RdBu',11))])
    BoSurfStatColLim([-2 2])
    saveas(f, [RPATH 'att._p_lasso.5fold.brainZabs.png'])
    
    for social_gradient_1 = 1
        
        for i = 1:3
            to_brain_m = (to_brain_w(:,i))
            change_p.a = to_brain_m(find(ntw(1,:)>0));
            change_p.c = to_brain_m(find(ntw(2,:)>0));
            change_p.t = to_brain_m(find(ntw(3,:)>0));
            
            
            
            cl(1, :) = [0.8844    0.7828    0.0195];
            cl(2, :) = [0.9412    0.2314    0.1255];
            cl(3, :) = [0.1922    0.6392    0.3294];
            
            fig_position = [200 200 600 400]; % coordinates for figures
            
            d{1} = change_p.a';
            d{2} = change_p.c';
            d{3} = change_p.t';
            
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
            set(gca,'XLim', [-3 3], 'YLim', [-1 1.5]);
            box off
            exportfigbo(f,[RPATH 'F3.pred.att',num2str(i), '.png'],'png', 10)
        end
    end
    
end

 save('/Users/sofievalk/Documents/GitHub/micasoft/sandbox/sofie/social_gradients/F3.mat','F3')
   
