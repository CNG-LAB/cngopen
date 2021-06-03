
%gradient 1

for f2_a = 1
    % social task along the gradient (G1)
    for social_g1 = 1
        ntw2 = ntw.*mask'
        socialg.a = mean(G1(keeptest==1,find(ntw2(1,:)>0)));
        socialg.c = mean(G1(keeptest==1,find(ntw2(2,:)>0)));
        socialg.t = mean(G1(keeptest==1,find(ntw2(3,:)>0)));
   
        cl(1, :) = [0.8844    0.7828    0.0195];
        cl(2, :) = [0.9412    0.2314    0.1255];
        cl(3, :) = [0.1922    0.6392    0.3294];
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d{1} = socialg.a';
        d{2} = socialg.c';
        d{3} = socialg.t';
        
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
        legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
        title('Gradient GG');
        set(gca,'XLim', [0 0.12], 'YLim', [-30 42]);
        box off
        exportfigbo(f,[RPATH 'F2.social.G1.png'],'png', 10)
    end
    
    % change in eccentricity
    for i = 1
        Xn            = G1_last;
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0));
        keep    = mintersect(keep1,  keep2,keep3);
        
        Ck      = Xn(keep,:);
        ik      = id(keep,:);
        ink     = idnum(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        GN      = term(gNk);
        Sub     = term(ink);
        
         Ck(Ck==0) = 1;
      
        M        = 1 + A + S + GN + random(Sub) + I;
        
        slm      = SurfStatLinMod(Ck,M,SW);
        
        
        for social_gradient_1 = 1
            for i = 1:3
                
                keep_presence = (find(strcmp(groupN(keep),'Presence')))
                keep_affect = (find(strcmp(groupN(keep),'Affect')))
                keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                keep_control   = (find(strcmp(groupN(keep),'Control1')))
                presence.rcc = mean(Ck(keep_control,find(ntw(i,:)>0)));
                presence.a = mean(Ck(keep_presence,find(ntw(i,:)>0)));
                presence.c = mean(Ck(keep_affect,find(ntw(i,:)>0)));
                presence.t = mean(Ck(keep_perspective,find(ntw(i,:)>0)));
                            
                cl(2, :) = [0.8844    0.7828    0.0195];
                cl(3, :) = [0.9412    0.2314    0.1255];
                cl(4, :) = [0.1922    0.6392    0.3294];
                cl(1, :) = [0.5 0.8 0.9];
                
                fig_position = [200 200 600 400]; % coordinates for figures
                
                d{1} = presence.rcc';
                d{2} = presence.a';
                d{3} = presence.c';
                d{4} = presence.t';
                
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
                    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
                %legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
                %title('Gradient 1');
                set(gca,'XLim', [-0.01 0.01], 'YLim', [-350 400]);
                box off
                exportfigbo(f,[RPATH 'F2.change.G1.network',num2str(i), '.png'],'png', 10)
            end
        end
  
        for perspective = 1
            slm      = SurfStatT(slm,-(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-perspective.g1.png'],'png', 10)
            
            slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.perspective-others.g1.png'],'png', 10)
            
            F2.G1perspective_slm = slm;
        end
        
        
        for affective = 1
            slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.affect-others.g1.png'],'png', 10)
            
            F2.G1affect_slm = slm;
             
            slm      = SurfStatT(slm,-(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-affect.g1.png'],'png', 10)
            
        end
        
        for presence = 1
            
            slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.presence-others.g1.png'],'png', 10)
            
            F2.G1presence_slm = slm;
            
            slm      = SurfStatT(slm,-(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-presence.g1.png'],'png', 10)         
        end
        
    % yeo networks
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.yeo_changeG1(i,1) = slm.t
        F2.yeo_changeG1(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.yeo_changeG1(i,3) = slm.t
        F2.yeo_changeG1(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.yeo_changeG1(i,5) = slm.t
        F2.yeo_changeG1(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
    % Social networks
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.sn_changeG1(i,1) = slm.t
        F2.sn_changeG1(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.sn_changeG1(i,3) = slm.t
        F2.sn_changeG1(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.sn_changeG1(i,5) = slm.t
        F2.sn_changeG1(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
    
        
    end
end

%gradient 2

for f2_b = 1
    % social task G1
    for social_g2 = 1
        ntw2 = ntw.*mask'
        socialg.a = mean(G2(keeptest==1,find(ntw2(1,:)>0)));
        socialg.c = mean(G2(keeptest==1,find(ntw2(2,:)>0)));
        socialg.t = mean(G2(keeptest==1,find(ntw2(3,:)>0)));
   
        cl(1, :) = [0.8844    0.7828    0.0195];
        cl(2, :) = [0.9412    0.2314    0.1255];
        cl(3, :) = [0.1922    0.6392    0.3294];
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d{1} = socialg.a';
        d{2} = socialg.c';
        d{3} = socialg.t';
        
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
        legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
        title('Gradient GG');
        set(gca,'XLim', [0 0.12], 'YLim', [-30 42]);
        box off
        exportfigbo(f,[RPATH 'F2.social.G2.png'],'png', 10)
    end
    
    % change in eccentricity
    for i = 1
        Xn            = G2_last;
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control1')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0));
        keep    = mintersect(keep1,  keep2,keep3);
        
        Ck      = Xn(keep,:);
        ik      = id(keep,:);
        ink     = idnum(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        GN      = term(gNk);
        Sub     = term(ink);
        
        Ck(Ck==0) = 1;
        M        = 1 + A + S + GN + random(Sub) + I;
        
        slm      = SurfStatLinMod(Ck,M,SW);
        
        
        for social_gradient_1 = 1
            for i = 1:3
                
                keep_presence = (find(strcmp(groupN(keep),'Presence')))
                keep_affect = (find(strcmp(groupN(keep),'Affect')))
                keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                keep_control   = (find(strcmp(groupN(keep),'Control')))
                presence.rcc = mean(Ck(keep_control,find(ntw(i,:)>0)));
                presence.a = mean(Ck(keep_presence,find(ntw(i,:)>0)));
                presence.c = mean(Ck(keep_affect,find(ntw(i,:)>0)));
                presence.t = mean(Ck(keep_perspective,find(ntw(i,:)>0)));
                            
                cl(2, :) = [0.8844    0.7828    0.0195];
                cl(3, :) = [0.9412    0.2314    0.1255];
                cl(4, :) = [0.1922    0.6392    0.3294];
                cl(1, :) = [0.5 0.8 0.9];
                
                fig_position = [200 200 600 400]; % coordinates for figures
                
                d{1} = presence.rcc';
                d{2} = presence.a';
                d{3} = presence.c';
                d{4} = presence.t';
                
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
                    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
                %legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
                %title('Gradient 1');
                set(gca,'XLim', [-0.01 0.01], 'YLim', [-350 400]);
                box off
                exportfigbo(f,[RPATH 'F2.change.g2network',num2str(i), '.png'],'png', 10)
            end
        end
  
        for perspective = 1
            slm      = SurfStatT(slm,-(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-perspective.g2.png'],'png', 10)
            
            slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.perspective-others.g2.png'],'png', 10)
            
            F2.G2perspective_slm = slm;
        end
        
        
        for affective = 1
            slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.affect-others.g2.png'],'png', 10)
            
            F2.G2affect_slm = slm;
            
            slm      = SurfStatT(slm,-(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-affect.g2.png'],'png', 10)
            
        end
        
        for presence = 1
            
            slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.presence-others.g2.png'],'png', 10)
            
            F2.G2presence_slm = slm;
            
            slm      = SurfStatT(slm,-(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-presence.g2.png'],'png', 10)
            
        end
        
    % yeo networks
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.yeo_changeG2(i,1) = slm.t
        F2.yeo_changeG2(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.yeo_changeG2(i,3) = slm.t
        F2.yeo_changeG2(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.yeo_changeG2(i,5) = slm.t
        F2.yeo_changeG2(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
    % Social networks
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.sn_changeG2(i,1) = slm.t
        F2.sn_changeG2(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.sn_changeG2(i,3) = slm.t
        F2.sn_changeG2(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.sn_changeG2(i,5) = slm.t
        F2.sn_changeG2(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
        
    end
end

%gradient 3

for f2_c = 1
    % social task G3
    for social_g3 = 1
        ntw2 = ntw.*mask'
        socialg.a = mean(G3(keeptest==1,find(ntw2(1,:)>0)));
        socialg.c = mean(G3(keeptest==1,find(ntw2(2,:)>0)));
        socialg.t = mean(G3(keeptest==1,find(ntw2(3,:)>0)));
   
        cl(1, :) = [0.8844    0.7828    0.0195];
        cl(2, :) = [0.9412    0.2314    0.1255];
        cl(3, :) = [0.1922    0.6392    0.3294];
        
        fig_position = [200 200 600 400]; % coordinates for figures
        
        d{1} = socialg.a';
        d{2} = socialg.c';
        d{3} = socialg.t';
        
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
        legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
        title('Gradient GG');
        set(gca,'XLim', [0 0.12], 'YLim', [-30 42]);
        box off
        exportfigbo(f,[RPATH 'F2.social.G3.png'],'png', 10)
    end
    
    % change in eccentricity
    for i = 1
        Xn            = G3_last;
        Xn(isnan(Xn)) = 0;
        Xn(isinf(Xn)) = 0;
        keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
        keep1   = union(keep1,find(strcmp(groupN,'Presence')))
        keep1   = union(keep1,find(strcmp(groupN,'Control')))
        keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
        keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0));
        keep    = mintersect(keep1,  keep2,keep3);
        
        Ck      = Xn(keep,:);
        ik      = id(keep,:);
        ink     = idnum(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = cellstr(sex(keep));
        
        A       = term(ak);
        S       = term(sk);
        GN      = term(gNk);
        Sub     = term(ink);
        
        %CHk = Days_last(keep);
        Ck(Ck==0) = 1;
        M        = 1 + A + S + GN + random(Sub) + I ;
        
        slm      = SurfStatLinMod(Ck,M,SW);
        
        
        for social_gradient_1 = 1
            for i = 1:3
                
                keep_presence = (find(strcmp(groupN(keep),'Presence')))
                keep_affect = (find(strcmp(groupN(keep),'Affect')))
                keep_perspective = (find(strcmp(groupN(keep),'Perspective')))
                keep_control   = (find(strcmp(groupN(keep),'Control')))
                presence.rcc = mean(Ck(keep_control,find(ntw(i,:)>0)));
                presence.a = mean(Ck(keep_presence,find(ntw(i,:)>0)));
                presence.c = mean(Ck(keep_affect,find(ntw(i,:)>0)));
                presence.t = mean(Ck(keep_perspective,find(ntw(i,:)>0)));
                            
                cl(2, :) = [0.8844    0.7828    0.0195];
                cl(3, :) = [0.9412    0.2314    0.1255];
                cl(4, :) = [0.1922    0.6392    0.3294];
                cl(1, :) = [0.5 0.8 0.9];
                
                fig_position = [200 200 600 400]; % coordinates for figures
                
                d{1} = presence.rcc';
                d{2} = presence.a';
                d{3} = presence.c';
                d{4} = presence.t';
                
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
                    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
                %legend([h1{1} h2{1} h3{1}], {'Attention', 'Affect','ToM'});
                %title('Gradient 1');
                set(gca,'XLim', [-0.01 0.01], 'YLim', [-350 400]);
                box off
                exportfigbo(f,[RPATH 'F2.change.g3network',num2str(i), '.png'],'png', 10)
            end
        end
  
        for perspective = 1
            slm      = SurfStatT(slm,-(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-perspective.g3.png'],'png', 10)
            
            slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.perspective-others.g3.png'],'png', 10)
            
            F2.G3perspective_slm = slm;
            
        end
       
        for affective = 1
            slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.affect-others.g3.png'],'png', 10)
            
            F2.G3affect_slm = slm;
            
            slm      = SurfStatT(slm,-(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-affect.g3.png'],'png', 10)
            
        end
        
        for presence = 1
            
            slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.presence-others.g3.png'],'png', 10)
            
            F2.G2presence_slm = slm;
            
            slm      = SurfStatT(slm,-(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
            [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
            
            f = figure,
            BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
            exportfigbo(f,[RPATH 'F2.others-presence.g3.png'],'png', 10)
            
        end
        
        
         % yeo networks
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.yeo_changeG3(i,1) = slm.t
        F2.yeo_changeG3(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.yeo_changeG3(i,3) = slm.t
        F2.yeo_changeG3(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.yeo_changeG3(i,5) = slm.t
        F2.yeo_changeG3(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
    % Social networks
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F2.sn_changeG3(i,1) = slm.t
        F2.sn_changeG3(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F2.sn_changeG3(i,3) = slm.t
        F2.sn_changeG3(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F2.sn_changeG3(i,5) = slm.t
        F2.sn_changeG3(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
      
    save('/Users/sofievalk/Documents/GitHub/micasoft/sandbox/sofie/social_gradients/F2.mat','F2')
    
    end
end
