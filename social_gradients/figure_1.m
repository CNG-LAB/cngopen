
%% Main analysis for figure 1

% figure 1c
%load the networks
for load_networks = 1
    att_network(1:10242)     = SurfStatReadData1([datapath, '/DATA/clean_conjunction_fwe05_k10_attention_3mm_fsa5_lh.mgh']);
    att_network(10243:20484) = SurfStatReadData1([datapath, '/DATA/clean_conjunction_fwe05_k10_attention_3mm_fsa5_rh.mgh']);
    
    tom_network(1:10242)     = SurfStatReadData1([datapath, '/DATA/qtom_qntom_3mm_fsa5_lh.mgh']);
    tom_network(10243:20484) = SurfStatReadData1([datapath, '/DATA/qtom_qntom_3mm_fsa5_rh.mgh']);
    
    com_network(1:10242)     = SurfStatReadData1([datapath, '/DATA/vemo-vneu_3mm_fsa5_lh.mgh']);
    com_network(10243:20484) = SurfStatReadData1([datapath, '/DATA/vemo-vneu_3mm_fsa5_rh.mgh']);
end

ntw = zeros(3,20484);
ntw(1,:) = att_network>0;
ntw(2,:) = com_network>0;
ntw(3,:) = tom_network>0;

f = figure,
BoSurfStatViewData(mean(G1(:,keeptest==1),2),SN,'')
colormap(flipud([cbrewer('seq','Greys',99);0,0,0]))
exportfigbo(f,[RPATH 'F1.mean.G1.png'],'png', 10)

% project the data on the vertices


% load the gradients

% visualize the gradients
for gradients_project = 1
    % gradient 1
    keep1 = find(squeeze(G1(2,:))~=0);
    keep3 = mintersect(keep1, find(keepallrest==1));
    keep_bv = keep3;
    %try to regress movement out before difference score
    keeptest = zeros(length(id),1);
    keeptest(keep3) = 1;
    
    f = figure,
    BoSurfStatViewData(mean(G1(:,keeptest==1),2),SN,'')
    colormap(flipud([cbrewer('seq','Greys',99);0,0,0]))
    exportfigbo(f,[RPATH 'F1.mean.G1.png'],'png', 10)
    
    % gradient 2
    keep1 = find(squeeze(G2(:,2,1))~=0);
    keep3 = mintersect(keep1, find(keepallrest==1));
    keep_bv = keep3;
    %try to regress movement out before difference score
    keeptest = zeros(length(id),1);
    keeptest(keep3) = 1;
    
    f = figure,
    BoSurfStatViewData(mean(G2(:,keeptest==1),2),SN,'')
    colormap(flipud([cbrewer('seq','Greys',99);0,0,0]))
    exportfigbo(f,[RPATH 'F1.mean.G2.png'],'png', 10)
    
    % gradient 3
    keep1 = find(squeeze(G3(:,2,1))~=0);
    keep3 = mintersect(keep1, find(keepallrest==1));
    keep_bv = keep3;
    %try to regress movement out before difference score
    keeptest = zeros(length(id),1);
    keeptest(keep3) = 1;
    
    f = figure,
    BoSurfStatViewData(mean(G3(:,keeptest==1),2),SN,'')
    colormap(flipud([cbrewer('seq','Greys',99);0,0,0]))
    exportfigbo(f,[RPATH 'F1.mean.G3.png'],'png', 10)
end

% visualize the eccentricity
for gradient_eccentricity_viz = 1
    f = figure,
    BoSurfStatViewData(mean(GGG(keeptest==1,:)),SN,'')
    colormap([0,0,0;flipud(batlow)])
    exportfigbo(f,[RPATH 'F1.mean.GGG.png'],'png', 10)
    
end

% social task-eccentricity 
for social_eccentricity = 1
    ntw2 = ntw.*mask'
    socialg.a = mean(GGG(keeptest==1,find(ntw2(1,:)>0)));
    socialg.c = mean(GGG(keeptest==1,find(ntw2(2,:)>0)));
    socialg.t = mean(GGG(keeptest==1,find(ntw2(3,:)>0)));
    
    
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
    exportfigbo(f,[RPATH 'F1.bar.GGG.png'],'png', 10)
end

% change in eccentricity
for i = 1
    Xn            = Z_last;
    Xn(isnan(Xn)) = 0;
    Xn(isinf(Xn)) = 0;
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,find(strcmp(groupN,'Presence')))
    keep1   = union(keep1,find(strcmp(groupN,'Control1')))
    keep2   = intersect(find(abs(mean(Xn,2)) <666), find(sum(Xn,2)~=0));
    keep3   = mintersect(find(~strcmp(group4,'Group3')),find(tpnum>0)); %Affect tc3 is excluded
    % in order to not include unbalanced training groups in the main
    % analysis.
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
    Ck(Ck==0)=1;
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
            exportfigbo(f,[RPATH 'F1.change.network',num2str(i), '.png'],'png', 10)
        end
    end
   
    % training versus the two other trainings
    for perspective = 1
        slm      = SurfStatT(slm,-(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
  
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.others-perspective.png'],'png', 10)
        
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
        
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.perspective-others.png'],'png', 10)

        F1.perspective_slm = slm;
    end
       
    for affective = 1
        slm      = SurfStatT(slm,-(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
        
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.others-affect.png'],'png', 10)
        
        
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.025)
        
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.affect-others.png'],'png', 10)

        F1.affect_slm = slm;
    end
    
    for presence = 1
        
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
        
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.presence-others.png'],'png', 10)

        F1.presence_slm = slm;
        
        slm      = SurfStatT(slm,-(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        [pval, peak, clus, clusid] = SurfStatP(slm,mask==1,0.0025)
        
        f = figure,
        BoSurfStatViewData(pval.C,SN,''), SurfStatColLim([0 0.025])
        exportfigbo(f,[RPATH 'F1.others-presence.png'],'png', 10)
        
    end
    
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F1.yeo_change(i,1) = slm.t
        F1.yeo_change(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F1.yeo_change(i,3) = slm.t
        F1.yeo_change(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:7
        slm      = SurfStatLinMod(mean(Ck(:,find(yeo_networks==i)),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F1.yeo_change(i,5) = slm.t
        F1.yeo_change(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end
    % Social networks
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Presence-(0.5*(GN.Affect)+(GN.Perspective))));
        
        F1.sn_change(i,1) = slm.t
        F1.sn_change(i,2) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Affect-(0.5*(GN.Perspective)+(GN.Presence))));
        
        F1.sn_change(i,3) = slm.t
        F1.sn_change(i,4) = (1 - tcdf(slm.t,slm.df))
        
    end
    for i = 1:3
        slm      = SurfStatLinMod(mean(Ck(:,find(ntw(i,:))),2),M);
        slm      = SurfStatT(slm,(GN.Perspective-(0.5*(GN.Affect)+(GN.Presence))));
        
        F1.sn_change(i,5) = slm.t
        F1.sn_change(i,6) = (1 - tcdf(slm.t,slm.df))
        
    end

    save('/Users/sofievalk/Documents/GitHub/micasoft/sandbox/sofie/social_gradients/F1.mat','F1')
    
end




