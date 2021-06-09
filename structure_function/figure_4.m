% bringing it all together

% two-D projection:

% yeo
yeo_test = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164;...
            230 148 34; 205 62 78];
f =figure,
scatter(c1-f1,fc_mpc,100,yeo_post_hoc,'filled','MarkerEdgeColor',[0 .5 .5]),
colormap(yeo_test./256)
   
% mesulam
mes_color = (cbrewer('div','Spectral',5));
mes_color =[[0 0 0] ;mes_color([1:2,4:5],:)];

f =figure,
scatter(c1-f1,fc_mpc,100,mesulam400,'filled'),
colormap(mes_color)
exportfigbo(f,[RPATH 'F4.scatter.mesulam.png'],'png', 10)

% macaque yeo
f =figure,
scatter(mc_diff,fc_mpc_mac,100,yeo_mcc2,'filled','MarkerEdgeColor',[0 .5 .5]),
colormap(yeo_test./256)

% macaque mesulam
f =figure,
scatter(mc_diff,fc_mpc_mac,100,mesul_m,'filled','MarkerEdgeColor',[0 .5 .5]),
colormap(mes_color)

% dual origin

f=figure, scatter(c1-f1,fc_mpc,100,dual400,'filled','MarkerEdgeColor',[0 .5 .5]), colormap( [flipud(cm3);cm2]),
exportfigbo(f,[RPATH 'F4.scatter.dual.png'],'png', 10)

% cortical expansion
f=figure, scatter(c1-f1,fc_mpc,100,HMS400,'filled','MarkerEdgeColor',[0 .5 .5]), colormap( jet),
exportfigbo(f,[RPATH 'F4.scatter.hms.png'],'png', 10)

% functional decoding
for functional_decoding = 1
    measure=fc_mpc;
    for j = 1:24
        
        L = gifti(['neurosynth_z_values/', files_lh(j).name]);
        R = gifti(['neurosynth_z_values/', files_rh(j).name]);
        funcmap= [L.cdata;R.cdata];
        
        bins1   = [];
        bins1   = quantileranks(measure,20);
        
        heri_ct = zeros(1,64984);
        for i = 1:400
            heri_ct(:,find(HCP400_7.cdata==i)) = bins1(i);
        end
        
        
        for i = 1:20
            tmp(j,i) = mean(funcmap(find(heri_ct==i)));
        end
        
    end
    
    for i = 1:20
        fc_mpc_binned(i) = mean(fc_mpc(bins1==i))
    end
    
    tmpz2= zscore(tmp,[],2)
    tmpz22 = tmpz2;
    tmpz2(tmpz2<0.5) = nan;
    weighted_mean = nanmean(tmpz2.*(fc_mpc_binned),2);
    [~,b,c] = unique(weighted_mean)
    
    f=figure,
    imagesc(tmpz22(b,:),[-4 4])
    colormap(vik)
    colorbar
    exportfigbo(f,[RPATH 'F4.coupling.png'],'png', 10)
    
    
    decode_map=c1-f1;
    for j = 1:24
        
        L = gifti(['/neurosynth_z_values/', files_lh(j).name]);
        R = gifti(['/neurosynth_z_values/', files_rh(j).name]);
        funcmap= [L.cdata;R.cdata];
        
        bins1   = [];
        bins1   = quantileranks(decode_map,20);
        
        heri_ct = zeros(1,64984);
        for i = 1:400
            heri_ct(:,find(HCP400_7.cdata==i)) = bins1(i);
        end
        
        
        for i = 1:20
            tmp2(j,i) = mean(funcmap(find(heri_ct==i)));
        end
    end
    
    
    for i = 1:20
        grd_binned(i) = mean(decode_map(bins1==i))
    end
    
    tmpz= zscore(tmp2,[],2)
    tmpzz = tmpz;
    tmpz(tmpz<0.5) = nan;
    weighted_mean = nanmean(tmpz.*(grd_binned),2);
    [~,d,e] = unique(weighted_mean)
    
    f=figure,
    imagesc(tmpzz(d,:),[-4 4])
    colormap(vik)
    colorbar
    exportfigbo(f,[RPATH 'F4.gradient_diff.png'],'png', 10)
    
    j=0
    for i = 1:24
        j=j+1;
        words{j} = [files_lh(i).name]
        worrs3{j} = words{j}(1:end-32);
        wrd{j} = num2str(j);
    end
    
    eucd= sqrt(c.^2+d.^2)
    f=figure,
    x = e; y = e; scatter(x,y,70,'filled');
    txty = worrs3;
    colormap(flipud(cbrewer('div','Spectral',24)))
    dx = -1; dy = 0.7; % displacement so the text does not overlay the data points
    text(x+dx, y+dy, txty,'FontSize',20);
    exportfigbo(f,[RPATH 'F4.2d.func.png'],'png', 10)
    
end

% genetic decoding
for genetic_decoding = 1
    % use abagen to extract gene lists for the Schaefer 400 parcellation
    % that are consistent across donors
    [~,tot_gene_list400,~] = xlsread('/tot_gene_list.xlsx');
    g = load('/tot_gene_consis_lh400.mat');
    %pcor = fdr_bh(g1.tot_gene_consis_lh(:,2),0.001);
    %sig_idx = find(pcor<0.05);
    sig_idx = find(g.tot_gene_consis_lh(:,1)>0.5);
    sig_gene = tot_gene_list400(sig_idx);
    %sig_gene = tot_gene_list400(pcor==1);
    
    %% coupling
    [~,cp_gene_list,~] = xlsread('coupling.xlsx','coupling');
    
    A = cp_gene_list(:,1);
    B = cp_gene_list(:,2);
    D = sig_gene;
    
    cp_1 = (intersect(B,D));
    cp_2 = (intersect(A,D));
    
    
    [~,g_gene_list,~] = xlsread('/Users/sofievalk/Dropbox/Heritability_gradients/Transcriptome_decoding/redo_updated_1_unal_2_al.xlsx','aligned');
    
    E = g_gene_list(:,1);
    F = g_gene_list(:,2);
    D = sig_gene;
    
    align_1  = (intersect(E,D));
    align_2  = (intersect(F,D));
    
    %uncoupling
    char(intersect(setdiff(cp_2,align_2),cp_2))
    %mpc
    char(intersect(setdiff(align_2,cp_2),cp_2))
    
    %coupling
    char(intersect(setdiff(cp_1,align_1),cp_1))
    %rsfc
    char(intersect(setdiff(align_1,cp_1),align_1))
    
end
