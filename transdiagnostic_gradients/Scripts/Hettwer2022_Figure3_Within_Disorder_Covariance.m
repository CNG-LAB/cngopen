%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coordinated Cortical Thickness alterations Across Mental Disorders: A Transdiagnostic ENIGMA Project
% Script 3: Computes disorder-specific covariance profiles and embeds
% individual disorders within a transdiagnostic space.

%This script runs analyses presented in Figure 3 of the manuscript.
%Meike Hettwer 2022 - CNG Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add paths and Toolboxes
%ENIGMA-1.1.3huch
%cbrewer


%data to load from GitHub:

%Epicenters.mat
%CT_psych_gradients.mat
%corr_dis.mat

%% Load ENIGMA summary statistics
%Schizophrenia
sum_stats = load_summary_stats('schizophrenia');
CT = sum_stats.CortThick_case_vs_controls;
CT_d_sz = CT.d_icv;

%Bipolar
sum_stats = load_summary_stats('bipolar');
CT = sum_stats.CortThick_case_vs_controls_adult;
CT_d_bipo = CT.d_icv;

%ASD
sum_stats = load_summary_stats('asd');
CT = sum_stats.CortThick_case_vs_controls_meta_analysis;
CT_d_asd = CT.d_icv;

%ADHD
sum_stats = load_summary_stats('adhd');
CT = sum_stats.CortThick_case_vs_controls_adult;
CT_d_adhd = CT.d_icv;

%MDD
sum_stats = load_summary_stats('depression');
CT = sum_stats.CortThick_case_vs_controls_adult;
CT_d_mdd = CT.d_icv;

%OCD
sum_stats = load_summary_stats('ocd');
CT = sum_stats.CortThick_case_vs_controls_adult;
CT_d_ocd = CT.d_icv;

%combine
disorders = [CT_d_bipo, CT_d_adhd, CT_d_asd, CT_d_mdd, CT_d_ocd, CT_d_sz];
disorder_names = {'BD'; 'ADHD'; 'ASD'; 'MDD'; 'OCD'; 'SCZ'};

%% Disorder-specific covariance (*Figure 3B*)

% via differences in Cohen's d values between regions
for i = 1:6
dis = disorders(:,i);
dis_diff=nan(length(dis),length(dis));
for j = 1:length(dis)%columns
    for k = 1:length(dis)%rows
    dis_diff(k,j) = dis(j)-dis(k); %caution: resulting matrix is asymmetrical as signs flip between rows and columns
    end    
end

dis_diff_all(:,:,i) = abs(dis_diff(:,:));%take abs so it's symmetric

clear dis_diff dis dis_diff_threshbin
end

%within-disorder covariance plots
color = cbrewer('seq','BuPu',50);
color(color<0) =0;
clims = [-0.2 0];
for i = 1:6
 this_dis = dis_diff_all(:,:,i)*-1;
 fig=figure; imagesc(this_dis, clims); title(disorder_names{i}); colormap(color);colorbar
 yticks([]), xticks([]); axis square
 %close all
end

%correlation between disorder-specific and transdiagnostic covariance plotted on the brain
for i = 1:6
    this_dis = dis_diff_all(:,:,i);
    for column = 1:68
    r_temp(column,i) = corr(this_dis(:,column)*-1,corr_dis(:,column));
    end
    
   fig=figure; plot_cortical(parcel_to_surface(r_temp(:,i)),'label_text',disorder_names{i},'color_range', [-0.6 0.6])
end
%% Cross-disorder similarity matrix (*Figure 3C*)

cross_disorder_matrix = corr(disorders);
figure; imagesc(cross_disorder_matrix); set(gca, 'XTick',[1:6],'XTickLabel',disorder_names','YTick',[1:6],'YTickLabel',disorder_names,'clim',[-1 1]);colormap(flipud(cbrewer('div','RdBu',25)));colorbar
%% cluster disorders (*Figure 3D*)

tree = linkage(corr(disorders)' ,'average');
figure()
dendrogram(tree,'Labels',disorder_names);

%% Embedding of individual disorders in 2D space framed by transdiagnostic epicenters and hubs
%A) epicenters
%load HCP connectivity data

[fc_ctx, fc_ctx_labels, ~, ~] = load_fc();
[sc_ctx, sc_ctx_labels, ~, ~] = load_sc();
[~, ~, fc_sctx, fc_sctx_labels] = load_fc();
[~, ~, sc_sctx, sc_sctx_labels] = load_sc();

%load Epicenters.mat from Github !!
%cortical
epi_fc_trans = fc_ctx_trans_p_sig;
epi_fc_trans_bin = find(fc_ctx_trans_p_sig>0);

epi_sc_trans = sc_ctx_trans_p_sig;
epi_sc_trans_bin = find(sc_ctx_trans_p_sig>0);
%subcortical
epi_fc_sctx_trans = fc_sctx_epi_p_sig;
epi_fc_sctx_trans_bin = find(fc_sctx_epi_p_sig>0);

epi_sc_sctx_trans = sc_sctx_epi_p_sig;
epi_sc_sctx_trans_bin = find(sc_sctx_epi_p_sig>0);

epicenter_overlap_fc_sctx = nan(1,6);
epicenter_overlap_sc_sctx = nan(1,6);
epicenter_overlap_fc_percent = nan(1,6);
epicenter_overlap_sc_percent = nan(1,6);
epicenter_overlap_fc_sctx = nan(1,6);
epicenter_overlap_sc_sctx = nan(1,6);

for i = 1:size(disorders,2)
      
% Thickness maps for each disorder
fc_ctx_epi              = zeros(size(fc_ctx, 1), 1);
fc_ctx_epi_p            = zeros(size(fc_ctx, 1), 1);
for seed = 1:size(fc_ctx, 1)
    seed_conn           = fc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn,  abs(disorders(:,i)));
    fc_ctx_epi(seed)    = r_tmp(1, 2);
    fc_ctx_epi_p(seed)  = spin_test(seed_conn,  abs(disorders(:,i)), 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Identify cortical epicenter values (from structural connectivity)
sc_ctx_epi              = zeros(size(sc_ctx, 1), 1);
sc_ctx_epi_p            = zeros(size(sc_ctx, 1), 1);
for seed = 1:size(sc_ctx, 1)
    seed_conn           = sc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn,  abs(disorders(:,i)));
    sc_ctx_epi(seed)    = r_tmp(1, 2);
    sc_ctx_epi_p(seed)  = spin_test(seed_conn,  abs(disorders(:,i)), 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Select only regions with p < 0.05 (functional epicenters)
fc_ctx_epi_p_sig = zeros(length(fc_ctx_epi_p), 1);
fc_ctx_epi_p_sig(fc_ctx_epi_p < 0.05) = fc_ctx_epi(fc_ctx_epi_p<0.05);

% Selecting only regions with p < 0.05 (structural epicenters)
sc_ctx_epi_p_sig = zeros(length(sc_ctx_epi_p), 1);
sc_ctx_epi_p_sig(sc_ctx_epi_p < 0.05) = sc_ctx_epi(sc_ctx_epi_p<0.05);
%$$$$$$$$$$
%subcortical

% Epicenter mapping - subcortical seeds but cortical map!
% Identify subcortical epicenters (from functional connectivity)
fc_sctx_epi             = zeros(size(fc_sctx, 1), 1);
fc_sctx_epi_p           = zeros(size(fc_sctx, 1), 1);
for seed = 1:size(fc_sctx, 1)
    seed_conn           = fc_sctx(seed, :);
    r_tmp               = corrcoef(seed_conn, abs(disorders(:,i)));
    fc_sctx_epi(seed)   = r_tmp(1, 2);
    fc_sctx_epi_p(seed) = spin_test(seed_conn, abs(disorders(:,i)), 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Identify subcortical epicenters (from structural connectivity)
sc_sctx_epi             = zeros(size(sc_sctx, 1), 1);
sc_sctx_epi_p           = zeros(size(sc_sctx, 1), 1);
for seed = 1:size(sc_sctx, 1)
    seed_conn           = sc_sctx(seed, :);
    r_tmp               = corrcoef(seed_conn, abs(disorders(:,i)));
    sc_sctx_epi(seed)   = r_tmp(1, 2);
    sc_sctx_epi_p(seed) = spin_test(seed_conn, abs(disorders(:,i)), 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Threshold and Project the results on the surface brain
% functional epicenters
fc_sctx_epi_p_sig = zeros(length(fc_sctx_epi_p), 1);
fc_sctx_epi_p_sig(fc_sctx_epi_p < 0.05) = fc_sctx_epi(fc_sctx_epi_p<0.05);

% structural epicenters
sc_sctx_epi_p_sig = zeros(length(sc_sctx_epi_p), 1);
sc_sctx_epi_p_sig(sc_sctx_epi_p < 0.05) = sc_sctx_epi(sc_sctx_epi_p<0.05);
                
%$$$$$$$$$$
% check overlap with transdiagnostic epicenters
%cortical
fc_ctx_epi_p_sig(fc_ctx_epi_p_sig<0)=0;
sc_ctx_epi_p_sig(sc_ctx_epi_p_sig<0)=0;
functional_epicenters(:,i)= fc_ctx_epi_p_sig;
structural_epicenters(:,i)= sc_ctx_epi_p_sig;     
%sub-cortical
fc_sctx_epi_p_sig(fc_sctx_epi_p_sig<0)=0;
sc_sctx_epi_p_sig(sc_sctx_epi_p_sig<0)=0;
sub_functional_epicenters(:,i)= fc_sctx_epi_p_sig;
sub_structural_epicenters(:,i)= sc_sctx_epi_p_sig;  

%check disorder_specific overlap with transdiagnostic epicenters
%cortical
position_fc = find(fc_ctx_epi_p_sig>0);
position_sc = find(sc_ctx_epi_p_sig>0);
[val_fc,pos_fc]=intersect(epi_fc_trans_bin,position_fc);
[val_sc,pos_sc]=intersect(epi_sc_trans_bin,position_sc);
%sub-cortical
position_fc_sctx = find(fc_sctx_epi_p_sig>0);
position_sc_sctx = find(sc_sctx_epi_p_sig>0);

[val_fc,pos_fc]=intersect(epi_fc_sctx_trans_bin,position_fc_sctx);
[val_sc,pos_sc]=intersect(epi_sc_sctx_trans_bin,position_sc_sctx);

epicenter_overlap_fc(i) = size(intersect(epi_fc_trans_bin,position_fc),1);
epicenter_overlap_sc(i) = size(intersect(epi_sc_trans_bin,position_sc),1);

epicenter_overlap_fc_sctx(i) = size(intersect(epi_fc_sctx_trans_bin,position_fc_sctx),1);
epicenter_overlap_sc_sctx(i) = size(intersect(epi_sc_sctx_trans_bin,position_sc_sctx),1);

epicenter_overlap_fc_percent(i) = ((epicenter_overlap_fc(i)+epicenter_overlap_fc_sctx(i))/(size(epi_fc_trans_bin,1)+size(epi_fc_sctx_trans_bin,1)))*100;
epicenter_overlap_sc_percent(i) = ((epicenter_overlap_sc(i)+epicenter_overlap_sc_sctx(i))/(size(epi_sc_trans_bin,1)+size(epi_sc_sctx_trans_bin,1)))*100;

end

T_epi_overlap = table(epicenter_overlap_fc_percent', epicenter_overlap_sc_percent','RowNames',disorder_names,'VariableNames',{'Functional','Structural'});

% plot disorder_specific epicenters
%functional
for i = 1:6
fig = figure; plot_subcortical(sub_functional_epicenters(:,i),'color_range',[0 0.5],'label_text',disorder_names{i}, 'ventricles', 'False');
fig = figure; plot_cortical(parcel_to_surface(functional_epicenters(:,i)),'color_range',[0 0.5],'label_text', disorder_names{i})
end

%structural
for i = 1:6
figure; plot_subcortical(sub_structural_epicenters(:,i),'color_range',[0 0.5],'label_text',disorder_names{i}, 'ventricles', 'False');
figure; plot_cortical(parcel_to_surface(structural_epicenters(:,i)),'color_range',[0 0.5],'label_text', disorder_names{i})
end
                               
%% Hubs vs disorder-specific maps

% Set up thresholded disorder correlation matrix and dc map
corr_dis_bin_thresh = bsxfun(@gt, corr_dis, prctile(corr_dis,80))'; 

% plot degree centrality map / cross-disorder covariance hub map 
dis_dc = sum(corr_dis_bin_thresh);

r=nan(1,6);
p=nan(1,6);
for i = 1:6
    r(i) = corr(dis_dc', abs(disorders(:,i))); %CHANGED TO ABSOLUTE IMPACT!
    p(i) = spin_test(abs(disorders(:,i)),  dis_dc', 'surface_name', 'fsa5', 'parcellation_name', ...
                                   'aparc', 'n_rot', 1000, 'type', 'pearson');
end
T_transhubs = table(r',p','RowNames',disorder_names,'VariableNames',{'r','p_spin'});
%% Correlate disorder maps with Gradients, 2D plot
Grad1_dis_cor=nan(1,6);
Grad2_dis_cor=nan(1,6);
for i = 1:size(disorders,2)
  Grad1_dis_cor(i) = corr(disorders(:,i),CT_psych_gradient1); 
  Grad2_dis_cor(i) = corr(disorders(:,i),CT_psych_gradient2);   
end
%% Plot how disorder effects are driven by transdiagnostic hubs vs epicenters (*Figure 3E*)

%col = cbrewer('seq','BuPu','50')
txty = {'BD', 'ADHD', 'ASD','MDD', 'OCD', 'SCZ'};
dx = 0.07;
dy = 0.001; 

figure;
subplot(1,2,1)

scatter(T_transhubs.r,T_epi_overlap.Functional,100, [0.5059    0.0588    0.4863],'filled');
set(gca,'xlim', [-1 1],'ylim', [0 100], 'XTick', [-1, -0.5 0 0.5 1],'YTick',[0 50 100],'TickDir','none','color','none','FontSize',10, 'FontName', 'Gill Sans MT');xlabel('Transdiagnostic hubs ({\itr})');ylabel('Epicenter overlap (%)'); axis square

hold on%plotting absolute correlations as direction of gradients is meaningless
% displacement so the text does not overlay the data points
%text(T_transhubs.r+dx, T_epi_overlap.Functional+dy, txty,'FontSize',10, 'FontName', 'Gill Sans MT');
hold on

scatter(T_transhubs.r,T_epi_overlap.Structural,100, [0.5451    0.4902    0.7294],'filled');
set(gca,'xlim', [-1 1],'ylim', [0 100], 'XTick', [-1, -0.5 0 0.5 1],'YTick',[0 50 100],'TickDir','none','color','none','FontSize',10, 'FontName', 'Gill Sans MT');xlabel('Transdiagnostic hubs ({\itr})');ylabel('Epicenter overlap (%)'); axis square
text(T_transhubs.r+dx, T_epi_overlap.Functional+dy, txty,'FontSize',10, 'FontName', 'Gill Sans MT');
text(T_transhubs.r+dx, T_epi_overlap.Structural+dy, txty,'FontSize',10, 'FontName', 'Gill Sans MT');

subplot(1,2,2)% gradients vs disorder maps

col = [162,181,205]/255;  
scatter(Grad1_dis_cor,Grad2_dis_cor,100, col, 'filled');
hold on%
text(Grad1_dis_cor+dx, Grad2_dis_cor+dy, txty,'FontSize',10, 'FontName', 'Gill Sans MT');
set(gca,'xlim', [-1 1],'ylim', [-1 1], 'XTick', [-1, -0.5 0 0.5 1],'YTick',[-1, -0.5 0 0.5 1],'TickDir','none','color','none','FontSize',10, 'FontName', 'Gill Sans MT');xlabel('G1 ({\itr})');ylabel('G2 ({\itr})'); axis square
hold off



