%% Longitudinal trajectories of resilient psychosocial functioning link to 
%% ongoing cortical myelination and functional reorganization during adolescence
% This script generates results presented in Figure 2 of the manuscript. It
% links longitudinal change in resilience scores to longitudinal
% trajectories of myelin-sensitive MT and contextualizes results with
% functional connectivity data.
% (c) Meike Hettwer 2024

% define and add paths
clc; clear; close all

% set up your paths
baseDir     = '';
dataDir     = '';
dataOutDir = '';
toolboxDir = '';
figDir      = '';

cd(dataOutDir)

% add required toolboxes to path
addpath(genpath(baseDir))
addpath(genpath(dataDir))
addpath(genpath(dataOutDir))
addpath(genpath([toolboxDir 'ENIGMA-1.1.3']))
addpath(genpath([toolboxDir 'fdr_bh']))
addpath(genpath([toolboxDir 'cbrewer/']))
addpath(genpath([toolboxDir 'Violinplot-Matlab-master/']))
addpath(genpath([toolboxDir 'RainCloudPlots-master']))

parpooling = 0;
% parallelize permutations?
if parpooling ==1
    parpool('Processes',20)
end

%% load data
load([dataOutDir 'subjects_to_use_for_delta.mat'])
load([dataOutDir 'included_parcels_fc.mat'])
load([dataOutDir '/NSPN_fc_mt_ct_subsample.mat'])%for 512 sessions from 295 individuals
load([dataOutDir 'Delta_10000_permutation_indices.mat'])
resilience_scores = readtable([dataOutDir '/resilience_based_on_pfactor/resilience_scores_prediction_repeatedmeasures_totalscores.csv']);

% atlases / parcellations
load('/yeo_glasser360.mat');
load('cortical_types_360.mat','types360','cmap_cortical_types')%
parc = readmatrix([toolboxDir 'ENIGMA-1.1.3/enigmatoolbox/datasets/parcellations/glasser_360_conte69.csv']);
% labels
labels_yeo = {'Visual','Somatomotor','Dorsal attention','Ventral attention','Limbic','Frontoparietal','Default mode'};
labels_types = {'Konicortex','Eulaminate3','Eulaminate2','Eulaminate1','Dysgranular','Agranular'};

% colormaps
load('yeo_colormap.mat'); yeo_cmap_surface = yeo_cmap; yeo_cmap(1,:)=[];%remove gray for 0 parcels/midline
cmap_resilience = [ones(1,3).*0.9; cbrewer('seq','BuPu',255)];cmap_resilience(cmap_resilience<0)=0;
color_positive_delta = cmap_resilience(240,:);
color_negative_delta = cmap_resilience(100,:);
RR = [6,117,111]/255;
RR3 = [150,150,150]/255;
RR4 = [0.8 0.8 0.8];

% set up variables
not_included_parcels = 1:360; not_included_parcels(included_parcels)=[];

subjects = subject_codes_to_use_for_delta; 
age_delta_mri = nan(numel(subjects),1);
mean_age = nan(numel(subjects),1);
sex_rep =cell(numel(subjects),1);
site_rep = cell(numel(subjects),1);
mt_delta_mean = nan(360,numel(subjects));
mt_delta_layers = nan(10,360,numel(subjects));
fc_delta = nan(330,330, numel(subjects));
subs_with_rep_mri_sessions = nan(numel(subjects),1);
sub_exist = zeros(length(subjects),1);
sub_repeated_behav = zeros(length(subjects),1);
mean_resilience_score = nan(length(subjects),1);
delta_resilience_score = nan(length(subjects),1);
baseline_resilience_score = nan(length(subjects),1);
age_delta_res = nan(length(subjects),1);
mt_baseline_mean = nan(360,length(subjects));

%% 1. compute deltas (i.e. simple T2 - T1)
% imaging
MTmean = squeeze(mean(MT,1));
FC = FC(17:end, 17:end,:);

for sub = 1:numel(subjects) %loop through 141 subjects
    idx = find(strcmp(subj,subjects(sub)));%find in 512 sessions      
    if length(idx)>1 
        subs_with_rep_mri_sessions(sub)=1; 
        age_delta_mri(sub) = age(idx(end))-age(idx(1));
        mean_age(sub) = mean([age(idx(end)) age(idx(1))]);
        sex_rep(sub) = sex(idx(1));
        site_rep(sub) = site(idx(1));

        mt_baseline_mean(:,sub) = MTmean(:,idx(1));

        % deltas
        mt_delta_layers(:,:,sub) = MT(:,:,idx(end))-MT(:,:,idx(1));
        mt_delta_mean(:,sub) = MTmean(:,idx(end))-MTmean(:,idx(1));
        fc_delta(:,:,sub) = FC(:,:,idx(end))-FC(:,:,idx(1));
    else
        subs_with_rep_mri_sessions(sub)=0;
    end
end

% resilience scores
for sub = 1:numel(subjects) % find 141 subjects in 1167 behavioral sessions
    if ~isnan(resilience_scores{str2double(subjects{sub}) == resilience_scores{:,1}, 7})
        scores = resilience_scores{str2double(subjects{sub}) == resilience_scores{:,1}, 7};
        age_res = resilience_scores{str2double(subjects{sub}) == resilience_scores{:,1}, 2};
        age_delta_res(sub) = age_res(end)-age_res(1);
        mean_resilience_score(sub) = mean(scores);
        baseline_resilience_score(sub) = scores(1);
        if length(scores)>1
            delta_resilience_score(sub) = (scores(end)-scores(1));
        end

        sub_exist(sub) = 1; %whether behavioral data is available for this subject
    else
        sub_exist(sub) = 0;
    end
end


%% 2. regress out time point difference
% a) for MT
clear slm
M1    = 1 + term(age_delta_mri);
slm  = SurfStatLinMod(mt_delta_mean',M1);
mt_delta_mean_agecorr = (mt_delta_mean' - (slm.X*slm.coef))';

%group_avg_mt_delta = mean(MTmean_agecorr,2);
%figure;plot_cortical(parcel_to_surface(group_avg_mt_delta,'glasser_360_conte69'),'surface_name','conte69','label_text','Delta Mean MT corrected','color_range',[min(group_avg_mt_delta),(min(group_avg_mt_delta))*-1])

% b) for MT layers
mt_delta_layers_agecorr = nan(size(mt_delta_layers,1),size(mt_delta_layers,2),size(mt_delta_layers,3));

for i = 1:size(mt_delta_layers,1) %10 layers

    this_layer = squeeze(mt_delta_layers(i,:,:))';
    clear slm
    M1    = 1 + term(age_delta_mri);
    slm  = SurfStatLinMod(this_layer,M1);
    mt_delta_layers_agecorr(i,:,:) = squeeze(mt_delta_layers(i,:,:)) - (slm.X*slm.coef)';
end

% c)for FC data
fc_delta_agecorr = nan(330,330,numel(subjects_to_use_for_delta));
for i = 1:size(fc_delta,2)%loop through rois
    this_roi = squeeze(fc_delta(:,i,:))';
    keep_row = ~any(isnan(this_roi), 1);
    this_roi = this_roi(:,keep_row); %excluded nans (here all in one row, the original diagonal)

    if any(isnan(this_roi))
        fprintf('spotted nans outside of diagonal!!')
    end
    clear slm
    M2    = 1 + term(age_delta_mri); 
    slm  = SurfStatLinMod(this_roi,M2);
    fc_delta_agecorr(i,keep_row,:) = ( this_roi - slm.X*slm.coef)' ;
end

% d) for resilience scores
clear slm
M3    = 1 + term(age_delta_res);
slm  = SurfStatLinMod(delta_resilience_score, M3);
delta_resilience_agecorr = delta_resilience_score - (slm.X*slm.coef);

%% winsorize deltas
% MT
mt_delta_winsorized = nan(size(mt_delta_mean_agecorr,1),size(mt_delta_mean_agecorr,2));
for roi = 1:360
    mt_delta_winsorized(roi,:) =  winsorize_by_sd(mt_delta_mean_agecorr(roi,:),3);
end

% MT layers
mt_delta_layers_winsorized = nan(size(mt_delta_layers_agecorr,1),size(mt_delta_layers_agecorr,2),size(mt_delta_layers_agecorr,3));
for lay = 1:size(mt_delta_layers_agecorr,1)
    for roi = 1:size(mt_delta_layers_agecorr,2)
        mt_delta_layers_winsorized(lay,roi,:) = winsorize_by_sd(squeeze(mt_delta_layers_agecorr(lay,roi,:)),3);
    end
end

%resilience delta
delta_resilience_score_winsorized = winsorize_by_sd(delta_resilience_agecorr,3);

% FC
delta_fc_winsorized = nan(size(fc_delta_agecorr,1),size(fc_delta_agecorr,2),size(fc_delta_agecorr,3));
for roi = 1:330
    for iroi = 1:330
        delta_fc_winsorized(roi,iroi,:) = winsorize_by_sd(squeeze(fc_delta_agecorr(roi,iroi,:)),3);
    end
end

% continue with winsorized data
delta_mt_mean = mt_delta_winsorized;
mt_delta_layers = mt_delta_layers_winsorized;
delta_resilience_score = delta_resilience_score_winsorized;
fc_delta = delta_fc_winsorized;

%% combine demographics
delta_table = [array2table(subjects,'VariableNames',{'subjects'}),...
    array2table(mean_age,'VariableNames',{'mean_age'}),...
    table(sex_rep,'VariableNames',{'sex'}),...
    array2table(site_rep,'VariableNames',{'site'}),...
    array2table(delta_resilience_score,'VariableNames',{'delta_resilience_score'}),...
    array2table(mean_resilience_score,'VariableNames',{'mean_resilience_score'})];

%% St up GLM res vs MRI
Age = delta_table.mean_age; Age = term(Age);
Sex = delta_table.sex; Sex = term(Sex);
Site = delta_table.site; Site = term(Site);
Mean_resilience = delta_table.mean_resilience_score; Mean_resilience = term(Mean_resilience);
Delta_resilience = delta_table.delta_resilience_score; Delta_resilience = term(Delta_resilience);

M_main_delta = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;

%% MT 
slm = SurfStatLinMod(delta_mt_mean', M_main_delta);
slm = SurfStatT( slm, delta_resilience_score);
mt_res_effect = slm.t;
figure;plot_cortical(parcel_to_surface(mt_res_effect,'glasser_360_conte69'),'surface_name','conte69','label_text','Delta Resilience * Delta Mean MT','color_range',[-5,5])

mt_res_effect_perm = nan(360,10000);
for perm = 1:10000
    permuted_delta_resilience = delta_table.delta_resilience_score(perm_idx(:,perm)); 
    permuted_mean_resilience = delta_table.mean_resilience_score(perm_idx(:,perm));
    Delta_resilience = term(permuted_delta_resilience);
    Mean_resilience = term(permuted_mean_resilience);

    M_perm = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;
    clear slm
    slm = SurfStatLinMod(delta_mt_mean', M_perm);
    slm = SurfStatT( slm, permuted_delta_resilience);
    mt_res_effect_perm(:,perm) = slm.t;

    fprintf(['Permutation ' num2str(perm) '\n'])
end

p_perm_fdr = get_perm_fdr_p(mt_res_effect,mt_res_effect_perm,360);

f = figure;meike_plot_cortical(parcel_to_surface(mt_res_effect.*p_perm_fdr','glasser_360_conte69'),'surface_name',...
    'conte69','cmap_custom',cmap_resilience,'label_text','Delta Resilience * Delta Mean MT','color_range',[0,4])
exportfigbo(f,[figDir, 'Delta_Resilience_Delta_MT_10000perm05.png'],'png', 12)
%unthresholded for supplementaries
f = figure;meike_plot_cortical(parcel_to_surface(mt_res_effect','glasser_360_conte69'),'surface_name',...
    'conte69','cmap_custom',cmap_resilience,'label_text','Delta Resilience * Delta Mean MT','color_range',[0,4])
exportfigbo(f,[figDir, 'Delta_Resilience_Delta_MT_unthresholded.png'],'png', 12)
%% contextualize MT results
%% Figure 2Aii) layer-specific effects within significant parcels
sig_mt_rois = find(p_perm_fdr);
layer_effect = nan(10,numel(sig_mt_rois)); layer_effectp = nan(10,numel(sig_mt_rois));
for lay=1:10
    for roi = 1:sum(p_perm_fdr)
        parcel_mean = squeeze(mt_delta_layers(lay,sig_mt_rois(roi),:));
        this_tbl = [array2table(parcel_mean,'VariableNames',{'mri'}) delta_table];
        lm = fitlm(this_tbl,'mri ~ delta_resilience_score + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
        layer_effect(lay,roi) =  lm.Coefficients.tStat(strcmp(lm.VariableNames,'delta_resilience_score'));
        layer_effectp(lay,roi) =  lm.Coefficients.pValue(strcmp(lm.VariableNames,'delta_resilience_score'));
    end
end

f_layers=figure('Position',[400 400 160 250]);
h = boxplot(layer_effect', 'Colors', RR3, 'Widths', 0.7,'Symbol','','orientation', 'horizontal');
set(h, 'LineWidth', 1.5);box off;xticks([0 2.5 5]), xlim([0 5]); yticks([])
exportfigbo(f_layers, [figDir 'layer_wise_mt_change_resilience.png'],'png', 12)

%% Figure 2B) Cortical types decoding
f_violin_decoding = figure('Position',[400 400 380 250]);
cyto_cmap_surface = cmap_cortical_types;
min_sig = min(mt_res_effect(sig_mt_rois));
plot([0 7],[min_sig min_sig],'Color','black');hold on
for class = 1:6
    mt_plot_this.(labels_types{class}) = mt_res_effect(types360==class)  ;
end
violinplot(mt_plot_this,labels_types,'ViolinColor',cmap_cortical_types(2:end,:),'MarkerSize',10, 'GroupOrder', labels_types); box off; ylim([-5 5]);
yticks([-5 0 5])
exportfigbo(f_violin_decoding, [figDir 'mt_res_effect_violin_cortical_types.png'],'png', 12);

%% 2C i) FC global connectivity change in PFC cluster
sig_mt_330 = p_perm_fdr(included_parcels); %reduce to 330 parcels
delta_mt_fc = squeeze(mean(fc_delta(logical(sig_mt_330),:,:),1,'omitnan')); 

fc_table = [table(mean(delta_mt_fc,1)', 'VariableNames',{'seed_fc'}) delta_table];
lm = fitlm(fc_table,'seed_fc ~ delta_resilience_score + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
t_seed_global =  lm.Coefficients.tStat(strcmp(lm.VariableNames,'delta_resilience_score'));
p_seed_global =  lm.Coefficients.pValue(strcmp(lm.VariableNames,'delta_resilience_score'));

%% 2C ii) FC connection from MT cluster to yeo networks
yeo_included = yeo_glasser(included_parcels);

% yeo_cmap_surface_masked = yeo_cmap_surface;
% yeo_cmap_surface_masked(9,:) = [0,0,0];
%% Figure 2C - seed based FC effect

clear slm
slm = SurfStatLinMod(delta_mt_fc', M_main_delta);
slm = SurfStatT( slm, delta_resilience_score);
seed_effect = slm.t;

%permutation
fc_res_effect_perm = nan(330,10000);
for perm = 1:10000
    permuted_delta_resilience = delta_table.delta_resilience_score(perm_idx(:,perm));
    permuted_mean_resilience = delta_table.mean_resilience_score(perm_idx(:,perm));
    Delta_resilience = term(permuted_delta_resilience);
    Mean_resilience = term(permuted_mean_resilience);

    M_perm = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;
    clear slm
    slm = SurfStatLinMod(delta_mt_fc', M_perm);
    slm = SurfStatT( slm, permuted_delta_resilience);
    fc_res_effect_perm(:,perm) = slm.t;
    fprintf(['Permutation ' num2str(perm) '\n'])
end

seed_fdr = get_perm_fdr_p(seed_effect,fc_res_effect_perm,330);

%prepare data to plot (project to 360)
to_plot = parcel_330_to_360(seed_effect.*seed_fdr','zer');
% add PFC mask
to_plot(logical(p_perm_fdr))=6;
%add mask for excluded parcels (due to low SNR)
to_plot(not_included_parcels) = -6;

cmap_seeds = [repmat([0.5 0.5 0.5],256,1); cmap_resilience;[0.2 0.2 0.2]/256];
f_seed_masked = figure;meike_plot_cortical(parcel_to_surface(to_plot,'glasser_360_conte69'),...
    'surface_name','conte69','color_range',[-4, 4],'cmap_custom', cmap_seeds,...
    'label_text','delta FC vs. delta resilience')
exportfigbo(f_seed_masked, [figDir 'seed_based_delta_fc_resilience_masked_against_10000perms.png'],'png', 12);

%% Figure 2D: contexturalize FC with yeo atlas
labels_yeo_adj = strrep(labels_yeo, ' ', '');

f_violin_decoding_yeo = figure('Position',[400 400 380 250]);

min_sig = min(seed_effect(logical(seed_fdr))); %add signiicance threshold to the plot
plot([0 7],[min_sig min_sig],'Color','black');hold on
for net = [1,2,3,4,6,7]%exclude limbic as limbic contains mostly missing parcels /low SNR
    fc_plot_this.(labels_yeo_adj{net}) = seed_effect(yeo_included==net);
end
violinplot(fc_plot_this,fliplr(labels_yeo_adj([1:4,6:7])),'GroupOrder',fliplr(labels_yeo_adj([1:4,6:7])),'ViolinColor',flipud(yeo_cmap([1:4,6:7],:)),'MarkerSize',10); box off; ylim([0 5]); %'MedianColor',[0 0 0],
yticks([0,2.5,5])
exportfigbo(f_violin_decoding_yeo, [figDir 'seed_based_delta_fc_yeo_violin.png'],'png', 12);
