%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coordinated Cortical Thickness alterations Across Mental Disorders: A
%Transdiagnostic ENIGMA Study
% Script 1: Loads ENIGMA summary statistics (Cohen's d maps), computes
% cross-disorder covariance hubs and epicenters

%This script runs analyses presented in Figure 1 of the manuscript.
%Meike Hettwer 2021 - CNG Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add Toolboxes 

%ENIGMA-1.1.3
%BrainSpace-0.1.2
%SurfStat
%cbrewer
%BCT
%GRETNA-2.0.0_release

%% load data from GitHub

%% load data
load('CT_6_disorders.mat'); % ENIGMA Cohen's d values for case-control differences in cortical thickness
load('corr_dis.mat'); %transdiagnostic covariance matrix

%load info which parcels/vertices are midbrain vertices to exclude during visualization/spin tests 
%(note: extra midbrain parcels when mapping from fsa5 to Desikan-Killiany are: 1 5 40)
load('midbrain_vertices_fsa5.mat');
load('DK_midbrain_parcels.mat');

%If data is loaded individually from the ENIGMA Toolbox
%{
%Schizophrenia
sum_stats = load_summary_stats('schizophrenia')
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
%}
%% Included disorders: BD, SZC, ADHD, ASD, MDD, OCD
disorders = [CT_d_bipo, CT_d_adhd, CT_d_asd, CT_d_mdd, CT_d_ocd, CT_d_sz];
disorder_names = {'Bipolar Disorder'; 'ADHD'; 'ASD'; 'Depression'; 'OCD'; 'Schizophrenia'};

%% Visualize Disorder-specific Cortical Thickness changes (*Figure 1A*)
range = [-0.35 0.35];
f1 = figure; plot_cortical(parcel_to_surface(CT_d_bipo), 'color_range', range, 'label_text', 'Bipolar Disorder Cohen´s d');
f2 = figure; plot_cortical(parcel_to_surface(CT_d_sz), 'color_range', range, 'label_text', 'Schizophrenia Cohen´s d');
f3 = figure; plot_cortical(parcel_to_surface(CT_d_adhd), 'color_range', range, 'label_text', 'ADHD Cohen´s d');
f4 = figure; plot_cortical(parcel_to_surface(CT_d_asd), 'color_range', range, 'label_text', 'ASD Cohen´s d');
f5 = figure; plot_cortical(parcel_to_surface(CT_d_ocd), 'color_range', range, 'label_text', 'OCD Cohen´s d');
f6 = figure; plot_cortical(parcel_to_surface(CT_d_mdd), 'color_range', range, 'label_text', 'Depression Cohen´s d');

%% Load cortico-cortical functional (rs-fMRI) and structural (DTI) connectivity data. 
%Data is accessed through the ENIGMA Toolbox and derived from 207 HCP young adults
[fc_ctx, fc_ctx_labels, ~, ~] = load_fc();
[sc_ctx, sc_ctx_labels, ~, ~] = load_sc();

% plot HCP connectivity matrices (*Figure 1B*)
reds = cbrewer('div','RdBu',100); reds = flipud(reds(1:(100/2),:)); reds(reds<0)=0;
figure; imagesc(fc_ctx, [0 1]); %plot functional connectivity matrix
axis square; colormap(reds);colorbar; set(gca, 'YTick', [], 'XTick', [])

purple = cbrewer('seq','BuPu',100); purple(purple<0) = 0;
figure; imagesc(sc_ctx, [0 10]); %plot structural connectivity matrix
axis square; colormap(purple);colorbar; set(gca, 'YTick', [], 'XTick', [])
 
%% Hubs (*Figure 1C*)
% Compute normative connectivity degree centrality (hubs) for top 20%
% connections, by taking the sum of strong connections per parcel
fc_ctx_bin_thresh = bsxfun(@gt, fc_ctx, prctile(fc_ctx,80))'; 
sc_ctx_bin_thresh = bsxfun(@gt, sc_ctx, prctile(sc_ctx,80))';  
fc_ctx_dc_thresh = sum(fc_ctx_bin_thresh);
sc_ctx_dc_thresh = sum(sc_ctx_bin_thresh);

plot_cortical(parcel_to_surface(fc_ctx_dc_thresh))
% Set up thresholded disorder correlation matrix and dc map
corr_dis_bin_thresh = bsxfun(@gt, corr_dis, prctile(corr_dis,80))'; 
figure; imagesc(corr_dis_bin_thresh);set(gca,'XTick',[],'YTick',[]);colormap(cbrewer('seq','Reds',50));colorbar; axis square;

% plot degree centrality map / cross-disorder covariance hub map 
dis_dc = sum(corr_dis_bin_thresh);
figure; plot_cortical(parcel_to_surface(dis_dc, 'aparc_fsa5'),'color_range',[0 40])

% Compute correlation between normative connectome hubs and transdiagnostic covariance hubs
% functional hubs
[r1,p1] = corr(dis_dc',fc_ctx_dc_thresh');
[p1_spin, r1_dist] = spin_test(dis_dc',fc_ctx_dc_thresh', 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
% structural hubs
[r2,p2] = corr(dis_dc',sc_ctx_dc_thresh');
[p2_spin, r2_dist] = spin_test(dis_dc',sc_ctx_dc_thresh', 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
                               
%% Epicenters (*Figure 1D*)
% Identify cortical epicenter values (from functional connectivity) by
% systematically correlating seed-based connectivity profiles with
% covariance hub map
fc_ctx_trans              = zeros(size(fc_ctx, 1), 1);
fc_ctx_trans_p            = zeros(size(fc_ctx, 1), 1);
for seed = 1:size(fc_ctx, 1)
    seed_conn           = fc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn, dis_dc);
    fc_ctx_trans(seed)    = r_tmp(1, 2);
    fc_ctx_trans_p(seed)  = spin_test(seed_conn, dis_dc, 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Identify cortical epicenter values (from structural connectivity)
sc_ctx_trans              = zeros(size(sc_ctx, 1), 1);
sc_ctx_trans_p            = zeros(size(sc_ctx, 1), 1);
for seed = 1:size(sc_ctx, 1)
    seed_conn           = sc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn, dis_dc);
    sc_ctx_trans(seed)    = r_tmp(1, 2);
    sc_ctx_trans_p(seed)  = spin_test(seed_conn, dis_dc, 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Select only regions with p < 0.05 (functional epicenters)
fc_ctx_trans_p_sig = zeros(length(fc_ctx_trans_p), 1);
fc_ctx_trans_p_sig(fc_ctx_trans_p < 0.05) = fc_ctx_trans(fc_ctx_trans_p<0.05);
figure; plot_cortical(parcel_to_surface(fc_ctx_trans_p_sig, 'aparc_fsa5'), ...
                'color_range', [0 0.5]');%Red_colorbar

% Selecting only regions with p < 0.05 (structural epicenters)
sc_ctx_trans_p_sig = zeros(length(sc_ctx_trans_p), 1);
sc_ctx_trans_p_sig(sc_ctx_trans_p < 0.05) = sc_ctx_trans(sc_ctx_trans_p<0.05);
figure; plot_cortical(parcel_to_surface(sc_ctx_trans_p_sig, 'aparc_fsa5'), ...
                'color_range', [0 0.5], 'cmap', 'RdBu_r')%cbrewer('seq','BuPu',100)

%% Compute subcortical epicenters in the same manner
% Load subcortico-cortical functional connectivity data
[~, ~, fc_sctx, fc_sctx_labels] = load_fc();

% Load subcortico-cortical structural connectivity data
[~, ~, sc_sctx, sc_sctx_labels] = load_sc();

% Epicenter mapping - subcortical seeds but cortical map!
% Identify subcortical epicenters (from functional connectivity)
fc_sctx_epi             = zeros(size(fc_sctx, 1), 1);
fc_sctx_epi_p           = zeros(size(fc_sctx, 1), 1);
for seed = 1:size(fc_sctx, 1)
    seed_conn           = fc_sctx(seed, :);
    r_tmp               = corrcoef(seed_conn, dis_dc);
    fc_sctx_epi(seed)   = r_tmp(1, 2);
    fc_sctx_epi_p(seed) = spin_test(seed_conn, dis_dc, 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Identify subcortical epicenters (from structural connectivity)
sc_sctx_epi             = zeros(size(sc_sctx, 1), 1);
sc_sctx_epi_p           = zeros(size(sc_sctx, 1), 1);
for seed = 1:size(sc_sctx, 1)
    seed_conn           = sc_sctx(seed, :);
    r_tmp               = corrcoef(seed_conn, dis_dc);
    sc_sctx_epi(seed)   = r_tmp(1, 2);
    sc_sctx_epi_p(seed) = spin_test(seed_conn, dis_dc, 'surface_name', 'fsa5', 'parcellation_name', ...
                                    'aparc', 'n_rot', 1000, 'type', 'pearson');
end

% Threshold and Project the results on the surface brain
% functional epicenters
fc_sctx_epi_p_sig = zeros(length(fc_sctx_epi_p), 1);
fc_sctx_epi_p_sig(fc_sctx_epi_p < 0.05) = fc_sctx_epi(fc_sctx_epi_p<0.05);
figure; plot_subcortical(fc_sctx_epi_p_sig, 'ventricles', 'False', ...
                    'color_range', [0 0.5], 'cmap', 'RdBu_r')

% structural epicenters
sc_sctx_epi_p_sig = zeros(length(sc_sctx_epi_p), 1);
sc_sctx_epi_p_sig(sc_sctx_epi_p < 0.05) = sc_sctx_epi(sc_sctx_epi_p<0.05);
figure; plot_subcortical(sc_sctx_epi_p_sig, 'ventricles', 'False', ...
                    'color_range', [0 0.5], 'cmap', 'RdBu_r')

               
%% END OF ANALYSES VISUALIZED IN FIGURE 1                
               
