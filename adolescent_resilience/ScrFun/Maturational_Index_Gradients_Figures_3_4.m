%% Maturational index and Gradients
% This script computes MPC and FC maturational indices, compares groups of
% individuals showing an increase vs. decrease in resilience capacities,
% and compares this group difference to nullmodels of random/permuted groups.
% It also computes the principal axis of MPC age effects and plots group
% differences

%% define your paths
baseDir     = '';
dataDir     = '';
dataOutDir = '';
toolboxDir = '';
figDir      = '';
cd(dataOutDir)

%% add toolboxes
addpath(genpath(baseDir))
addpath(genpath(dataDir))
addpath(genpath(dataOutDir))
addpath(genpath([toolboxDir 'micaopen-master/MPC/']))
addpath(genpath([toolboxDir 'surfstat/']))
addpath(genpath([toolboxDir 'ENIGMA-1.1.3']))
addpath(genpath([toolboxDir 'BrainSpace-0.1.2']))
addpath(genpath([toolboxDir 'fdr_bh']))
addpath(genpath([toolboxDir 'cbrewer/']))
addpath(genpath([toolboxDir 'micaopen-master/surfstat/surfstat_addons/']))
addpath(genpath([toolboxDir 'compare_correlation_coefficients/']))

%% load data
load([dataOutDir '/NSPN_fc_mt_subsample.mat']);%for 512 sessions
load([dataOutDir 'subjects_to_use_for_delta.mat']);
load([dataOutDir 'delta_table_extracted.mat']);
load([dataOutDir 'included_parcels_fc.mat']); not_included_parcels = 1:360; not_included_parcels(included_parcels)=[];

%% colormaps, labels, atlases
load([toolboxDir  'ScientificColourMaps7/roma/roma.mat'])
roma = flipud(roma);
cmap_resilience = [ones(1,3).*0.9; cbrewer('seq','BuPu',255)];
color_positive_delta = cmap_resilience(250,:);
color_negative_delta = cmap_resilience(70,:);
cmap_MI_diff = flipud(cbrewer('div','RdBu',256));

RR = [6,117,111]/255;
RR2 = [2,189,184]/255;
RR3 = [150,150,150]/255;
RR4 = [0.8 0.8 0.8];
RR_colors = [RR; RR2; RR3; RR4];
category_cmap = [1 1 1; [106 29 77;145, 20, 64; 149, 176, 182 ;46, 96, 107]/256];
load('yeo_colormap.mat'); yeo_cmap_surface = yeo_cmap; yeo_cmap(1,:)=[];%includes gray for 0 parcels/midline, not needed here
mode_cmap = [RR3;cmap_resilience(10,:); cmap_resilience(30,:); cmap_resilience(90,:); cmap_resilience(70,:)];

%labels and atlases
load('yeo_glasser360.mat');
labels_yeo = {'Visual','Somatomotor','Dorsal attention','Ventral attention','Limbic','Frontoparietal','Default mode'};
yeo_included = yeo_glasser(included_parcels);
load('/data/p_02792/Resources/yeo_glasser360.mat');
labels_yeo_adj = strrep(labels_yeo, ' ', '');
load('cortical_types_360.mat','types360','cmap_cortical_types');
cmap_curv = [cmap_MI_diff(1:127,:);(ones(2,3).*0.9); cmap_MI_diff(130:end,:)];

%% which figures to plot and which analyses to run

plot_f_demographics_perm = 1;
plot_f_demographics = 1;
run_group_MIs = 1;
run_permutation_indices = 1;
run_permutations = 1;
recompute_normative_MI = 1;
compute_overlap = 1;
parpooling =1;
compute_gradients = 1;

if parpooling == 1
    parpool('Processes',30)
end
%% start analyses

MPC = nan(size(MT,2),size(MT,2),size(MT,3));
for sub = 1:length(subj)
    MPC(:,:,sub) = build_mpc(MT(:,:,sub), []); %uses code from MICA githup repository micaopen-master/MPC/
end
FC = FC(17:end,17:end,:); %exclude subcortical regions (1:16)

%% check true demographics to match permutation set-up for later permutation (balanced by sex and age strata)

subjects = unique(subj);%295
baseline_age = nan(numel(subjects),1); %baseline_sex = nan(numel(subjects),1);
for sub = 1:numel(subjects)%295 loops through subjects,note measurement timepoint!
    ages = age(strcmp(subj,subjects(sub)));%this goes into 512 subj variable to match age
    baseline_age(strcmp(subjects,subjects(sub))) = ages(1);
    sexes = sex(strcmp(subj,subjects(sub)));%this goes into 512 subj variable to match age
    baseline_sex(strcmp(subjects,subjects(sub))) = sexes(1);
end
baseline_age = baseline_age(subjects_to_use_for_delta);
baseline_sex = baseline_sex(subjects_to_use_for_delta);

% histograms for sex and age distributions
if plot_f_demographics == 1
    f_demographics = figure('Position',[400 400 480 400]);
    subplot(1,2,1)
    % Define custom age bin edges
    ageBinEdges = [14, 16, 18, 20, 22, 24.99];
    % Extract age data
    agePositiveDelta = baseline_age(delta_table_extracted.delta_resilience > 0);
    ageNegativeDelta = baseline_age(delta_table_extracted.delta_resilience < 0);
    % Plot histograms
    histogram(agePositiveDelta, ageBinEdges, 'FaceAlpha', 0.5,'FaceColor',RR2,'FaceAlpha',0.5,'LineWidth',1); hold on
    histogram(ageNegativeDelta, ageBinEdges, 'FaceAlpha', 0.5,'FaceColor',RR3,'FaceAlpha',0.5,'LineWidth',1);
    box off
    set(gca,'FontName','Seaford');yticks([0 10 20]);ylim([0 20]);xticklabels({'14-15','16-17','18-19','20-21','>22'});xticks(15:2:23)
    %legend({'+delta', '-delta'});
    title('Age');set(gca,'FontName','Seaford','FontSize',10);

    subplot(1,2,2)
    histogram(strcmp(delta_table_extracted.sex(delta_table_extracted.delta_resilience > 0),'female'),'FaceColor',RR2,'FaceAlpha',0.5,'LineWidth',1);hold on;
    histogram(strcmp(delta_table_extracted.sex(delta_table_extracted.delta_resilience < 0),'female'),'FaceColor',RR3,'FaceAlpha',0.5,'LineWidth',1);hold on;
    legend({' + \Delta SRS',' - \Delta SRS'}); title('Sex');
    xticks([0 1]);xticklabels({'Male', 'Female'}); legend boxoff
    box off;set(gca,'FontName','Seaford','FontSize',10);yticks(0:10:50);ylim([0 50])
    exportfigbo(f_demographics,[figDir 'demographics_true_delta_groups.png'],'png',10)
end
% how many females in each group
percent_female_resilient = sum(strcmp(delta_table_extracted.sex(delta_table_extracted.delta_resilience > 0),'female')) / sum(delta_table_extracted.delta_resilience > 0);
percent_female_vulnerable = sum(strcmp(delta_table_extracted.sex(delta_table_extracted.delta_resilience < 0),'female')) / sum(delta_table_extracted.delta_resilience < 0);

%% set up demographics tables
resilience_scores = readtable([dataOutDir '/resilience_based_on_pfactor/resilience_scores_prediction_repeatedmeasures_totalscores.csv']);
% extract group membership for each subject. As we study longitudinal
% change, 'resilient' and 'vulnerable' refer to 'becoming more
% resilient/vulnerable with age'

resilience_group = nan(numel(subj),1);
for sub = 1:numel(subj) %here we need all sessions for the MI computation
    sub_idx = find(strcmp(subj,subj(sub)));
    if ~isnan(delta_table.delta_resilience(strcmp(delta_table.subjects,subj(sub))))
        resilience_group(sub_idx,:) = delta_table.delta_resilience(strcmp(delta_table.subjects,subj(sub))) > 0;
    end
end

% create demographic tables for each roup
demo_tbl = [table(subj) table(age) table(sex) table(site','VariableNames',{'site'}) table(resilience_group,'VariableNames',{'resilience_group'})];
demo_tbl_resilient = demo_tbl(demo_tbl.resilience_group==1,:);
demo_tbl_vulnerable = demo_tbl(demo_tbl.resilience_group==0,:);
% extract MPC matrices per group
MPC_resilience = MPC(:,:,demo_tbl.resilience_group==1); %193/512 sessions
MPC_vulnerable = MPC(:,:,demo_tbl.resilience_group==0); %153/512 session
% true indices
resilient_indices = find(demo_tbl.resilience_group==1);
vulnerable_indices = find(demo_tbl.resilience_group==0);
% save
save([dataOutDir 'delta_resilience_group_indices.mat'],'resilient_indices','vulnerable_indices','MPC_resilience','MPC_vulnerable')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normative MI all subjects %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. maturational index - normative development
% the general / full sample MI is based on all individuals that have at
% least 2 measurement timepoints.
% skip if already computed
all_demo_tbl = [table(age,'VariableNames', {'age'}), table(sex,'VariableNames',{'sex'}), table(site','VariableNames',{'site'}),...
    table(subj,'VariableNames',{'subj'})];

%  filter out only individuals with at least two measurement timepoints
%  (n=199)
counts = tabulate(all_demo_tbl.subj);
counts_numeric = cellfun(@double, counts(:, 2));
subjects_to_keep = ismember(all_demo_tbl.subj, counts(counts_numeric > 1, 1));
all_demo_tbl_repeated = all_demo_tbl(subjects_to_keep,:);

MPC_repeated = MPC(:,:,subjects_to_keep);
FC_repeated = FC(:,:,subjects_to_keep);

%compute MI for MPC and FC (Figure 3 B, D, F)
if recompute_normative_MI == 1
    [MI_MPC, p_MPC_MI, MI_MPC_ratio] = compute_MI_lme(MPC_repeated, all_demo_tbl_repeated);
    p_MPC_fdr = fdr_bh(p_MPC_MI', 0.05);
    f1=figure; plot_cortical(parcel_to_surface(MI_MPC.*(p_MPC_fdr'),'glasser_360_conte69'),'surface_name','conte69','color_range',[-1 1],'label_text','MI MPC (whole sample; FDR < 0.05)');
    save([dataOutDir 'MI_MPC.mat'],'MI_MPC','p_MPC_MI','p_MPC_fdr','MI_MPC_ratio');
    exportfigbo(f1, [figDir 'MI_MPC_normative.png'],'png',10)
    [MI_FC, p_FC_MI, MI_FC_ratio] = compute_MI_lme(FC_repeated, all_demo_tbl_repeated);
    p_FC_fdr = fdr_bh(p_FC_MI', 0.05)';
    save([dataOutDir 'MI_FC.mat'],'MI_FC','p_FC_MI','p_FC_fdr','MI_FC_ratio');

    %mask missing parcels for FC (those with low SNR)
    mask_fc = parcel_330_to_360(MI_FC.*(p_FC_fdr),'zer');
    mask_fc(not_included_parcels) = -6;
    cmap_mask_fc = [ 0.6 0.6 0.6; flipud(cbrewer('div','RdBu',98)); 0.6 0.6 0.6]; cmap_mask_fc(cmap_mask_fc<0)=0;

    f2=figure; plot_cortical(parcel_to_surface(mask_fc,'glasser_360_conte69'),'surface_name','conte69','color_range',[-1 1],'label_text','MI FC (whole sample; FDR < 0.05)');
    exportfigbo(f2, [figDir 'MI_FC_normative.png'],'png',10)
end

% if skipped, load pre-computed MIs
load([dataOutDir 'MI_MPC.mat']);
load([dataOutDir 'MI_FC.mat']);

%% visualize overlap between FC and MPC MIs, each thresholded at p<0.05 FDR
if compute_overlap == 1
    p_MPC_extracted_fdr = fdr_bh(p_MPC_MI(included_parcels),0.05);

    MI_MPC_extracted_fdr = MI_MPC(included_parcels).*p_MPC_extracted_fdr;
    MI_FC_fdr = MI_FC.*p_FC_fdr;

    ov = zeros(length(MI_MPC_extracted_fdr),1);
    pos_pos = ov; pos_pos(MI_MPC_extracted_fdr>0 & MI_FC_fdr>0)=1;
    neg_neg = ov; neg_neg(MI_MPC_extracted_fdr<0 & MI_FC_fdr<0)=1;
    pos_neg = ov; pos_neg(MI_MPC_extracted_fdr>0 & MI_FC_fdr<0)=1;
    neg_pos = ov; neg_pos(MI_MPC_extracted_fdr<0 & MI_FC_fdr>0)=1;

    % define overlap in cross-modal maturational modes, pos = conservative,
    % neg = disruptive
    directions = ov;
    directions(logical(pos_pos))=1; directions(logical(pos_neg))=2; directions(logical(neg_pos))=3; directions(logical(neg_neg))=4;
    directions_360 = parcel_330_to_360(directions,'zer');

    plot_directions = directions_360;
    plot_directions(not_included_parcels) = -6; %grey mask

    category_cmap = [1 1 1;  [148 17 6; 239 138 110 ;85 155 179 ; 33 80 101]/256];
    cmap_mask_fc = [ 0.6 0.6 0.6; category_cmap];

    f3=figure; plot_cortical(parcel_to_surface(plot_directions,'glasser_360_conte69'),'surface_name','conte69', 'label_text','maturational categories','color_range',[-1 4]);
    exportfigbo(f3, [figDir 'MI_maturational_categories.png'],'png',10)
    close all

    % visualize proportion of overlap modes in each yeo network in a bar
    % graph (Figure 3G)
    dir_yeo_prob = zeros(5,7);
    for net = 1:7 
        idx = unique(directions_360(yeo_glasser == net)) +1; %indexing into networks)
        counts = histcounts(directions_360(yeo_glasser == net),'Normalization', 'probability' );
        dir_yeo_prob(idx,net) = counts(idx)*100;
    end
    
    %exclude limbic, as many parcels here are missing due to low SNR
    dir_yeo_prob(:,5)=[];

    f_categories = figure('position',[400 400 450 350]);
    bh = bar(dir_yeo_prob','stacked','FaceColor','flat','LineWidth',0.75);
    for k = 1:5
        bh(k).CData = category_cmap(k,:);
    end
    xticks(1:6);xticklabels(labels_yeo([1:4 6:7]));yticks([0 50 100]);ylim([0 100]);box off;
    ax = gca;
    ax.XRuler.TickLength = [0 0];ax.YRuler.TickLength = [0 0];

    exportfigbo(f_categories, [figDir 'MI_Categories_Yeo.png'],'png',10)
    save([dataOutDir 'MI_categories.mat'],'directions','directions_360')
end
load([dataOutDir 'MI_categories.mat'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%
%% true group differences %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%

if run_group_MIs ==1

    %% MPC
    %resilient group
    this_age = age(resilient_indices); Age         = term(this_age);
    this_sex = sex(resilient_indices); Sex         = term(this_sex);
    this_subj = subj(resilient_indices); Subj        = term(var2fac(str2double(this_subj)));
    this_site = site(resilient_indices)';Site        = term(this_site);
    M_resilient          = 1 + Age + Sex + Site + random(Subj) + I;

    clear slm_mpc slm_mpc_model MI_MPC_resilient p_MI_MPC_resilient  MI_MPC_vulnerable p_MI_MPC_vulnerable_perm

    parfor roi = 1:size(MPC_resilience,1) 
        this_roi_mpc = squeeze(MPC_resilience(roi,:,:))';
        this_roi_mpc(:,roi)=[]; %remove nans to fit GLM
        this_roi_mpc(:,sum(this_roi_mpc)==0) = []; %check if any columns are all zero and remove to fit GLM (should be none)

        slm_mpc_model(roi) = SurfStatLinMod(this_roi_mpc, M_resilient);
        slm_mpc(roi) = SurfStatT( slm_mpc_model(roi), age(resilient_indices));
    end
    for roi = 1:size(MPC_resilience,1)
        MPC_baseline_resilient(roi,1:length(slm_mpc(roi).coef(1,:))) = (slm_mpc(roi).coef(1,:) + slm_mpc(roi).coef(2,:)*14  + (slm_mpc(roi).coef(3,:)*1/2) + (slm_mpc(roi).coef(4,:)*1/2) +(slm_mpc(roi).coef(5,:)*1/3)+ (slm_mpc(roi).coef(6,:)*1/3) + (slm_mpc(roi).coef(7,:)*1/3))';
        MPC_change_resilient(roi,1:length(slm_mpc(roi).coef(1,:))) = slm_mpc(roi).coef(2,:)';

        %compute MI
        [MI_MPC_resilient(roi), p_MI_MPC_resilient(roi)] = corr(MPC_baseline_resilient(roi,1:length(slm_mpc(roi).coef(1,:)))',MPC_change_resilient(roi,1:length(slm_mpc(roi).coef(1,:)))','type','Spearman','rows','complete');
    end

    %vulnerable group
    this_age = age(vulnerable_indices); Age         = term(this_age);
    this_sex = sex(vulnerable_indices); Sex         = term(this_sex);
    this_subj = subj(vulnerable_indices); Subj        = term(var2fac(str2double(this_subj)));
    this_site = site(vulnerable_indices)';Site        = term(this_site);
    M_vulnerable          = 1 + Age + Sex + Site + random(Subj) + I;

    clear slm_mpc slm_mpc_model
    parfor roi = 1:size(MPC_vulnerable,1) %takes approx 30 mins; with parfor 20 workers: 2-3 min

        this_roi_mpc = squeeze(MPC_vulnerable(roi,:,:))';
        this_roi_mpc(:,roi)=[]; %remove nans to fit GLM
        this_roi_mpc(:,sum(this_roi_mpc)==0) = []; %check if any columns are all zero and remove to fit GLM

        slm_mpc_model(roi) = SurfStatLinMod(this_roi_mpc, M_vulnerable);
        slm_mpc(roi) = SurfStatT( slm_mpc_model(roi), age(vulnerable_indices));
    end
    for roi = 1:size(MPC_vulnerable,1)
        MPC_baseline_vulnerable(roi,1:length(slm_mpc(roi).coef(1,:))) = (slm_mpc(roi).coef(1,:) + slm_mpc(roi).coef(2,:)*14  + (slm_mpc(roi).coef(3,:)*1/2) + (slm_mpc(roi).coef(4,:)*1/2) +(slm_mpc(roi).coef(5,:)*1/3)+ (slm_mpc(roi).coef(6,:)*1/3) + (slm_mpc(roi).coef(7,:)*1/3))';
        MPC_change_vulnerable(roi,1:length(slm_mpc(roi).coef(1,:))) = slm_mpc(roi).coef(2,:)';

        %compute MI
        [MI_MPC_vulnerable(roi), p_MI_MPC_vulnerable(roi)] = corr(MPC_baseline_vulnerable(roi,1:length(slm_mpc(roi).coef(1,:)))',MPC_change_vulnerable(roi,1:length(slm_mpc(roi).coef(1,:)))','type','Spearman','rows','complete');
    end
    save([dataOutDir 'MI_delta_groups_MPC_res_groups.mat'],'MI_MPC_vulnerable','p_MI_MPC_vulnerable','MI_MPC_resilient','p_MI_MPC_resilient','MPC_baseline_resilient','MPC_change_resilient','MPC_baseline_vulnerable','MPC_change_vulnerable');
end

% if skipped, load precomputed MIs
load([dataOutDir 'MI_delta_groups_MPC_res_groups.mat']);

%% permutation of random groups
% permutation set up:
% repeated sessions of same subject are allocated together
% permutation is done within age bins
% permutation keeps similar sex distribution as observed in true groups
% group size imbalances are considered within age bins

% Set the random seed for reproducibility
rng('default');
rng(42);

% Define age bins and the number of subjects in each age bin
baseline_age = baseline_table.baseline_age(subjects_to_use_for_delta);
ageBins = [14, 16, 18, 20, 22, 27];
numBins = length(ageBins) - 1;

age_bin_counts = zeros(1, numBins);
delta_table_extracted_female = delta_table_extracted(strcmp(delta_table_extracted.sex,'female'),:);
delta_table_extracted_male = delta_table_extracted(strcmp(delta_table_extracted.sex,'male'),:);

% Count the true number of subjects in each age bin
for i = 1:numBins
    lowerAge = ageBins(i);
    upperAge = ageBins(i+1);
    binSubjects = sum(lowerAge <= baseline_age & baseline_age < upperAge);
    these_subs = find(lowerAge <= baseline_age & baseline_age < upperAge);
    age_bin_counts(i) = binSubjects;
    % Count the number of subjects in each age bin within groups
    binSubjectsResilientAge = sum(lowerAge <= baseline_age(delta_table_extracted.delta_resilience>0) & baseline_age(delta_table_extracted.delta_resilience>0) < upperAge);
    binSubjectsVulnerableAge = sum(lowerAge <= baseline_age(delta_table_extracted.delta_resilience<0) & baseline_age(delta_table_extracted.delta_resilience<0) < upperAge);
    age_bin_counts_resvul(i,:) = [binSubjectsResilientAge binSubjectsVulnerableAge];

    binSubjectsResilientFemale = sum(strcmp(delta_table_extracted.sex(...
        lowerAge <= baseline_age & baseline_age <= upperAge & delta_table_extracted.delta_resilience > 0), 'female'));

    binSubjectsVulnerableFemale = sum(strcmp(delta_table_extracted.sex(...
        lowerAge <= baseline_age & baseline_age <= upperAge & delta_table_extracted.delta_resilience < 0), 'female'));

    sex_age_bin_fem_counts_resvul(i,:) = [binSubjectsResilientFemale' binSubjectsVulnerableFemale'];
end

% Calculate the number of subjects per group and per age bin
num_resilient = sum(delta_table_extracted.delta_resilience>0);
num_vulnerable = sum(delta_table_extracted.delta_resilience<0);
num_total = num_resilient + num_vulnerable;
prop_resilient = age_bin_counts_resvul(:,1) ./ sum(age_bin_counts_resvul,2);
prop_vulnerable = age_bin_counts_resvul(:,2) ./ sum(age_bin_counts_resvul,2);

% Calculate the number of females in each group per age bin
prop_female_resilient = sex_age_bin_fem_counts_resvul(:,1)./age_bin_counts_resvul(:,1); 
prop_female_vulnerable = sex_age_bin_fem_counts_resvul(:,2)./age_bin_counts_resvul(:,2);

if run_permutation_indices == 1
    resilient_permuted = [];
    vulnerable_permuted = [];
    %loop through permutations
    for perm = 1:10000
        perm;
        uniquePermutation = false;
        while ~uniquePermutation
            % Initialize resilient and vulnerable bins for this permutation
            resilient_bin = [];
            vulnerable_bin = [];

            for bin = 1:numBins
                lowerAge = ageBins(bin);
                upperAge = ageBins(bin+1);
                binSubjects = find( lowerAge <= baseline_age & baseline_age < upperAge); %find subjects of that age bin across groups
                resilientGroupSize = round(numel(binSubjects) *prop_resilient(bin));
                %sex
                these_sexes = delta_table_extracted.sex(binSubjects);
                find_females = find(strcmp(these_sexes,'female')); %indexing out of selected subjects
                find_males = find(strcmp(these_sexes,'male'));
                %how many females per group?
                femaleToRes = round( resilientGroupSize*prop_female_resilient(bin)); %out of all sunjects from that bin, how many go to the resilient group, and out of thise, how many should be female?

                %randomize female indices and allocate to groups
                fem_rand_idx = randperm(numel(find_females));
                fem_res = binSubjects(find_females(fem_rand_idx(1:femaleToRes))); %how many females go into the resilient group 
                fem_vul = binSubjects(find_females(fem_rand_idx(femaleToRes+1:end)));

                %randomize males and fill groups until group size is reached
                male_rand_idx = randperm(numel(find_males));
                malesToRes = resilientGroupSize-femaleToRes;
                malesToVul = numel(binSubjects) - resilientGroupSize - numel(fem_vul);

                male_res = binSubjects(find_males(male_rand_idx(1:malesToRes)));
                male_vul = binSubjects(find_males(male_rand_idx(malesToRes+1:malesToRes+malesToVul)));

                resilient_bin = [resilient_bin; fem_res; male_res];
                vulnerable_bin = [vulnerable_bin; fem_vul; male_vul];

            end

            % Convert combined_permutation to cell arrays for comparison
            bin_resilient_cell = num2cell(resilient_bin);
            bin_vulnerable_cell = num2cell(vulnerable_bin);
            resilient_cell = num2cell(resilient_permuted, 1);
            vulnerable_cell = num2cell(vulnerable_permuted, 1);

            % Check if the combined permutation is unique
            if ~any(cellfun(@(x) isequal(x, bin_resilient_cell), resilient_cell)) && ...
                    ~any(cellfun(@(x) isequal(x, bin_vulnerable_cell), vulnerable_cell))
                % Store the unique permutation
                resilient_permuted(:, perm) = resilient_bin;
                vulnerable_permuted(:, perm) = vulnerable_bin;
                uniquePermutation = true;
            end
        end
    end
    save([dataOutDir 'delta_permuted_group_indices_age_sex_balanced.mat'],'resilient_permuted','vulnerable_permuted')
end
% if skipped, load pre-computed
load([dataOutDir 'delta_permuted_group_indices_age_sex_balanced.mat'])

%plot histograms of random permutations
rand_idx = randi(10000,3,1);

if plot_f_demographics_perm == 1
    f_demographics_perm = figure('Position',[400 400 1000 700]);
    for i = 1:3
        subplot(2,3,i)
        histogram(baseline_age(resilient_permuted(:,rand_idx(i))), ageBinEdges, 'FaceAlpha', 0.5,'FaceColor',RR2,'FaceAlpha',0.5,'LineWidth',1); hold on
        histogram(baseline_age(vulnerable_permuted(:,rand_idx(i))), ageBinEdges, 'FaceAlpha', 0.5,'FaceColor',RR3,'FaceAlpha',0.5,'LineWidth',1);
        box off
        set(gca,'FontName','Seaford');yticks([0 10 20]);ylim([0 20]);xticklabels({'14-15','16-17','18-19','20-21','>22'});xticks(15:2:23)
        title(['Permutation ',num2str(rand_idx(i)) ]);set(gca,'FontName','Seaford','FontSize',10);

        subplot(2,3,i+3)
        histogram(strcmp(baseline_sex(resilient_permuted(:,rand_idx(i))),'female'),'FaceColor',RR2,'FaceAlpha',0.5,'LineWidth',1);hold on;
        histogram(strcmp(baseline_sex(vulnerable_permuted(:,rand_idx(i))),'female'),'FaceColor',RR3,'FaceAlpha',0.5,'LineWidth',1);hold on;
        %legend({' + \Delta SRS',' - \Delta SRS'});
        xticks([0 1]);xticklabels({'Male', 'Female'}); %legend boxoff
        box off;set(gca,'FontName','Seaford','FontSize',10);yticks(0:10:50);ylim([0 50])
    end
    exportfigbo(f_demographics_perm,[figDir 'demographics_permuted_delta_groups.png'],'png',10)
end

% extract session indices for each group based on participant allocation
% resilient permutation groups, as each participant has multiple
% measurement timepoints which were permuted together
for perm = 1:10000
    perm_indices_resilient_sessions = [] ;
    for sub = 1:size(resilient_permuted,1)
        this_sub =subjects(subjects_to_use_for_delta(resilient_permuted(sub,perm)));
        perm_indices_resilient_sessions = [perm_indices_resilient_sessions; find(strcmp(subj,this_sub))]; %append sessions
    end
    perm_indices_resilient.(['perm' num2str(perm)]) = perm_indices_resilient_sessions;
end

% vulnerable permutation groups
for perm = 1:10000
    perm_indices_vulnerable_sessions = [] ;
    for sub = 1:size(vulnerable_permuted,1)
        this_sub =subjects(subjects_to_use_for_delta(vulnerable_permuted(sub,perm)));
        perm_indices_vulnerable_sessions = [perm_indices_vulnerable_sessions; find(strcmp(subj,this_sub))]; %append sessions
    end
    perm_indices_vulnerable.(['perm' num2str(perm)]) = perm_indices_vulnerable_sessions;
end

%% run permutation of Maturational index
if run_permutations == 1
    %set up GLM
    clear slm
    MPC_baseline_res_perm = nan(360,360);
    MPC_baseline_vul_perm = nan(360,360);
    MPC_change_res_perm = nan(360,360);
    MPC_change_vul_perm = nan(360,360);
    MI_MPC_resilience_perm = nan(360,10000);
    MI_MPC_vulnerable_perm = nan(360,10000);
    MPC_baseline_res_perm_save = nan(360,360,10000);
    MPC_change_res_perm_save = nan(360,360,10000);
    group_diff = nan(360,10000);

    for perm = 1:10000%
        fprintf("\n start Permutation %d", perm)

        %use indices generated above
        indices_resilient_perm = perm_indices_resilient.(['perm' num2str(perm)]);
        indices_vulnerable_perm = perm_indices_vulnerable.(['perm' num2str(perm)]);

        %group 1
        this_age = age(indices_resilient_perm); Age         = term(this_age);
        this_sex = sex(indices_resilient_perm); Sex         = term(this_sex);
        this_subj = subj(indices_resilient_perm); Subj        = term(var2fac(str2double(this_subj)));
        this_site = site(indices_resilient_perm)';Site        = term(this_site);
        M1          = 1 + Age + Sex + Site + random(Subj) + I;

        parfor roi = 1:size(MPC,1)

            this_roi_mpc = squeeze(MPC(roi,:,indices_resilient_perm))';
            this_roi_mpc(:,roi)=[]; %remove nans to fit GLM
            this_roi_mpc(:,sum(this_roi_mpc)==0) = [];

            slm_mpc_model(roi) = SurfStatLinMod(this_roi_mpc, M1);
            slm_mpc(roi) = SurfStatT( slm_mpc_model(roi), age(indices_resilient_perm));
        end
        for roi = 1:size(MPC,1)
            MPC_baseline_res_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = (slm_mpc(roi).coef(1,:) + slm_mpc(roi).coef(2,:)*14  + (slm_mpc(roi).coef(3,:)*1/2) + (slm_mpc(roi).coef(4,:)*1/2) +(slm_mpc(roi).coef(5,:)*1/3)+ (slm_mpc(roi).coef(6,:)*1/3) + (slm_mpc(roi).coef(7,:)*1/3))';
            MPC_change_res_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = slm_mpc(roi).coef(2,:)';

            %compute MI
            [MI_MPC_resilience_perm(roi,perm), p_MI_MPC_resilience_perm(roi,perm)] = corr(MPC_baseline_res_perm(roi,1:length(slm_mpc(roi).coef(1,:)))',MPC_change_res_perm(roi,1:length(slm_mpc(roi).coef(1,:)))','type','Spearman','rows','complete');

            MPC_baseline_res_perm_save(roi,1:length(slm_mpc(roi).coef(1,:)),perm) =  MPC_baseline_res_perm(roi,1:length(slm_mpc(roi).coef(1,:)));
            MPC_change_res_perm_save(roi,1:length(slm_mpc(roi).coef(1,:)),perm) =  MPC_change_res_perm(roi,1:length(slm_mpc(roi).coef(1,:)));

        end

        clear M1 slm_mpc slm_mpc_model
        this_age = age(indices_vulnerable_perm); Age         = term(this_age);
        this_sex = sex(indices_vulnerable_perm); Sex         = term(this_sex);
        this_subj = subj(indices_vulnerable_perm); Subj        = term(var2fac(str2double(this_subj)));
        this_site = site(indices_vulnerable_perm)';Site        = term(this_site);
        M1          = 1 + Age + Sex + Site + random(Subj) + I;

        parfor roi = 1:size(MPC,1)

            this_roi_mpc = squeeze(MPC(roi,:,indices_vulnerable_perm))';
            this_roi_mpc(:,roi)=[]; %remove nans to fit GLM
            this_roi_mpc(:,sum(this_roi_mpc)==0) = [];

            slm_mpc_model(roi) = SurfStatLinMod(this_roi_mpc, M1);
            slm_mpc(roi) = SurfStatT( slm_mpc_model(roi), age(indices_vulnerable_perm));
        end
        for roi =  1:size(MPC,1)
            if numel(unique(this_site)) == 3
                MPC_baseline_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = (slm_mpc(roi).coef(1,:) + slm_mpc(roi).coef(2,:)*14  + (slm_mpc(roi).coef(3,:)*1/2) + (slm_mpc(roi).coef(4,:)*1/2) +(slm_mpc(roi).coef(5,:)*1/3)+ (slm_mpc(roi).coef(6,:)*1/3) + (slm_mpc(roi).coef(7,:)*1/3))';
                MPC_change_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = slm_mpc(roi).coef(2,:)';
            elseif numel(unique(this_site)) == 2 %in some permutations, one group only has subjects from 2 sites
                MPC_baseline_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = (slm_mpc(roi).coef(1,:) + slm_mpc(roi).coef(2,:)*14  + (slm_mpc(roi).coef(3,:)*1/2) + (slm_mpc(roi).coef(4,:)*1/2) +(slm_mpc(roi).coef(5,:)*1/2)+ (slm_mpc(roi).coef(6,:)*1/2))';
                MPC_change_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:))) = slm_mpc(roi).coef(2,:)';
            end

            %compute MI
            [MI_MPC_vulnerable_perm(roi,perm), p_MI_MPC_vulnerable_perm(roi,perm)] = corr( MPC_baseline_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:)))',MPC_change_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:)))','type','Spearman','rows','complete');

            MPC_baseline_vul_perm_save(roi,1:length(slm_mpc(roi).coef(1,:)),perm) =  MPC_baseline_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:)));
            MPC_change_vul_perm_save(roi,1:length(slm_mpc(roi).coef(1,:)),perm) =  MPC_change_vul_perm(roi,1:length(slm_mpc(roi).coef(1,:)));

        end
        group_diff(:,perm) = MI_MPC_resilience_perm(:,perm)-MI_MPC_vulnerable_perm(:,perm); %group 1 should be the bigger one resembling the bigger resilient group in res-vul

    end
    MPC_change_vul_perm_mean_save = squeeze(mean(MPC_change_vul_perm_save,2,'omitnan'));
    MPC_change_res_perm_mean_save = squeeze(mean(MPC_change_res_perm_save,2,'omitnan'));
    MPC_baseline_vul_perm_mean_save = squeeze(mean(MPC_baseline_vul_perm_save,2,'omitnan'));
    MPC_baseline_res_perm_mean_save = squeeze(mean(MPC_baseline_res_perm_save,2,'omitnan'));

    save([dataOutDir 'MI_delta_groups_MPC_balanced_null_models_surfstat_10000perms.mat'],'group_diff','MI_MPC_resilience_perm','MI_MPC_vulnerable_perm','perm_indices_resilient','perm_indices_vulnerable','MPC_change_res_perm','MPC_change_vul_perm','MPC_baseline_res_perm','MPC_baseline_vul_perm','MPC_baseline_res_perm_mean_save','MPC_change_res_perm_mean_save','MPC_baseline_vul_perm_mean_save','MPC_change_vul_perm_mean_save');

end
%% if skipped, load permuted and true group diff
%% assess significance of group differences
load([dataOutDir 'MI_delta_groups_MPC_balanced_null_models_surfstat_10000perms.mat'],'group_diff');
load([dataOutDir 'MI_delta_groups_MPC_res_groups.mat'],'MI_MPC_vulnerable','MI_MPC_resilient');

% group difference in
mpc_diff_resilient_vulnerable = MI_MPC_resilient - MI_MPC_vulnerable;

for i = 1:360
    alpha = 0.05; % significance level (5%)
    lower_quantile = quantile(group_diff(i,:), alpha / 2, 2); % lower quantile for each row
    upper_quantile = quantile(group_diff(i,:), 1 - alpha / 2, 2); % upper quantile for each row

    % Check if mpc_diff_resilient_vulnerable is in the top or bottom 2.5%
    below_bottom_quantile(i) = mpc_diff_resilient_vulnerable(i) < lower_quantile;
    above_top_quantile(i) = mpc_diff_resilient_vulnerable(i) > upper_quantile;
end
% combine for 2-sided testing based on 10.000 permutations
p_prop_05 = below_bottom_quantile+above_top_quantile;

% compute difference via z test (as has previously been done to test for
% group differences in MI, see Dorfschmidt et al. Science Advances
mpc_diff_resilient_vulnerable_p = nan(360,1);z_mpc_diff_resilient_vulnerable = nan(360,1);
for i = 1:360
    [mpc_diff_resilient_vulnerable_p(i),z_mpc_diff_resilient_vulnerable(i)] = compare_correlation_coefficients(MI_MPC_resilient(i),MI_MPC_vulnerable(i),360,360);
end
mpc_diff_resilient_vulnerable_p_fdr = fdr_bh(mpc_diff_resilient_vulnerable_p,0.05);
% only the overlap between both significance assessments is used as significance
mpc_diff_combined_threshold = mpc_diff_resilient_vulnerable_p_fdr'.*p_prop_05;

% plot group difference on the surface (FIGURE 4A)
figure_MI_mpc_diff = figure;  plot_cortical(parcel_to_surface(mpc_diff_resilient_vulnerable.*mpc_diff_combined_threshold,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','MI MPC res-vul')
exportfigbo(figure_MI_mpc_diff,[figDir 'MI_MPC_group_diff_combined_thresh_folded.png'],'png',10)

%% plot what directions reflect (Figure 4A bottom)
% resilient more disruptive
resilient_more_disruptive = mpc_diff_resilient_vulnerable'.* (MI_MPC.*p_MPC_fdr' < 0) .* (mpc_diff_resilient_vulnerable<0)' .* mpc_diff_combined_threshold';
directions_mpc_diff(resilient_more_disruptive~=0) = 1;
figure_resilient_more_disruptive = figure;  plot_cortical(parcel_to_surface(resilient_more_disruptive,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','resilient_more_disruptive')
exportfigbo(figure_resilient_more_disruptive,[figDir 'MI_MPC_resilient_more_disruptive_folded.png'],'png',10)

% resilient less conservative
resilient_less_conservative = mpc_diff_resilient_vulnerable'.* (MI_MPC.*p_MPC_fdr' > 0) .* (mpc_diff_resilient_vulnerable<0)' .* mpc_diff_combined_threshold';
figure_resilient_less_conservative = figure;  plot_cortical(parcel_to_surface(resilient_less_conservative,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','resilient_less_conservative')
exportfigbo(figure_resilient_less_conservative,[figDir 'MI_MPC_resilient_less_conservative_folded.png'],'png',10)

%differences at tipping points
resilient_tipping_points = mpc_diff_resilient_vulnerable'.* (p_MPC_fdr == 0)' .* (mpc_diff_resilient_vulnerable<0)' .* mpc_diff_combined_threshold';
figure_resilient_tipping_points = figure;  plot_cortical(parcel_to_surface(resilient_tipping_points,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','resilient_tipping_pointse')
exportfigbo(figure_resilient_tipping_points,[figDir 'MI_MPC_resilient_tipping_points.png'],'png',10)

%% plot scatter to show differences

these_parcels = find(resilient_less_conservative);

f=figure('Position',[50 50 2000 1000]);
for i = 1:numel(these_parcels)
    subplot(5,8,i)
    %roi = mpc_plot_examples(i)
    roi = these_parcels(i);s
    scatter(MPC_baseline_resilient(roi,:),MPC_change_resilient(roi,:),4,'filled','MarkerFaceColor',color_positive_delta)
    hold on
    scatter(MPC_baseline_vulnerable(roi,:),MPC_change_vulnerable(roi,:),4,'filled','MarkerFaceColor',color_negative_delta)

    p1 = polyfit(MPC_baseline_resilient(roi,:),MPC_change_resilient(roi,:), 1); % Linear trendline
    f1 = polyval(p1, MPC_baseline_resilient(roi,:), 1);
    plot(MPC_baseline_resilient(roi,:), f1, 'Color','black', 'LineWidth', 3,'Color',color_positive_delta);hold on
    p2 = polyfit(MPC_baseline_vulnerable(roi,:),MPC_change_vulnerable(roi,:), 1); % Linear trendline
    f2 = polyval(p2, MPC_baseline_vulnerable(roi,:), 1);
    plot(MPC_baseline_vulnerable(roi,:), f2, 'Color','black', 'LineWidth', 3,'Color',color_negative_delta);hold on
    xlim([0 2]); xticks([0 1 2]); ylim([-0.07 0.07]); yticks([-0.07 0 0.07])
    hold on
    title(nmmt{roi})

end

find(resilient_tipping_points~=0)
these_rois = [35 134 332];

f_selected_rois= figure('Position',[400 400 620 160]);
for i = 1:3
    subplot(1,3,i)
    %roi = mpc_plot_examples(i)
    roi = these_rois(i);
    scatter(MPC_baseline_resilient(roi,:),MPC_change_resilient(roi,:),3,'filled','MarkerFaceColor',color_positive_delta)
    hold on
    scatter(MPC_baseline_vulnerable(roi,:),MPC_change_vulnerable(roi,:),3,'filled','MarkerFaceColor',color_negative_delta)

    p1 = polyfit(MPC_baseline_resilient(roi,:),MPC_change_resilient(roi,:), 1); % Linear trendline
    f1 = polyval(p1, MPC_baseline_resilient(roi,:), 1);
    plot(MPC_baseline_resilient(roi,:), f1, 'Color','black', 'LineWidth', 3,'Color',color_positive_delta);hold on
    p2 = polyfit(MPC_baseline_vulnerable(roi,:),MPC_change_vulnerable(roi,:), 1); % Linear trendline
    f2 = polyval(p2, MPC_baseline_vulnerable(roi,:), 1);
    plot(MPC_baseline_vulnerable(roi,:), f2, 'Color','black', 'LineWidth', 3,'Color',color_negative_delta);hold on
    x_lim = [1.8 1.8 1.8];
    xlim([0 x_lim(i)]); xticks([0 x_lim(i)]); ylim([-0.07 0.07]); yticks([-0.07 0 0.07])
    title(nmmt{roi})

    hold on
end
exportfigbo(f_selected_rois,[figDir 'MI_diff_mpc_selected_rois.png'],'png',10)

%% Contextualize with cortical types and Yeo networks (Figure 4B)
%% yeo
yeo_counts = histcounts(yeo_glasser(logical(mpc_diff_combined_threshold)))
f_yeo_mpc = figure('Position',[400 400 500 500]);
donut(yeo_counts,labels_yeo , yeo_cmap);legend boxoff;legend('Location','eastoutside');axis off
exportfigbo(f_yeo_mpc, [figDir 'donut_mpc_yeo.png'],'png',10)

%% cyto
labels_cyto = {'Konicortex','Eulaminate-III','Eulaminate-II','Eulaminate-I','Dysgranular','Agranular'}
cyto_counts = histcounts(types360(logical(mpc_diff_combined_threshold)));
f_cyto_mpc = figure('Position',[400 400 500 500]);
donut(cyto_counts,labels_cyto , cmap_cortical_types(2:end,:));legend boxoff;legend('Location','eastoutside');axis off
exportfigbo(f_cyto_mpc, [figDir 'donut_mpc_cyto.png'],'png',10)

%% cross-modal maturational modes
mpc_diff_combined_threshold_extr = mpc_diff_combined_threshold(included_parcels);
mode_counts = histcounts(directions(logical(mpc_diff_combined_threshold_extr)));
f_modes_mpc = figure('Position',[400 400 500 500]);
donut(mode_counts,{'none', '+mpc/+fc','+mpc/-fc','-mpc/+fc','-mpc/-fc'},category_cmap);legend boxoff;legend('Location','eastoutside');axis off
exportfigbo(f_modes_mpc, [figDir 'donut_mpc_modes.png'],'png',10)

%% export table with directions for all significant parcels
diff = mpc_diff_resilient_vulnerable;
z = z_mpc_diff_resilient_vulnerable;
p_perm = mpc_diff_resilient_vulnerable_p;

%not yet computed directions:
resilient_less_disruptive = mpc_diff_resilient_vulnerable'.* (MI_MPC.*p_MPC_fdr' < 0) .* (mpc_diff_resilient_vulnerable>0)' .* mpc_diff_combined_threshold';
resilient_more_conservative = mpc_diff_resilient_vulnerable'.* (MI_MPC.*p_MPC_fdr' > 0) .* (mpc_diff_resilient_vulnerable>0)' .* mpc_diff_combined_threshold';

directions_mpc_diff = zeros(360,1);
directions_mpc_diff(resilient_more_disruptive~=0) = 1;
directions_mpc_diff(resilient_less_conservative~=0) = 2;
directions_mpc_diff(resilient_tipping_points~=0) = 3;
directions_mpc_diff(resilient_less_disruptive~=0) = 4;
directions_mpc_diff(resilient_more_conservative~=0) = 5;

vulnerable_more_disruptive = (MI_MPC.*p_MPC_fdr' < 0) .* (mpc_diff_resilient_vulnerable>0)' .* mpc_diff_combined_threshold';

T_directions = [array2table(nmmt,'VariableNames',{'ROI'}) array2table([diff' z p_perm directions_mpc_diff],'VariableNames',{'delta','z','pperm','direction'})];
T_export = T_directions(logical(mpc_diff_combined_threshold),:);
writetable(T_export,[dataOutDir 'mpc_diff_ROIs.csv'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       FC                      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_group_MIs == 1
    %% true group difference
    %% FC
    %resilient group
    this_age = age(resilient_indices); Age         = term(this_age);
    this_sex = sex(resilient_indices); Sex         = term(this_sex);
    this_subj = subj(resilient_indices); Subj        = term(var2fac(str2double(this_subj)));
    this_site = site(resilient_indices)';Site        = term(this_site);
    M_resilient          = 1 + Age + Sex + Site + random(Subj) + I;

    FC_resilience = FC(:,:,resilient_indices);
    FC_vulnerable = FC(:,:,vulnerable_indices);

    clear slm_FC slm_FC_model MI_FC_resilient p_MI_FC_resilient  MI_FC_vulnerable p_MI_FC_vulnerable
    parfor roi = 1:size(FC_resilience,1) %takes approx 30 mins; with parfor 20 workers: 2-3 min

        this_roi_FC = squeeze(FC_resilience(roi,:,:))';
        this_roi_FC(:,roi)=[]; %remove nans to fit GLM
        this_roi_FC(:,sum(this_roi_FC)==0) = []; %check if any columns are all zero and remove to fit GLM

        slm_FC_model(roi) = SurfStatLinMod(this_roi_FC, M_resilient);
        slm_FC(roi) = SurfStatT( slm_FC_model(roi), age(resilient_indices));
    end
    for roi = 1:size(FC_resilience,1)
        FC_baseline_resilient(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/3)+ (slm_FC(roi).coef(6,:)*1/3) + (slm_FC(roi).coef(7,:)*1/3))';
        FC_change_resilient(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';

        %compute MI
        [MI_FC_resilient(roi), p_MI_FC_resilient(roi)] = corr(FC_baseline_resilient(roi,1:length(slm_FC(roi).coef(1,:)))',FC_change_resilient(roi,1:length(slm_FC(roi).coef(1,:)))','type','Spearman','rows','complete');
    end

    %vulnerable group
    this_age = age(vulnerable_indices); Age         = term(this_age);
    this_sex = sex(vulnerable_indices); Sex         = term(this_sex);
    this_subj = subj(vulnerable_indices); Subj        = term(var2fac(str2double(this_subj)));
    this_site = site(vulnerable_indices)';Site        = term(this_site);
    M_vulnerable          = 1 + Age + Sex + Site + random(Subj) + I;

    clear slm_FC slm_FC_model
    parfor roi = 1:size(FC_vulnerable,1) %takes approx 30 mins; with parfor 20 workers: 2-3 min

        this_roi_FC = squeeze(FC_vulnerable(roi,:,:))';
        this_roi_FC(:,roi)=[]; %remove nans to fit GLM
        this_roi_FC(:,sum(this_roi_FC)==0) = []; %check if any columns are all zero and remove to fit GLM

        slm_FC_model(roi) = SurfStatLinMod(this_roi_FC, M_vulnerable);
        slm_FC(roi) = SurfStatT( slm_FC_model(roi), age(vulnerable_indices));
    end
    for roi = 1:size(FC_vulnerable,1)
        FC_baseline_vulnerable(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/3)+ (slm_FC(roi).coef(6,:)*1/3) + (slm_FC(roi).coef(7,:)*1/3))';
        FC_change_vulnerable(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';

        %compute MI
        [MI_FC_vulnerable(roi), p_MI_FC_vulnerable(roi)] = corr(FC_baseline_vulnerable(roi,1:length(slm_FC(roi).coef(1,:)))',FC_change_vulnerable(roi,1:length(slm_FC(roi).coef(1,:)))','type','Spearman','rows','complete');
    end
    save([dataOutDir 'MI_delta_groups_FC_res_groups.mat'],'MI_FC_vulnerable','p_MI_FC_vulnerable','MI_FC_resilient','p_MI_FC_resilient','FC_baseline_resilient','FC_change_resilient','FC_baseline_vulnerable','FC_change_vulnerable');
end
% load precomputed data
load([dataOutDir 'MI_delta_groups_FC_res_groups.mat']);

%% permutation of random groups
if run_permutations == 1
    clear slm
    for perm = 1:10000%
        fprintf("\n start Permutation %d", perm)

        %use indices generated above
        indices_resilient_perm = perm_indices_resilient.(['perm' num2str(perm)]);
        indices_vulnerable_perm = perm_indices_vulnerable.(['perm' num2str(perm)]);

        %group 1 (+delta SRS)
        this_age = age(indices_resilient_perm); Age         = term(this_age);
        this_sex = sex(indices_resilient_perm); Sex         = term(this_sex);
        this_subj = subj(indices_resilient_perm); Subj        = term(var2fac(str2double(this_subj)));
        this_site = site(indices_resilient_perm)';Site        = term(this_site);
        M1          = 1 + Age + Sex + Site + random(Subj) + I;

        parfor roi = 1:size(FC,1)
            this_roi_FC = squeeze(FC(roi,:,indices_resilient_perm))';
            this_roi_FC(:,roi)=[]; %remove nans to fit GLM
            this_roi_FC(:,sum(this_roi_FC)==0) = [];

            slm_FC_model(roi) = SurfStatLinMod(this_roi_FC, M1);
            slm_FC(roi) = SurfStatT( slm_FC_model(roi), age(indices_resilient_perm));
        end
        for roi = 1:size(FC,1)
            if numel(unique(this_site)) == 3
                FC_baseline_res_perm(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/3)+ (slm_FC(roi).coef(6,:)*1/3) + (slm_FC(roi).coef(7,:)*1/3))';
                FC_change_res_perm(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';
            elseif numel(unique(this_site)) == 2
                FC_baseline_res_perm(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/2)+ (slm_FC(roi).coef(6,:)*1/2))';
                FC_change_res_perm(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';
            end
            %compute MI
            [MI_FC_res_perm(roi,perm), p_MI_FC_res_perm(roi,perm)] = corr(FC_baseline_res_perm(roi,1:length(slm_FC(roi).coef(1,:)))',FC_change_res_perm(roi,1:length(slm_FC(roi).coef(1,:)))','type','Spearman','rows','complete');

        end
        
        % group 2 (=-delta SRS)
        this_age = age(indices_vulnerable_perm); Age         = term(this_age);
        this_sex = sex(indices_vulnerable_perm); Sex         = term(this_sex);
        this_subj = subj(indices_vulnerable_perm); Subj        = term(var2fac(str2double(this_subj)));
        this_site = site(indices_vulnerable_perm)';Site        = term(this_site);
        M1          = 1 + Age + Sex + Site + random(Subj) + I;

        clear slm_FC slm_FC_model
        parfor roi = 1:size(FC,1)

            this_roi_FC = squeeze(FC(roi,:,indices_vulnerable_perm))';
            this_roi_FC(:,roi)=[]; %remove nans to fit GLM
            this_roi_FC(:,sum(this_roi_FC)==0) = [];

            slm_FC_model(roi) = SurfStatLinMod(this_roi_FC, M1);
            slm_FC(roi) = SurfStatT( slm_FC_model(roi), age(indices_vulnerable_perm));
        end
        for roi =  1:size(FC,1)
            if numel(unique(this_site)) == 3
                FC_baseline_vul_perm(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/3)+ (slm_FC(roi).coef(6,:)*1/3) + (slm_FC(roi).coef(7,:)*1/3))';
                FC_change_vul_perm(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';
            elseif numel(unique(this_site)) == 2
                FC_baseline_vul_perm(roi,1:length(slm_FC(roi).coef(1,:))) = (slm_FC(roi).coef(1,:) + slm_FC(roi).coef(2,:)*14  + (slm_FC(roi).coef(3,:)*1/2) + (slm_FC(roi).coef(4,:)*1/2) +(slm_FC(roi).coef(5,:)*1/2)+ (slm_FC(roi).coef(6,:)*1/2))';
                FC_change_vul_perm(roi,1:length(slm_FC(roi).coef(1,:))) = slm_FC(roi).coef(2,:)';
            end
            %compute MI
            [MI_FC_vulnerable_perm(roi,perm), p_MI_FC_vulnerable_perm(roi,perm)] = corr( FC_baseline_vul_perm(roi,1:length(slm_FC(roi).coef(1,:)))',FC_change_vul_perm(roi,1:length(slm_FC(roi).coef(1,:)))','type','Spearman','rows','complete');

        end
        group_diff_fc(:,perm) = MI_FC_res_perm(:,perm)-MI_FC_vulnerable_perm(:,perm); %group 1 should be the bigger one resembling the bigger resilient group in res-vul

    end
    save([dataOutDir 'MI_delta_groups_FC_null_models_surfstat_10000perms.mat'],'group_diff_fc','MI_FC_res_perm','MI_FC_vulnerable_perm','perm_indices_resilient','perm_indices_vulnerable','FC_baseline_res_perm','FC_baseline_vul_perm','FC_change_res_perm','FC_change_vul_perm');
end

%% if skipped, load permuted and true group diff
%% assess significance of group differences in FC-MI
load([dataOutDir 'MI_delta_groups_FC_null_models_surfstat_10000perms.mat'],'group_diff_fc');
load([dataOutDir 'MI_delta_groups_FC_res_groups.mat'],'MI_FC_vulnerable','MI_FC_resilient');

FC_diff_resilient_vulnerable = MI_FC_resilient - MI_FC_vulnerable;

for i = 1:330
    alpha = 0.05; % significance level (5%)
    lower_quantile = quantile(group_diff_fc(i,:), alpha / 2, 2); % lower quantile for each row
    upper_quantile = quantile(group_diff_fc(i,:), 1 - (alpha / 2), 2); % upper quantile for each row

    % Check if mpc_diff_resilient_vulnerable is in the top or bottom 2.5%
    below_bottom_quantile_fc(i) = FC_diff_resilient_vulnerable(i) < lower_quantile;
    above_top_quantile_fc(i) = FC_diff_resilient_vulnerable(i) > upper_quantile;
end
% combine for 2-sided testing
p_prop_fc_05 = below_bottom_quantile_fc+above_top_quantile_fc;

% z-test
for i = 1:330
    [fc_diff_resilient_vulnerable_p(i),z_fc_diff_resilient_vulnerable(i)] = compare_correlation_coefficients(MI_FC_resilient(i),MI_FC_vulnerable(i),330,330);
end

% only show parcels as significant if they survive both z-test and
% permutation
fc_group_diff_p_combined = p_prop_fc_05.*fdr_bh(fc_diff_resilient_vulnerable_p,0.05); %only use parcels that are significant both in z-test and permutation

cmap_fc_mi = [0.5 0.5 0.5; cmap_curv;[0.5 0.5 0.5]/256];
to_plot = parcel_330_to_360(FC_diff_resilient_vulnerable,'zer');

to_plot(not_included_parcels) = -6;%add mask for excluded parcels
%unthresholded map
to_plot_thresh = parcel_330_to_360(FC_diff_resilient_vulnerable.*fc_group_diff_p_combined,'zer');
to_plot_thresh(not_included_parcels) = -6;%add mask for excluded parcels

figure_fc_folded = figure; plot_cortical(parcel_to_surface(to_plot,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','MI FC diff unthresholded')
exportfigbo(figure_fc_folded,[figDir 'MI_FC_group_diff_unthresholded.png'],'png',10)
% thresholded map
figure_fc_thresh = figure; plot_cortical(parcel_to_surface(to_plot_thresh,'glasser_360_conte69'),'color_range',[-0.8 0.8],'surface_name','conte69','label_text','MI FC diff thresholded')
exportfigbo(figure_fc_thresh,[figDir 'MI_FC_group_diff_combined_thresh_folded.png'],'png',10)

%% Donuts
%FC
fc_categ = [sum(fc_group_diff_p_combined(directions==0)) sum(fc_group_diff_p_combined(directions==1)) sum(fc_group_diff_p_combined(directions==2)) sum(fc_group_diff_p_combined(directions==3)) sum(fc_group_diff_p_combined(directions==4))];
%overlay: A) regions with significant group difference and B) (significant)
%positive and negative MIs from the whole-sample MI
FC_diff_thresh = FC_diff_resilient_vulnerable.*fc_group_diff_p_combined;
ov = zeros(330,1);
conservative_MI_fc_pos = ov; conservative_MI_fc_pos( logical((FC_diff_thresh>0)' .* (MI_FC_fdr>0)) )=1;
disruptive_MI_fc_pos = ov; disruptive_MI_fc_pos( logical((FC_diff_thresh>0)' .* (MI_FC_fdr<0)) )=1;
conservative_MI_fc_neg = ov; conservative_MI_fc_neg( logical((FC_diff_thresh<0)' .* (MI_FC_fdr>0)) )=1;
disruptive_MI_fc_neg = ov; disruptive_MI_fc_neg( logical((FC_diff_thresh<0)' .* (MI_FC_fdr<0)) )=1;
MI_fc_none = ov; MI_fc_none( logical((FC_diff_thresh~=0)' .* (MI_FC_fdr==0)) )=1;

fc_modes = [sum(MI_fc_none) sum(conservative_MI_fc_pos) sum(disruptive_MI_fc_pos) sum(conservative_MI_fc_neg) sum(disruptive_MI_fc_neg)];

f_category_modes = figure('Position',[400 400 500 500]);
donut(fc_modes, {'none','+ \Delta more conservative','- \Delta more disruptive','- \Delta more conservative','+ \Delta more disruptive'},mode_cmap);legend boxoff;legend('Location','eastoutside');axis off
exportfigbo(f_category_modes, [figDir 'donut_fc_modes.png'],'png',10)

f_category_donut = figure('Position',[400 400 500 500]);
donut(fc_categ, {'none','+mpc/+fc','+mpc/-fc','-mpc/+fc','-mpc/-fc'},category_cmap);legend boxoff;legend('Location','eastoutside');xticks([]);yticks([]);axis off
exportfigbo(f_category_donut, [figDir 'donut_fc_categories.png'],'png',10)

%% %%%%%%%%%%%%%%% %%
%%    gradients    %%
%% %%%%%%%%%%%%%%% %%

if compute_gradients == 1
    %mean grad to allign everyone to
    meanMPC2 = squeeze(mean(MPC_repeated,3));
    meanMPC2(eye(size(meanMPC2))==1) = 0;

    %% age gradient
    % Set up model
    this_age = all_demo_tbl_repeated.age; Age         = term(this_age);
    this_sex = all_demo_tbl_repeated.sex; Sex         = term(this_sex);
    this_subj = all_demo_tbl_repeated.subj; Subj        = term(this_subj);
    this_site = all_demo_tbl_repeated.site; Site        = term(this_site);
    M1          = 1 + Age + Sex + Site + random(Subj) + I;
    % run for positive and negative age effects
    age_effect_pos_neg = [1, -1];
    age_effect_label = {'positive','negative'};

    MPC2_long = [];
    parfor sub = 1:size(MPC_repeated,3)
        MPC2_long(sub,:) = squareform(MPC_repeated(:,:,sub)); %from matrix to vector
    end
    idx = sum(MPC2_long)==0; % identify zero columns

    clear slm
    slm = SurfStatLinMod(MPC2_long(:,~idx), M1);
    slm = SurfStatT( slm, this_age);

    age_t_MPC = zeros(1, size(MPC_repeated,2));
    age_t_MPC(~idx) = slm.t;
    age_t_MPC_all = squareform(age_t_MPC); %use for alignment

    gm = GradientMaps('kernel','na','approach','dm'); % initializes gradient computation
    gm = gm.fit(age_t_MPC_all);

    mpc_age_grad  = gm.gradients{1}(:,1); %first gradient of person i, procrustes

    clear gm
    plot_this = parcel_to_surface(mpc_age_grad+10,'glasser_360_conte69');
    cmap_grad_mid = [0.9 0.9 0.9; cmap_grad];
    f_agegrad = figure;  plot_cortical(plot_this,'surface_name','conte69','label_text','positive age effect whole sample','color_range',[10-0.18 10+0.12]);
    exportfigbo(f_agegrad, [figDir 'MPC_age_grad_all.png'],'png',8)

    clear age_t_MPC

    %% compute age gradients for each group
    for i = 1:2%length(groups_names)
        for age_effect = 1%:2 %positive or negative
            group_idx = group_names(i);
            this_age = age(indices_res_groups.(group_names{i})); Age         = term(this_age);
            this_sex = sex(indices_res_groups.(group_names{i})); Sex         = term(this_sex);
            this_subj = subj(indices_res_groups.(group_names{i})); Subj        = term(this_subj);
            this_site = site(indices_res_groups.(group_names{i}))'; Site        = term(this_site);
            M1          = 1 + Age + Sex + Site + random(Subj) + I;

            MPC2_long = [];
            parfor sub = 1:size(MPC,3)
                MPC2_long(sub,:) = squareform(MPC(:,:,sub)); %from matrix to vector
            end
            MPC2_long = MPC2_long(indices_res_groups.(group_names{i}),:);
            idx = sum(MPC2_long)==0; % identify zero columns
            clear slm
            slm = SurfStatLinMod(MPC2_long(:,~idx), M1);
            slm = SurfStatT( slm, this_age*age_effect_pos_neg(age_effect) );

            pos_q = SurfStatQ( slm);
            age_t_MPC_tmp = zeros(1, size(MPC,2));
            age_t_MPC_tmp(~idx) = slm.t;
            age_t_MPC(:,:,age_effect) = squareform(age_t_MPC_tmp);

            clear gm

            gm = GradientMaps('kernel','na','approach','dm','align','pa'); % initializes gradient computation
            gm = gm.fit({age_t_MPC_all(:,:,age_effect),age_t_MPC(:,:,age_effect)});

            mpc_age_grad(:,i,age_effect) = gm.aligned{2}(:,1);
        end
    end
    for age_effect = 1%:2 %positive or negative
        %plot age change gradients
        f_agegrad_resilient = figure;  plot_cortical(parcel_to_surface(mpc_age_grad(:,1,age_effect)+10,'glasser_360_conte69'),'surface_name','conte69','label_text',[age_effect_label{age_effect} ' age effect resilient group'],'color_range',[10-0.17 10+0.12]);
        exportfigbo(f_agegrad_resilient, [figDir 'MPC_' age_effect_label{age_effect} '_age_effect_resilient_group.png'],'png',8)

        f_agegrad_vulnerable = figure; plot_cortical(parcel_to_surface(mpc_age_grad(:,2,age_effect)+10,'glasser_360_conte69'),'surface_name','conte69','label_text',[age_effect_label{age_effect} ' age effect vulnerable group'],'color_range',[10-0.17 10+0.12]);
        exportfigbo(f_agegrad_vulnerable, [figDir 'MPC_' age_effect_label{age_effect} '_age_effect_vulnerable_group.png'],'png',8)
        close all
    end

    f_grad = figure('Position',[400 400 600 400]);
    subplot(2,2,1);
    scatter(baseline(1,:)',mpc_age_grad(:,1,1),10,color_positive_delta,'filled') % 1=resilient group
    xlim([-0.2 0.2]);ylim([-0.2 0.2]);axis square
    xlabel('baseline MPC gradient');ylabel([age_effect_label{1} ' age effect gradient'])
    title('+ delta')
    subplot(2,2,2);

    scatter(baseline(2,:)',mpc_age_grad(:,2,1),10,color_negative_delta,'filled') % 2 = vulnerable group
    xlim([-0.2 0.2]);ylim([-0.2 0.2]); axis square
    title('- delta');
    xlabel('baseline MPC gradient')

    exportfigbo(f_grad, [figDir 'MPC_gradients_age_effect_resilience_scatter.png'],'png',8)

    % Plot the density function
    f_dist = figure('position',[400 400 300 280]);
    [f_res, x_res] = ksdensity(mpc_age_grad(:,1,1));
    [f_vul, x_vul] = ksdensity(mpc_age_grad(:,2,1));

    plot(x_res, f_res, 'LineWidth', 2, 'Color', color_positive_delta);
    hold on
    plot(x_vul, f_vul, 'LineWidth', 2, 'Color', color_negative_delta);
    box off;set(gca,'tickdir','none','YTick',[0 5])

    %xlabel('Gradient loading');
    %ylabel('Density');

    exportfigbo(f_dist, [figDir 'distribution_gradient.png'],'png',10)

end

%% Supplementary analyses: testing MPC-MI group differences via interactions 
%% MPC

% 1) interaction model
combined_groups = ~isnan(demo_tbl.resilience_group);
resilience_group_coded = repmat(cellstr('none'),numel(combined_groups),1);
resilience_group_coded(find(demo_tbl.resilience_group==1)) = cellstr('resilient');
resilience_group_coded(find(demo_tbl.resilience_group==0)) = cellstr('vulnerable');

demo_tbl_interaction = [table(age,'VariableNames', {'age'}), table(sex,'VariableNames',{'sex'}), table(site','VariableNames',{'site'}),...
    table(subj,'VariableNames',{'subj'}), table(resilience_group_coded,'VariableNames',{'Resilience'})];
demo_tbl_interaction = demo_tbl_interaction(combined_groups,:);
MPC_interaction = MPC(:,:,combined_groups);
[MI_mpc_interaction, p_MI_mpc_interaction] = compute_MI_lme_interaction(MPC_interaction, demo_tbl_interaction, 'Resilience_vulnerable');

diff_int = MI_mpc_interaction.Resilience_vulnerable_0-MI_mpc_interaction.Resilience_vulnerable_1;
z_diff_int_p = nan(360,1);z_diff_int = nan(360,1);
for i = 1:360
    [z_diff_int_p(i),z_diff_int(i)] = compare_correlation_coefficients(MI_mpc_interaction.Resilience_vulnerable_0(i),MI_mpc_interaction.Resilience_vulnerable_1(i),360,360);
end
z_diff_int_p_fdr = fdr_bh(z_diff_int_p,0.05);

% visualize results
figure; plot_cortical(parcel_to_surface(diff_int.*z_diff_int_p_fdr,'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8])
figure; plot_cortical(parcel_to_surface(diff_int,'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8],'label_text','interaction model')
figure; plot_cortical(parcel_to_surface(diff_int,'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8],'label_text','original model')

%correlate this with map in main manuscript = 0.966

% Scatter plot
f_int = figure('Position',[200 200 280 280]);
scatter(mpc_diff_resilient_vulnerable, diff_int, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7]);

% Linear fit
p = polyfit(mpc_diff_resilient_vulnerable, diff_int, 1);
x_fit = linspace(min(mpc_diff_resilient_vulnerable), max(mpc_diff_resilient_vulnerable), 100);
y_fit = polyval(p, x_fit);

% Plotting the linear fit line
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
xlim([-1.2 1.2]);ylim([-1.2 1.2])
xticks([-1 0 1]);yticks([-1 0 1])

% Labels and legend
xlabel('Separate Models');
ylabel('Combined Model');
exportfigbo(f_int,[figDir 'interaction_model_scatter.png'],'png',10)

% permutation of random groups
demo_tbl_interaction_perm = [table(age,'VariableNames', {'age'}), table(sex,'VariableNames',{'sex'}), table(site','VariableNames',{'site'}),...
    table(subj,'VariableNames',{'subj'}), table(resilience_group_coded,'VariableNames',{'Resilience'})];

if run_permutations == 1
    for perm = 1:1000%do only 1000 for supplementaries
        fprintf("\n start Permutation %d", perm)

        %use indices generated above
        indices_resilient_perm = perm_indices_resilient.(['perm' num2str(perm)]);
        indices_vulnerable_perm = perm_indices_vulnerable.(['perm' num2str(perm)]);
        indices_combined = [indices_resilient_perm; indices_vulnerable_perm];

        data = MPC(:,:,indices_combined);
        demographics_tbl = demo_tbl_interaction_perm(indices_combined,:);

        [MI_mpc_interaction_perm, p_MI_mpc_interaction] = compute_MI_lme_interaction(data, demographics_tbl, 'Resilience_vulnerable');

        diff_int_perm(:,perm) = MI_mpc_interaction.Resilience_vulnerable_0-MI_mpc_interaction.Resilience_vulnerable_1;
    end
    save([dataOutDir 'Supplementaries_MI_delta_groups_interaction_MPC_null_models_surfstat_10000perms.mat'],'diff_int_perm');
end

%% compute MIs via pearson's correlation
load([dataOutDir 'MI_delta_groups_MPC_res_groups.mat'],'MI_MPC_vulnerable','p_MI_MPC_vulnerable','MI_MPC_resilient','p_MI_MPC_resilient','MPC_baseline_resilient','MPC_change_resilient','MPC_baseline_vulnerable','MPC_change_vulnerable')

for roi = 1:360
    MI_res_pearson(roi) = corr(MPC_baseline_resilient(roi,:)',MPC_change_resilient(roi,:)','type','Pearson');
    MI_vul_pearson(roi) = corr(MPC_baseline_vulnerable(roi,:)',MPC_change_vulnerable(roi,:)','type','Pearson');
end

f_pearson = figure;plot_cortical(parcel_to_surface(MI_res_pearson-MI_vul_pearson,'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8])
exportfigbo(f_pearson,[figDir 'MI_diff_pearson.png'],'png',10)

% Scatter plot
f_pear = figure('Position',[200 200 280 280]);
scatter(mpc_diff_resilient_vulnerable, MI_res_pearson-MI_vul_pearson, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7]);

% Linear fit
p = polyfit(mpc_diff_resilient_vulnerable, MI_res_pearson-MI_vul_pearson, 1);
x_fit = linspace(min(mpc_diff_resilient_vulnerable), max(mpc_diff_resilient_vulnerable), 100);
y_fit = polyval(p, x_fit);

% Plotting the linear fit line
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
xlim([-1.2 1.2]);ylim([-1.2 1.2])
xticks([-1 0 1]);yticks([-1 0 1])

% Label and legend
xlabel('MI Spearman');
ylabel('MI Pearson');
exportfigbo(f_pear,[figDir 'pearson_MI_diff_scatter.png'],'png',10)

%% plot association between MPC MI and MPC age gradient
f_pear = figure('Position',[200 200 280 280]);
scatter(mpc_diff_resilient_vulnerable, MI_res_pearson-MI_vul_pearson, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7]);

% Linear fit
p = polyfit(mpc_diff_resilient_vulnerable, MI_res_pearson-MI_vul_pearson, 2);
x_fit = linspace(min(mpc_diff_resilient_vulnerable), max(mpc_diff_resilient_vulnerable), 100);
y_fit = polyval(p, x_fit);

% Plotting the linear fit line
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
xlim([-1.2 1.2]);ylim([-1.2 1.2])
xticks([-1 0 1]);yticks([-1 0 1])

%% FC interaction model
load([dataOutDir 'MI_delta_groups_FC_res_groups.mat'])
FC_interaction = FC(:,:,combined_groups);
[MI_fc_interaction, p_MI_fc_interaction] = compute_MI_lme_interaction(FC_interaction, demo_tbl_interaction, 'Resilience_vulnerable');

diff_int_fc = MI_fc_interaction.Resilience_vulnerable_0-MI_fc_interaction.Resilience_vulnerable_1;%from interaction
fc_diff_resilient_vulnerable = MI_FC_resilient - MI_FC_vulnerable;

figure_int = figure; plot_cortical(parcel_to_surface(parcel_330_to_360(diff_int_fc,'zer'),'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8],'label_text','interaction model')
figure_orig = figure; plot_cortical(parcel_to_surface(parcel_330_to_360(fc_diff_resilient_vulnerable,'zer'),'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8],'label_text','original model')

% Scatter plot
f_int = figure('Position',[200 200 280 280]);
scatter(fc_diff_resilient_vulnerable, diff_int_fc, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7]);
[r_comb,p_comb]= corr(fc_diff_resilient_vulnerable', diff_int_fc);
% Linear fit
p = polyfit(fc_diff_resilient_vulnerable, diff_int_fc, 1);
x_fit = linspace(min(fc_diff_resilient_vulnerable), max(fc_diff_resilient_vulnerable), 100);
y_fit = polyval(p, x_fit);
% Plotting the linear fit line
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
xlim([-1.2 1.2]);ylim([-1.2 1.2])
xticks([-1 0 1]);yticks([-1 0 1])
% Labels and legend
xlabel('Separate Models');
ylabel('Combined Model');
exportfigbo(f_int,[figDir 'interaction_model_scatter_fc.png'],'png',10)

% compute MI via pearson correlaton
for roi = 1:330
    MI_res_pearson_fc(roi) = corr(FC_baseline_resilient(roi,:)',FC_change_resilient(roi,:)','type','Pearson');
    MI_vul_pearson_fc(roi) = corr(FC_baseline_vulnerable(roi,:)',FC_change_vulnerable(roi,:)','type','Pearson');
end

f_pearson_fc = figure;plot_cortical(parcel_to_surface(parcel_330_to_360(MI_res_pearson_fc-MI_vul_pearson_fc,'zer'),'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.8 0.8])
exportfigbo(f_pearson_fc,[figDir 'MI_diff_pearson_fc.png'],'png',10)

% Scatter plot
f_pear = figure('Position',[200 200 280 280]);
scatter(fc_diff_resilient_vulnerable, MI_res_pearson_fc-MI_vul_pearson_fc, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7]);
[r_sp, p_sp] = corr(fc_diff_resilient_vulnerable', (MI_res_pearson_fc-MI_vul_pearson_fc)');
% Linear fit
p = polyfit(fc_diff_resilient_vulnerable, MI_res_pearson_fc-MI_vul_pearson_fc, 1);
x_fit = linspace(min(fc_diff_resilient_vulnerable), max(fc_diff_resilient_vulnerable), 100);
y_fit = polyval(p, x_fit);

% Plotting the linear fit line
hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
xlim([-1.2 1.2]);ylim([-1.2 1.2])
xticks([-1 0 1]);yticks([-1 0 1])
xlabel('MI Spearman');
ylabel('MI Pearson');
exportfigbo(f_pear,[figDir 'pearson_MI_diff_scatter_fc.png'],'png',10)

