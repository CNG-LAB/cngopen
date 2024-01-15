%% Linking longitudinal change in susceptibility/resilience 
%% to psychosocial stressors to myeloarchitectural trajectories
% This script consists of 3 main parts:
% 1) characterizing MT change in an adolescent sample of n = 199
% 2) linking MT change to change in resilience capacities
% 3) linking functional connectivity change to change in resilience
% capacities

% define and add your paths here
baseDir     = '';
dataDir     = '';
dataOutDir  = '';
toolboxDir  = '';
figDir      = '';
cd(dataOutDir)

addpath(genpath(baseDir))
addpath(genpath(dataDir))
addpath(genpath(dataOutDir))
% required toolboxes. 
% MPC and surfstat addons were downloaded from: https://github.com/MICA-MNI/micaopen
addpath(genpath('/data/pt_02792/Toolboxes/Resources/micaopen-master/MPC/'))
addpath(genpath([toolboxDir 'surfstat/']))
addpath(genpath([toolboxDir 'ENIGMA-1.1.3']))
addpath(genpath([toolboxDir 'BrainSpace-0.1.2']))
addpath(genpath([toolboxDir 'fdr_bh']))
addpath(genpath([toolboxDir 'cbrewer/']))
addpath(genpath([toolboxDir 'Resources/micaopen-master/surfstat/surfstat_addons/']))
addpath(genpath('/data/p_02792/Resources/'))
addpath(genpath([toolboxDir 'Violinplot-Matlab-master/']))

%% load data
% data
load([dataOutDir 'included_parcels_fc.mat'])
load([dataOutDir '/NSPN_fc_mt_ct_subsample.mat'])%for 512 sessions from 295 individuals
load([dataOutDir 'Delta_10000_permutation_indices.mat'])

%% load atlases, colormaps, labels, and other useful things
% atlases / parcellations
load('yeo_glasser360.mat');
load('cortical_types_360.mat');
parc = readmatrix([toolboxDir 'ENIGMA-1.1.3/enigmatoolbox/datasets/parcellations/glasser_360_conte69.csv']);
% labels
labels_yeo = {'Visual','Somatomotor','Dorsal attention','Ventral attention','Limbic','Frontoparietal','Default mode'};
labels_types = {'Konicortex','Eulaminate3','Eulaminate2','Eulaminate1','Dysgranular','Agranular'};
% colormaps
RR = [6,117,111]/255;
RR3 = [150,150,150]/255;
RR4 = [0.8 0.8 0.8];
load([toolboxDir 'ScientificColourMaps7/roma/roma.mat'])
roma = flipud(roma);
cmap_resilience = [ones(1,3).*0.9; cbrewer('seq','BuPu',255)];cmap_resilience(cmap_resilience<0)=0;
color_positive_delta = cmap_resilience(240,:);
color_negative_delta = cmap_resilience(100,:);

cmap_curv = [roma(2:127,:);(ones(1,3).*0.95); (ones(1,3).*0.8); roma(130:end,:)];
load('yeo_colormap.mat'); yeo_cmap_surface = yeo_cmap; yeo_cmap(1,:)=[];%remove gray for 0 parcels/midline
cmap_norm = cbrewer('seq','BuPu',255);cmap_norm(cmap_norm<0)=0;
cmap_grad1 = cbrewer('seq','BuPu',100); cmap_grad1(cmap_grad1<0)=0;

%% settings
winsorize = 1;
plot_general_principles = 1;
run_supplementaries= 1;
parpooling = 0;

if parpooling ==1
    parpool('Processes',30)
end

%% extract delta change in structure & function
MTmean = squeeze(mean(MT,1));

subjects = unique(subj);%295
sex_rep =cell(numel(subjects),1);
age_delta = nan(numel(subjects),1);
age_first_session = nan(numel(subjects),1);
site_rep = cell(numel(subjects),1);
delta_mt_mean = nan(360,numel(subjects));
fc_baseline = nan(346,346,numel(subjects));
fc_second_session = nan(346,346,numel(subjects));
delta_mt = nan(10,360,numel(subjects));
mt_baseline = nan(10,360,numel(subjects));
mt_second_session = nan(10,360,numel(subjects));
mt_mean_baseline = nan(360,numel(subjects));
subs_with_rep_mri_sessions = nan(numel(subjects),1);
delta_fc = nan(346,346, numel(subjects));
delta_ct = nan(360, numel(subjects));
mean_age = nan(numel(subjects),1);

for sub = 1:length(subjects) %loop through 295
    idx = find(strcmp(subj,subjects(sub)));%find subject index in 512 sessions
    sex_rep(sub)=sex(idx(1));
    site_rep(sub) = site(idx(1));
    age_first_session(sub) = age(idx(1));
    mt_baseline(:,:,sub) = MT(:,:,idx(1));
    mt_mean_baseline(:,sub) = MTmean(:,idx(1));
    fc_baseline(:,:,sub) = FC(:,:,idx(1));

    if length(idx)>1
        subs_with_rep_mri_sessions(sub)=1; %which of the 295 subjects have a follow-up?
        age_delta(sub) = age(idx(end))-age(idx(1));
        mean_age(sub) = mean([age(idx(end)) age(idx(1))]);
        % extract 2nd session for plotting
        mt_second_session(:,:,sub) = MT(:,:,idx(end));
        fc_second_session(:,:,sub) = FC(:,:,idx(end));

        % deltas
        delta_mt(:,:,sub) = (MT(:,:,idx(end))-MT(:,:,idx(1)))/ (age(idx(end))-age(idx(1)));
        delta_mt_mean(:,sub) = (MTmean(:,idx(end))-MTmean(:,idx(1))) / (age(idx(end))-age(idx(1)));
        delta_fc(:,:,sub) = (FC(:,:,idx(end))-FC(:,:,idx(1))) / (age(idx(end))-age(idx(1)));
    else
        subs_with_rep_mri_sessions(sub)=0;
    end
end

%% Extract resilience scores - baseline & change
resilience_scores = readtable([dataOutDir '/resilience_based_on_pfactor/resilience_scores_prediction_repeatedmeasures_totalscores.csv']);
sub_exist = zeros(length(subjects),1);
sub_repeated_behav = zeros(length(subjects),1);
mean_resilience_score = nan(length(subjects),1);
baseline_resilience_score = nan(length(subjects),1);
delta_resilience_score = nan(length(subjects),1);
end_resilience_score = nan(length(subjects),1);

for sub = 1:numel(subjects)
    if ~isnan(resilience_scores{str2double(subjects{sub}) == resilience_scores{:,1}, 7})
        scores = resilience_scores{str2double(subjects{sub}) == resilience_scores{:,1}, 7};
        mean_resilience_score(strcmp(subjects,subjects(sub))) = mean(scores);
        baseline_resilience_score(strcmp(subjects,subjects(sub))) = scores(1);
        if length(scores)>1
            delta_resilience_score(strcmp(subjects,subjects(sub))) = scores(end)-scores(1);
            end_resilience_score(strcmp(subjects,subjects(sub))) = scores(end);
            sub_repeated_behav(sub) = 1;
        end
        sub_exist(sub) = 1; %whether behavioral data is available for this individual
    else
        sub_exist(sub) = 0;
    end
end
delta_resilience_group = double(delta_resilience_score>0);

%% winsorize deltas
delta_mt_winsorized = nan(size(delta_mt_mean));
delta_fc_winsorized = nan(size(delta_fc));
if winsorize ==1
    for roi = 1:360
        delta_mt_winsorized(roi,:) =  delta_mt_mean(roi,:);
        sd = std(delta_mt_winsorized(roi,:),'omitnan');
        m = mean(delta_mt_winsorized(roi,:),'omitnan');
        delta_mt_winsorized(roi,delta_mt_winsorized(roi,:)<(m-3*sd)) = m-3*sd;
        delta_mt_winsorized(roi,delta_mt_winsorized(roi,:)>(m+3*sd)) = m+3*sd;
    end
    for roi = 1:346
        for iroi = 1:346
            delta_fc_winsorized(roi,iroi,:) = squeeze(delta_fc(roi,iroi,:));
            sd = std(squeeze(delta_fc_winsorized(roi,iroi,:)),'omitnan');
            m = mean(squeeze(delta_fc_winsorized(roi,iroi,:)),'omitnan');
            if sum(squeeze(delta_fc_winsorized(roi,iroi,:))<(m-3*sd))
                delta_fc_winsorized(roi,iroi,logical(squeeze(delta_fc_winsorized(roi,iroi,:))<(m-3*sd))) = m-3*sd;
            elseif sum(squeeze(delta_fc_winsorized(roi,iroi,:))>(m+3*sd))
                delta_fc_winsorized(roi,iroi,logical(squeeze(delta_fc_winsorized(roi,iroi,:))>(m+3*sd))) = m+3*sd;
            end
        end
    end
    %resilience delta
    sd = std(delta_resilience_score,'omitnan');
    m = mean(delta_resilience_score,'omitnan');
    delta_resilience_score_winsorized = delta_resilience_score;
    delta_resilience_score_winsorized(delta_resilience_score_winsorized < (m-3*sd)) = m-3*sd;
    delta_resilience_score_winsorized(delta_resilience_score_winsorized > (m+3*sd)) = m+3*sd;
    % continue with winsorized data
    delta_mt_mean = delta_mt_winsorized;
    delta_fc = delta_fc_winsorized;
    delta_resilience_score = delta_resilience_score_winsorized;
end

baseline_table = [array2table(subjects,'VariableNames',{'subjects'}) array2table(age_first_session,'VariableNames',{'baseline_age'}) table(sex_rep,'VariableNames',{'sex'}) array2table(site_rep,'VariableNames',{'site'}) array2table(baseline_resilience_score,'VariableNames',{'baseline_resilience'})];
delta_table = [array2table(subjects,'VariableNames',{'subjects'}) array2table(mean_age,'VariableNames',{'mean_age'}) table(sex_rep,'VariableNames',{'sex'}) array2table(site_rep,'VariableNames',{'site'})  array2table(delta_resilience_score,'VariableNames',{'delta_resilience'})];

% extract data of individuals who have both repeated MRIs and repeated
% behavioral data available
subjects_to_use_for_delta = intersect(find(sub_repeated_behav), find(subs_with_rep_mri_sessions)');
subject_codes_to_use_for_delta = subjects(subjects_to_use_for_delta);
delta_table_extracted = [delta_table(subjects_to_use_for_delta,:) array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'}),...
    array2table(delta_resilience_group(subjects_to_use_for_delta),'VariableNames',{'delta_resilience_group'})];
delta_mt_mean_extracted=delta_mt_mean(:,subjects_to_use_for_delta);
delta_resilience = delta_table_extracted.delta_resilience;

%save data
save([dataOutDir 'delta_table_extracted.mat'],'baseline_table','delta_table','delta_table_extracted')
save([dataOutDir 'subjects_to_use_for_delta.mat'],'subjects_to_use_for_delta','subject_codes_to_use_for_delta')
%% Plot delta MT average and STD
% Including the total sample (i.e., all individuals with repeated mri sessions) to show general principles in a larger sample

if plot_general_principles ==1
    %mean delta mt
    f_mt = figure;meike_plot_cortical(parcel_to_surface(mean(delta_mt_mean,2,'omitnan'),'glasser_360_conte69'),'surface_name','conte69','label_text','Delta MT','cmap_custom',cmap_norm,'color_range',[0,20])
    exportfigbo(f_mt,[figDir, 'Delta_MT_mean.png'],'png', 12)
    %std delta mt
    f_mt_sd = figure;meike_plot_cortical(parcel_to_surface(std(delta_mt_mean,[],2,'omitnan'),'glasser_360_conte69'),'surface_name','conte69','label_text','Delta MT std','cmap_custom',cmap_norm,'color_range',[0,100])
    exportfigbo(f_mt_sd,[figDir, 'Delta_MT_std.png'],'png', 12)

    % plot 3 age strata
    strata14_17 = age_first_session<17';
    strata17_21 = age_first_session>17' & age_first_session<21';
    strata21_24 =age_first_session >21';
    age_group = [strata14_17 strata17_21 strata21_24];
    age_label = {'<17y', '17-21y', '>21y' };
    age_strata_colors = cbrewer('seq','BuPu',9);age_strata_colors = [age_strata_colors(9,:); age_strata_colors(5,:); age_strata_colors(3,:)];

    f_agestrata = figure('Position',[400 400 480 280]);
    for i = 1:3
        subplot(1,3,i)
        tp1 = mean(squeeze(mean(mt_baseline(:,:,age_group(:,i)),1)),1)';
        tp2 = mean(squeeze(mean(mt_second_session(:,:,age_group(:,i)),1)),1)';
        tp1(isnan(tp2))=[]; tp2(isnan(tp2))=[];
        for sub=1:size(tp1,1)
            scatter([1,2],[tp1(sub),tp2(sub)],20,RR4,'filled','d'); xlim([0.8, 2.2]); ylim([760, 960])
            plot([1,2],[tp1(sub),tp2(sub)],'Color',RR4); xlim([0.8, 2.2]); ylim([750, 950]);set(gca,'XTick',[1,2],'XTickLabel',{'T1','T2'});xlabel(age_label{i});box off
            hold on
        end
        plot([1,2],[mean(tp1),mean(tp2)],'Color',age_strata_colors(i,:),'LineWidth',3);
        scatter([1,2],[mean(tp1),mean(tp2)],40,age_strata_colors(i,:),'filled','d'); xlim([0.8, 2.2]); ylim([760, 960])

        if i>1
            set(gca,'YTick',[],'FontSize',8,'FontName','Seaford')

        elseif i<2
            ylabel('MT','FontName','Seaford','FontSize',10)
        end
        hold on
    end
    hold off
    exportfigbo(f_agestrata,[figDir, 'MT_age_strata.png'],'png', 12)

    % plot layer-wise normative MT change
    f_layers = figure('Position',[400 400 450 280]);
    subplot(1,2,1)
    for i = 1:3
        tp2 = squeeze(mean(mt_second_session(:,:,age_group(:,i)),2))';
        out = isnan(tp2(:,1));
        tp2 = mean(tp2(~out,:),1)';
        tp1 = squeeze(mean(mt_baseline(:,:,age_group(:,i)),2))';
        tp1 = mean(tp1(~out,:),1)';
        hold on
        plot(tp2-tp1,1:10,'Color',age_strata_colors(i,:),'LineWidth',2)
        xlabel('Mean \Delta MT','FontSize',10,'FontName','Seaford')
        ylabel('Cortical depth','FontSize',10,'FontName','Seaford')
        set(gca, 'XLim',[1 25], 'XTick',[1 25],'FontSize',8,'FontName','Seaford')
        set(gca, 'YLim',[1 10],'YTick',[])%, 'YTickLabel',{'White', 'Pial'},'FontSize',8,'FontName','Seaford')
    end
    exportfigbo(f_layers,[figDir, 'MT_age_strata_layers.png'],'png', 12)

    % plot layer-wise normative MT change std
    subplot(1,2,2)
    for i = 1:3
        tp2 = squeeze(mean(mt_second_session(:,:,age_group(:,i)),2,'omitnan'));
        out = isnan(tp2(1,:));
        tp2 = tp2(:,~out);
        tp1 = squeeze(mean(mt_baseline(:,:,age_group(:,i)),2,'omitnan'));
        tp1 = tp1(:,~out);
        hold on
        plot(std(tp2'-tp1'),1:10,'Color',age_strata_colors(i,:),'LineWidth',2)
        xlabel('SD \Delta MT','FontSize',10,'FontName','Seaford')
        %ylabel('Cortical depth','FontSize',10,'FontName','Seaford')
        set(gca, 'XLim',[20 40], 'XTick',[20 40],'FontSize',8,'FontName','Seaford')
        set(gca, 'YLim',[1 10],'YTick',[]);
    end
    exportfigbo(f_layers,[figDir, 'MT_age_strata_layers.png'],'png', 12)

    %% Normative MT change gradient
    % captures a large scale axis of coordinated MT change patterns,
    % vomputed via diffusion map embedding as implemented in BrainSpace

    delta_mt_r = corr(delta_mt_mean(:,logical(subs_with_rep_mri_sessions))');
    gm_mt = GradientMaps('kernel','na','approach','dm','n_components', 10);
    gm_mt = gm_mt.fit(delta_mt_r);
    mt_gradient = gm_mt.gradients{1}(:,1);
    plot_grad = parcel_to_surface(mt_gradient,'glasser_360_conte69');plot_grad(parc==0)=-100;%-100 to mask midbrain
    f_mtgrad = figure;meike_plot_cortical(plot_grad,'surface_name','conte69','label_text','MT change gradient','color_range',[-0.1 0.1],'cmap_custom',cmap_grad1)
    exportfigbo(f_mtgrad, [figDir 'MT_change_gradient.png'],'png', 12)
    %variance explained:
    %     var_ex = gm_mt.lambda{1}(1) / sum(gm_mt.lambda{1});
    %     handles = scree_plot(gm_mt.lambda{1});

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             MT Change * Resilience      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 2Ai) Delta MT and change in resilience
% Set up model for the true effect
Age = delta_table_extracted.mean_age; Age = term(Age);
Sex = delta_table_extracted.sex; Sex = term(Sex);
Site = delta_table_extracted.site; Site = term(Site);
Mean_resilience = delta_table_extracted.mean_resilience_score; Mean_resilience = term(Mean_resilience);
Delta_resilience = delta_table_extracted.delta_resilience; Delta_resilience = term(Delta_resilience);

M1 = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;

slm = SurfStatLinMod(delta_mt_mean_extracted', M1);
slm = SurfStatT( slm, delta_resilience);
mt_res_effect = slm.t;
figure;plot_cortical(parcel_to_surface(mt_res_effect,'glasser_360_conte69'),'surface_name','conte69','label_text','Delta Resilience * Delta Mean MT','color_range',[-5,5])

% create permutation indices or load 'Delta_10000_permutation_indices' for
% exact reproducibility
% for perm = 1:10000
%     perm_idx(:,perm) = randperm(numel(subjects_to_use_for_delta),numel(subjects_to_use_for_delta));
% end
% save([dataOutDir 'Delta_10000_permutation_indices.mat'], 'perm_idx')

mt_res_effect_perm = nan(360,10000);
for perm = 1:10000
    permuted_delta_resilience = delta_resilience(perm_idx(:,perm));
    permuted_mean_resilience = delta_table_extracted.mean_resilience_score(perm_idx(:,perm));
    Delta_resilience = term(permuted_delta_resilience);
    Mean_resilience = term(permuted_mean_resilience);

    M1 = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;
    clear slm
    slm = SurfStatLinMod(delta_mt_mean_extracted', M1);
    slm = SurfStatT( slm, permuted_delta_resilience);
    mt_res_effect_perm(:,perm) = slm.t;

    %fprintf(['Permutation ' num2str(perm)])
end

save([dataOutDir 'Delta_permutations_10000.mat'], 'mt_res_effect_perm','perm_idx','mt_res_effect')
load([dataOutDir 'Delta_permutations_10000.mat']) %if computation is skipped

p_perm_fdr = get_perm_fdr_p(mt_res_effect,mt_res_effect_perm,360);

f = figure;meike_plot_cortical(parcel_to_surface(mt_res_effect.*p_perm_fdr','glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_resilience,'label_text','Delta Resilience * Delta Mean MT','color_range',[0,4])
exportfigbo(f,[figDir, 'Delta_Resilience_Delta_MT_10000perm05.png'],'png', 12)

% plot direction of effects and individual age trends
region_base_tp2 = [mean(squeeze(mean(mt_baseline(:,logical(p_perm_fdr),subjects_to_use_for_delta),1)),1)' mean(squeeze(mean(mt_second_session(:,logical(p_perm_fdr),subjects_to_use_for_delta),1)),1)'];
age_base_tp2 = [age_first_session (age_first_session+age_delta)];
age_base_tp2= age_base_tp2(subjects_to_use_for_delta,:);
timepoints = age_base_tp2 - mean(age_base_tp2,2);

extr_delta_resilience_score = delta_resilience_score(subjects_to_use_for_delta,:);

%winsorize data for plotting (as they are also winsorized for the model
mm1 = mean(region_base_tp2(:,1));
ss1 = std(region_base_tp2(:,1));
mm2 = mean(region_base_tp2(:,2));
ss2 = std(region_base_tp2(:,2));

f_lines = figure('Position',[400 400 450 280]);
subplot(1,2,1)
for sub = 1:length(age_base_tp2)
    %winsorize
    if region_base_tp2(sub,1) < (mm1-(3*ss1))
        region_base_tp2(sub,1) = mm1-(3*ss1);
    elseif region_base_tp2(sub,1) > (mm1+(3*ss1))
        region_base_tp2(sub,1) = mm1+(3*ss1);
    end
    if region_base_tp2(sub,2) < (mm2-(3*ss2))
        region_base_tp2(sub,2) = mm2-(3*ss2);
    elseif region_base_tp2(sub,2) > (mm2+(3*ss2))
        region_base_tp2(sub,2) = mm2+(3*ss2);
    end
    if extr_delta_resilience_score(sub)>0
        col = [224, 200, 217]/256;
    elseif extr_delta_resilience_score(sub)<0
        col = cmap_resilience(40,:);
    end
    scatter(timepoints(sub,:), region_base_tp2(sub,:),10,col,'filled');
    plot(timepoints(sub,:), region_base_tp2(sub,:), 'Color',col,'LineWidth',0.7);
    hold on
end
%plot group means:
p2 = polyfit(timepoints(extr_delta_resilience_score<0,:), region_base_tp2(extr_delta_resilience_score<0,:), 1); % Linear trendline
f2 = polyval(p2, [-1.2 1.2]);
plot([-1.2 1.2], f2, 'Color',[0 0 0], 'LineWidth', 2,'Color',color_negative_delta,'LineWidth',3);

p1 = polyfit(timepoints(extr_delta_resilience_score>0,:), region_base_tp2(extr_delta_resilience_score>0,:), 1); % Linear trendline
f1 = polyval(p1, [-1.2 1.2]);
plot([-1.2 1.2], f1, 'Color',[0 0 0], 'LineWidth', 2,'Color',color_positive_delta,'LineWidth',3);

xlabel({'Study visit','(centered, years)'},'FontName','Seaford','FontSize',1)
ylabel('Mean MT in PFC cluster','FontName','Seaford','FontSize',10)
set(gca,'FontSize',8,'FontName','Seaborn')
ylim([650 1050]);xlim([-1.2 1.2])
box off

%% Figure 2Aii) layer-specific effects within significant parcels
sig_mt_rois = find(p_perm_fdr);
layer_effect = nan(10,numel(sig_mt_rois)); layer_effectp = nan(10,numel(sig_mt_rois));
for lay=1:10
    for roi = 1:sum(p_perm_fdr)
        parcel_mean = squeeze(delta_mt(lay,sig_mt_rois(roi),subjects_to_use_for_delta));
        this_tbl = [array2table(parcel_mean,'VariableNames',{'mri'}) delta_table(subjects_to_use_for_delta,:), array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
        lme = fitlme(this_tbl,'mri ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
        layer_effect(lay,roi) =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
        layer_effectp(lay,roi) =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));
    end
end

subplot(1,2,2)
for roi = 1:numel(sig_mt_rois)
    plot(layer_effect(:,roi), 1:10, 'LineWidth', 2,'Color',RR4);
    xlabel('t','FontName','Seaford')
    %ylabel('Cortical depth','FontName','Arial')
    set(gca, 'XLim',[2 5], 'XTick',[2 5],'FontSize',8,'FontName','Seaborn')
    set(gca, 'YLim',[1 10],'YTick',[1 10], 'YTickLabel',{'White', 'Pial'},'FontName','Seaford','FontSize',10)
    hold on
    plot(mean(layer_effect,2), 1:10, 'LineWidth', 3,'Color',color_positive_delta);
    box off
end
exportfigbo(f_lines, [figDir 'layer_wise_mt_change_resilience.png'],'png', 12)

%or as boxplots
%h = boxplot(layer_effect', 'Colors', RR, 'Widths', 0.7,'Symbol','','orientation', 'horizontal');
%set(h, 'LineWidth', 2);box off;xticks([0 2.5 5]), xlim([0 5])

%% supplementary analyses:
% 1) are there cross-sectional baseline effects of change in SRS or SRS
%{
Age = delta_table_extracted.mean_age; Age = term(Age);
Sex = delta_table_extracted.sex; Sex = term(Sex);
Site = delta_table_extracted.site; Site = term(Site);
Mean_resilience = delta_table_extracted.mean_resilience_score; Mean_resilience = term(Mean_resilience);
Delta_resilience = delta_table_extracted.delta_resilience; Delta_resilience = term(Delta_resilience);
Baseline_resilience = baseline_resilience_score(subjects_to_use_for_delta); Baseline_resilience = term(Baseline_resilience);

M1 = 1+ Delta_resilience + Mean_resilience + Age + Sex + Site;
M2 = 1+ Baseline_resilience + Age + Sex + Site;

mt_baseline_extracted = squeeze(mean(mt_baseline(:,:,subjects_to_use_for_delta),1));
mt_res_effect_perm = nan(360,10000,2);

for models = 1:2
    if models == 1
        model = M1;
        contrast = delta_table_extracted.delta_resilience;
    elseif models == 2
        model = M2;
        contrast = baseline_resilience_score(subjects_to_use_for_delta);
    end
    %true effect
    slm = SurfStatLinMod(mt_baseline_extracted', model);
    slm = SurfStatT( slm, contrast);
    mt_res_effect_base(:,models) = slm.t;

    for perm = 1:10000
        permuted_delta_resilience = delta_resilience(perm_idx(:,perm)); %load perm_idx from above for exact reproduction of figures
        permuted_mean_resilience = delta_table_extracted.mean_resilience_score(perm_idx(:,perm));
        base_tmp = baseline_resilience_score(subjects_to_use_for_delta);
        permuted_baseline_resilience = base_tmp(perm_idx(:,perm));
        Delta_resilience = term(permuted_delta_resilience);
        Mean_resilience = term(permuted_mean_resilience);
        Baseline_resilience = term(permuted_baseline_resilience);

         clear slm
        if models == 1
            model = 1+ Delta_resilience + Mean_resilience + Age + Sex + Site;
            slm = SurfStatLinMod(mt_baseline_extracted', model);
            slm = SurfStatT( slm, permuted_delta_resilience);
        elseif models == 2
           model = 1+ Baseline_resilience + Age + Sex + Site;
            slm = SurfStatLinMod(mt_baseline_extracted', model);
            slm = SurfStatT( slm, permuted_baseline_resilience);
        end

        mt_res_effect_base_perm(:,perm,models) = slm.t;

    end

    p_perm_base_fdr(:,models) = get_perm_fdr_p(mt_res_effect_base(:,models),squeeze(mt_res_effect_base_perm(:,:,models)),360);
end

f_base = figure;meike_plot_cortical(parcel_to_surface(mt_res_effect_base(:,1).*p_perm_base_fdr(:,1),'glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_resilience,'label_text','Delta Resilience * Baseline Mean MT','color_range',[-5,5])
exportfigbo(f_base, [figDir 'delta_resilience_baseline_mt.png'],'png', 10)

f_base2 = figure;meike_plot_cortical(parcel_to_surface(mt_res_effect_base(:,2),'glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_resilience,'label_text','Baseline Resilience * Baseline Mean MT','color_range',[-5,5])
exportfigbo(f_base2, [figDir 'baseline_resilience_baseline_mt_unthresholded.png'],'png', 10)

2) are results robust to subsampling?

for perm = 1:1000
    perm_idx_sub(:,perm) = randperm(numel(subjects_to_use_for_delta),113);
    %113 = 80% of the data
end
mt_res_effect_sub = nan(360,1000);
for perm = 1:1000
    permuted_delta_resilience = delta_resilience(perm_idx_sub(:,perm)); %load perm_idx from above for exact reproduction of figures
    permuted_mean_resilience = delta_table_extracted.mean_resilience_score(perm_idx_sub(:,perm));
    permuted_age = delta_table_extracted.mean_age(perm_idx_sub(:,perm));
    permuted_sex = delta_table_extracted.sex(perm_idx_sub(:,perm));
    permuted_site = delta_table_extracted.site(perm_idx_sub(:,perm));
    Delta_resilience = term(permuted_delta_resilience);
    Mean_resilience = term(permuted_mean_resilience);
    Age_perm =  term(permuted_age);
    Sex_perm = term(permuted_sex);
    Site_perm = term(permuted_site);

    delta_mt_perm_sub = delta_mt_mean_extracted(:,perm_idx_sub(:,perm));

    M1 = 1 + Delta_resilience + Mean_resilience + Age_perm + Sex_perm + Site_perm;
    clear slm
    slm = SurfStatLinMod(delta_mt_perm_sub', M1);
    slm = SurfStatT( slm, permuted_delta_resilience);
    mt_res_effect_sub(:,perm) = slm.t;

    %fprintf(['Permutation ' num2str(perm)])
end

%check correlation between original t-map and t-map derived from
sub-sampling
for perm = 1:size(mt_res_effect_sub,2)
    corr_subs(perm) = corr(mt_res_effect_sub(:,perm),mt_res_effect');
end

% Plot distribution of correlations
f_sub = figure('Position',[400 400 280 280]);
h=histogram(corr_subs)
h.NumBins = 30
h.FaceColor = [0.7, 0.7, 0.7]
h.LineWidth = 0.5
xlabel('r');
xlim([0 1])
xticks([0:0.2:1])
ylabel('Frequency');
box off
hold off;
exportfigbo(f_sub,[figDir 'distribution_subsampling_delta.png'],'png',10)

mt_effet_perm_avg = max(mt_res_effect_sub(logical(p_perm_fdr),:)); %mean
%across parcels within significant cluster

%f_base2 = figure;meike_plot_cortical(parcel_to_surface(mt_effet_perm_avg,'glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_resilience,'label_text','Baseline Resilience * Baseline Mean MT','color_range',[-5,5])
avg_true = max(mt_res_effect(logical(p_perm_fdr)));

% Plot distribution
f_sub = figure('Position',[400 400 280 280]);
[f,xi,bw] = ksdensity(mt_effet_perm_avg);
plot(xi,f,'LineWidth', 2, 'Color', [0.7, 0.7, 0.7])
hold on
yLimits = ylim;
xlim([1 7])
plot([avg_true, avg_true], [0, yLimits(2)], 'k-', 'LineWidth', 2); % Solid black line% Set font properties
set(gca, 'FontSize', 10, 'FontName', 'Roboto');% Add labels and legend
xlabel('Max t-value');
ylabel('Density');
box off
hold off;
exportfigbo(f_sub,[figDir 'distribution_subsampling_res_effect.png'],'png',10)

[h,p,ci,zval]=ztest(avg_true,mean(mt_effet_perm_avg),std(mt_effet_perm_avg))
%}
%% FC connection from MT cluster to yeo networks
sig_mt_330 = p_perm_fdr(included_parcels);
delta_fc_ctx = delta_fc(17:end,17:end,:);
delta_mt_fc = squeeze(mean(delta_fc_ctx(logical(sig_mt_330),:,subjects_to_use_for_delta),1,'omitnan')); %regions that are in the same yeo network will be nans
yeo_included = yeo_glasser(included_parcels);

yeo_cmap_surface_masked = yeo_cmap_surface;
yeo_cmap_surface_masked(9,:) = [0,0,0];
%% Figure 2C - seed based FC effect

% main analysis
Age = delta_table_extracted.mean_age; Age = term(Age);
Sex = delta_table_extracted.sex; Sex = term(Sex);
Site = delta_table_extracted.site; Site = term(Site);
Mean_resilience = delta_table_extracted.mean_resilience_score; Mean_resilience = term(Mean_resilience);
Delta_resilience = delta_table_extracted.delta_resilience; Delta_resilience = term(Delta_resilience);

%true effect
M1 = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;
clear slm
slm = SurfStatLinMod(delta_mt_fc', M1);
slm = SurfStatT( slm, delta_resilience);
seed_effect = slm.t;

%permutation
fc_res_effect_perm = nan(330,10000);
for perm = 1:10000
    permuted_delta_resilience = delta_resilience(perm_idx(:,perm));
    permuted_mean_resilience = delta_table_extracted.mean_resilience_score(perm_idx(:,perm));
    Delta_resilience = term(permuted_delta_resilience);
    Mean_resilience = term(permuted_mean_resilience);

    M1 = 1 + Delta_resilience + Mean_resilience + Age + Sex + Site;
    clear slm
    slm = SurfStatLinMod(delta_mt_fc', M1);
    slm = SurfStatT( slm, permuted_delta_resilience);
    fc_res_effect_perm(:,perm) = slm.t;

    %fprintf(num2str(perm))
end

fc_p_prop_pos = nan(330,1); fc_p_prop_neg = nan(330,1);
for i = 1:330
    a = seed_effect(i); %true t-values
    b = fc_res_effect_perm(i,:); %t-values from permutation
    fc_p_prop_pos(i) = sum(b<a) / size(fc_res_effect_perm(i,:),2);
    fc_p_prop_neg(i) = sum(b>a) / size(fc_res_effect_perm(i,:),2);
end
seed_fdr = fdr_bh(fc_p_prop_pos,0.025) + fdr_bh(fc_p_prop_neg,0.025); %two tailed testing

%plot thresholded seed effects
f_seed = figure;meike_plot_cortical(parcel_to_surface(parcel_330_to_360(seed_effect.*seed_fdr,'zer'),'glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_resilience,'color_range',[0, 4])
exportfigbo(f_seed, [figDir 'seed_based_delta_fc_resilience.png'],'png', 12);
%add PFC mask
to_plot = parcel_330_to_360(seed_effect.*seed_fdr','zer');
not_included_parcels = 1:360; not_included_parcels(included_parcels)=[];
to_plot(logical(p_perm_fdr))=6;
%add mask for excluded parcels
to_plot(not_included_parcels) = -6;

cmap_seeds = [repmat([0.5 0.5 0.5],256,1); cmap_resilience;[0.2 0.2 0.2]/256];
f_seed_masked = figure;meike_plot_cortical(parcel_to_surface(to_plot,'glasser_360_conte69'),'surface_name','conte69','color_range',[-4, 4],'cmap_custom', cmap_seeds)
exportfigbo(f_seed_masked, [figDir 'seed_based_delta_fc_resilience_masked_against_10000perms.png'],'png', 12);

%% supplementary figure: plot general delta FC for this cluster
cmap_fc = [0.5 0.5 0.5; flipud(cbrewer('div','RdBu',256));[0.2 0.2 0.2]/256];cmap_fc(cmap_fc<0)=0;cmap_fc(cmap_fc>1)=1;
to_plot = parcel_330_to_360(mean(delta_mt_fc,2),'zer');
to_plot(logical(p_perm_fdr))=6;
%add mask for excluded parcels
to_plot(not_included_parcels) = -6;
f_fc = figure;meike_plot_cortical(parcel_to_surface(to_plot,'glasser_360_conte69'),'surface_name','conte69','color_range',[-0.04 0.04],'cmap_custom', cmap_fc)
exportfigbo(f_fc, [figDir 'delta_fc-group_average.png'],'png', 12);


%% global FC change of the cluster in delta positive vs negative resilience (Figure 2i)
% extract global FC and resilience scores at first and last timepoints to plot age trends
region_base_tp2 = [mean(squeeze(mean(fc_baseline(:,logical(p_perm_fdr),subjects_to_use_for_delta),1,'omitnan')),1)' mean(squeeze(mean(fc_second_session(:,logical(p_perm_fdr),subjects_to_use_for_delta),1,'omitnan')),1)'];
delta_fc_PFC = mean(squeeze(mean(delta_fc_ctx(:,logical(p_perm_fdr),subjects_to_use_for_delta),1,'omitnan')),1)';

% statistical global resilience effect?
lme = fitlme([table(delta_fc_PFC) delta_table_extracted],'delta_fc_PFC ~ delta_resilience + mean_resilience_score + mean_age + sex + site');
FC_glob_p = lme.Coefficients.pValue(strcmp(lme.CoefficientNames,'delta_resilience'));

%plot lines
mm1 = mean(region_base_tp2(:,1));
ss1 = std(region_base_tp2(:,1));
mm2 = mean(region_base_tp2(:,2));
ss2 = std(region_base_tp2(:,2));

age_base_tp2 = [age_first_session (age_first_session+age_delta)];
age_base_tp2= age_base_tp2(subjects_to_use_for_delta,:);
timepoints = age_base_tp2 - mean(age_base_tp2,2);

f_lines_fc = figure('Position',[400 400 450 280]);
subplot(1,2,1)
for sub = 1:length(age_base_tp2)
    %winsorize
    if region_base_tp2(sub,1) < (mm1-(3*ss1))
        region_base_tp2(sub,1) = mm1-(3*ss1);
    elseif region_base_tp2(sub,1) > (mm1+(3*ss1))
        region_base_tp2(sub,1) = mm1+(3*ss1);
    end
    if region_base_tp2(sub,2) < (mm2-(3*ss2))
        region_base_tp2(sub,2) = mm2-(3*ss2);
    elseif region_base_tp2(sub,2) > (mm2+(3*ss2))
        region_base_tp2(sub,2) = mm2+(3*ss2);
    end

    if extr_delta_resilience_score(sub)>0
        col = [224, 200, 217]/256;
    elseif extr_delta_resilience_score(sub)<0
        col = cmap_resilience(40,:);
    end
    scatter(timepoints(sub,:), region_base_tp2(sub,:),10,col,'filled');
    plot(timepoints(sub,:), region_base_tp2(sub,:), 'Color',col,'LineWidth',0.7);
    hold on
    %plot group means:
    p1 = polyfit(timepoints(extr_delta_resilience_score>0,:), region_base_tp2(extr_delta_resilience_score>0,:), 1); % Linear trendline
    f1 = polyval(p1, [-1.2 1.2]);
    plot([-1.2 1.2], f1, 'Color',[0 0 0], 'LineWidth', 2,'Color',color_positive_delta,'LineWidth',3);
    hold on
    p2 = polyfit(timepoints(extr_delta_resilience_score<0,:), region_base_tp2(extr_delta_resilience_score<0,:), 1); % Linear trendline
    f2 = polyval(p2, [-1.2 1.2]);
    plot([-1.2 1.2], f2, 'Color',[0 0 0], 'LineWidth', 2,'Color',color_negative_delta,'LineWidth',3);

    hold on
    xlabel({'Study visit','(centered, years)'},'FontName','Seaford','FontSize',1)
    ylabel('Global FC in PFC cluster','FontName','Seaford','FontSize',10)
    ylim([-0.1 0.8])

    set(gca,'FontSize',8,'FontName','Seaborn')
    box off
end
exportfigbo(f_lines_fc,[figDir 'fc_resilience_globaldc.png'],'png', 12)

%% decoding figure panels, Figure 2D
% decoding of fc seed effect --> Yeo networks
f_violin_decoding = figure('Position',[400 400 1200 280]);
subplot(1,3,1)
% decoding of mt effect --> Cortical types
cyto_cmap_surface = cmap_cortical_types;
min_sig = min(mt_res_effect(sig_mt_rois));
subplot(1,3,1)
plot([0 7],[min_sig min_sig],'Color','black');hold on
for class = 1:6
    mt_plot_this.(labels_types{class}) = mt_res_effect(types360==class)  ;
end
violinplot(mt_plot_this,labels_types,'ViolinColor',cmap_cortical_types(2:end,:),'MarkerSize',10, 'GroupOrder', labels_types); box off; ylim([-5 5]);
yticks([-5 0 5])
subplot(1,3,2)
labels_yeo_adj = strrep(labels_yeo, ' ', '');
min_sig = min(seed_effect(logical(seed_fdr))); %add signiicance threshold to the plot
plot([0 7],[min_sig min_sig],'Color','black');hold on
for net = [1,2,3,4,6,7]%exclude limbic as limbic contains mostly missing parcels /low SNR
    fc_plot_this.(labels_yeo_adj{net}) = seed_effect(yeo_included==net);
end
violinplot(fc_plot_this,fliplr(labels_yeo_adj([1:4,6:7])),'GroupOrder',fliplr(labels_yeo_adj([1:4,6:7])),'ViolinColor',flipud(yeo_cmap([1:4,6:7],:)),'MarkerSize',10); box off; ylim([0 5]); %'MedianColor',[0 0 0],
yticks([0,2.5,5])

% MT change gradient links to general (unthresholded) resilience patterm
subplot(1,3,3)
norm_grad_r = corr(mt_gradient,mt_res_effect');
[norm_grad_p, norm_grad_d]   = spin_test(mt_gradient, mt_res_effect, 'surface_name', 'fsa5', ...
    'parcellation_name', 'glasser_360', 'n_rot', 10000, ...
    'type', 'pearson');

scatter(mt_gradient,mt_res_effect,10,RR3,'o','filled')
hold on
p1 = polyfit(mt_gradient,mt_res_effect, 1); % Linear trendline
f1 = polyval(p1, mt_gradient, 1);
plot(mt_gradient, f1, 'Color','black', 'LineWidth', 2)
%xlabel('\Delta MT gradient','FontName','Arial','FontSize',10)
ylim([-5 5]); xlim([min(mt_gradient)-0.01 max(mt_gradient)+0.01])
%ylabel('\Delta MT vs \Delta resilience','FontName','Arial','FontSize',10); axis square
xticks([]);yticks([])
exportfigbo(f_violin_decoding,[figDir 'delta_mt_decoding.png'],'png', 12)
exportgraphics(f_violin_decoding,[figDir 'delta_mt_decoding.pdf'],'ContentType','vector')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Supplementary analyses  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PFC sub-clusters --> are effects driven by sub-regions of the PFC? E.g. by regions athat are part of a specific network?

if run_supplementaries == 1
    % fc of DMN and Frontoparietal subparts of the significant clusters
    delta_fc_frontoparietal = squeeze(mean(delta_fc_ctx(intersect(find(sig_mt_330),find(yeo_included==6)),:,subjects_to_use_for_delta),1,'omitnan')); %regions that are in the same yeo network will be nans!
    seed_effectp = nan(330,1);seed_effect_SE = nan(330,1);
    for seed = 1:330
        this_seed = delta_fc_frontoparietal(seed,:)';
        this_tbl = [array2table(this_seed,'VariableNames',{'mri'}) delta_table(subjects_to_use_for_delta,:), array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
        lme = fitlme(this_tbl,'mri ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
        seed_effect_fp(seed) =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
        seed_effectp(seed) =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));
        seed_effect_SE(seed) = lme.Coefficients.SE(strcmp(lme.Coefficients.Name,'delta_resilience'));
    end
    seed_fdr = fdr_bh(seed_effectp',0.05);
    to_plot = seed_effect_fp.*seed_fdr;
    to_plot(intersect(find(sig_mt_330),find(yeo_included==6))) = 6; %show ROI in blue
    to_plot = parcel_330_to_360(to_plot,'zer')
    to_plot(not_included_parcels) = -6;
    f_fp = figure;meike_plot_cortical(parcel_to_surface(to_plot,'glasser_360_conte69'),'surface_name','conte69','cmap_custom',cmap_seeds,'color_range',[-4, 4])
    exportfigbo(f_fp,[figDir 'sub_cluster_frontoparietal.png'],'png', 12)

    delta_fc_dmn = squeeze(mean(delta_fc_ctx(intersect(find(sig_mt_330),find(yeo_included==7)),:,subjects_to_use_for_delta),1,'omitnan')); %regions that are in the same yeo network will be nans!
    for seed = 1:330
        this_seed = delta_fc_dmn(seed,:)';
        this_tbl = [array2table(this_seed,'VariableNames',{'mri'}) delta_table(subjects_to_use_for_delta,:), array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
        lme = fitlme(this_tbl,'mri ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
        seed_effect_dmn(seed) =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
        seed_effectp(seed) =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));
        seed_effect_SE(seed) = lme.Coefficients.SE(strcmp(lme.Coefficients.Name,'delta_resilience'));
    end
    seed_fdr = fdr_bh(seed_effectp',0.05);
    to_plot = seed_effect_dmn.*seed_fdr;
    to_plot(intersect(find(sig_mt_330),find(yeo_included==7))) = 6;
    to_plot = parcel_330_to_360(to_plot,'zer');
    to_plot(not_included_parcels) = -6;
    to_plot = parcel_to_surface(to_plot,'glasser_360_conte69');
    f_dmn = figure;meike_plot_cortical(to_plot,'surface_name','conte69','cmap_custom',cmap_seeds,'color_range',[-4, 4])
    exportfigbo(f_dmn,[figDir 'sub_cluster_defaultmode.png'],'png', 12)

    % run permutations to test for significance
    %     %for each parcel
    %     sig_mt_330_idx = find(sig_mt_330);
    %     for parcel = 1:numel(sig_mt_330_idx)
    %         this_data = squeeze(delta_fc_ctx(sig_mt_330_idx(parcel),:,subjects_to_use_for_delta));
    %         for seed = 1:330
    %             if seed~=sig_mt_330_idx(parcel) %this is the diagonal of nans
    %                 this_seed = this_data(seed,:)';
    %                 this_tbl = [array2table(this_seed,'VariableNames',{'mri'}) delta_table(subjects_to_use_for_delta,:), array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
    %                 lme = fitlme(this_tbl,'mri ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
    %                 seed_effect(parcel,seed) =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
    %                 seed_effectp(parcel,seed) =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));
    %             end
    %         end
    %         seed_fdr(parcel,:) = fdr_bh(seed_effectp(parcel,:)',0.05);
    %         to_plot = seed_effect(parcel,:)'.*seed_fdr(parcel,:)';
    %         to_plot(sig_mt_330_idx(parcel))=-2; %show ROI
    %         figure;plot_cortical(parcel_to_surface(parcel_330_to_360(to_plot,'zer'),'glasser_360_conte69'),'surface_name','conte69','color_range',[-4, 4],'label_text',nmfc{sig_mt_rois(parcel)})
    %end

    %% is the effect of MT cluster functional connectivity a global effect?
    fc_table = [table(mean(delta_mt_fc,1)', 'VariableNames',{'seed_fc'}) delta_table(subjects_to_use_for_delta,:) array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
    lme = fitlme(fc_table,'seed_fc ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
    t_seed_global =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
    p_seed_global =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));

    f_scat = figure;
    scatter(delta_table.delta_resilience(subjects_to_use_for_delta,:), mean(delta_mt_fc,1)',10,'filled','black')
    hold on
    p1 = polyfit(delta_table.delta_resilience(subjects_to_use_for_delta,:), mean(delta_mt_fc,1)', 1); % Linear trendline
    f1 = polyval(p1, delta_table.delta_resilience(subjects_to_use_for_delta,:));
    plot(delta_table.delta_resilience(subjects_to_use_for_delta,:), f1, 'Color',[162 11 39]/255, 'LineWidth', 2)
    xlabel('Delta Resilience')
    ylabel('Cluster vs FC change')
    exportfigbo(f_scat, [figDir 'global_fc_of_cluster.png'],'png', 12)
    % is this not specific to the cluster but actually true for the whole
    % brain? no!
    delta_fc_glob = mean(squeeze(mean(delta_fc,1,'omitnan')),1)';
    fc_table = [table(delta_fc_glob(subjects_to_use_for_delta), 'VariableNames',{'global_fc'}) delta_table(subjects_to_use_for_delta,:) array2table(mean_resilience_score(subjects_to_use_for_delta),'VariableName',{'mean_resilience_score'})];
    lme = fitlme(fc_table,'global_fc ~ delta_resilience + mean_age + mean_resilience_score + sex + site'); %0 = female & 0 = resilient, 1 = male & 1 = vulnerable
    t_global =  lme.Coefficients.tStat(strcmp(lme.Coefficients.Name,'delta_resilience'));
    p_global =  lme.Coefficients.pValue(strcmp(lme.Coefficients.Name,'delta_resilience'));
end