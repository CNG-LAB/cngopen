%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coordinated Cortical Thickness alterations Across Mental Disorders: A Transdiagnostic ENIGMA Project
% Script 1: Loads ENIGMA summary statistics (Cohen's d maps), computes
% cross-disorder similarity, derives transdiagnostic Gradients, and
% contextualizes these gradients with functional, macro- and
% microstructural and transcriptomic profiles.

%This script runs analyses presented in Figure 2 of the manuscript.
%Meike Hettwer 2022 - CNG Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add Toolboxes and data paths

%set your paths
%DataPath = '....';
%ScriptPath = '....';

%% add Toolboxes and data paths
%required toolboxes are:
%ENIGMA-1.1.3
%BrainSpace-0.1.2
%SurfStat
%cbrewer
%spiderPlot
%gifti-master
%quantileranks
%GAMBA

%load data from GitHub
%CT_psych_gradients.mat
%Valk2020_CT_structural_covariance_gradient.mat
%midbrain_vertices_fsa5.mat     %useful for conversion between atlases (note: extra midbrain - parcels when mapping from fsa5 to Desikan-Killiany are: 1 5 40)
%DK_midbrain_parcels.mat

%% 1. Load Summary statistics Data
%If data uploaded on Github is used:
load(strcat(DataPath,'CT_6_Disorders.mat'));

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
%% Disorders: BD, SZC, ADHD, ASD, MDD, OCD
disorders = [CT_d_bipo, CT_d_adhd, CT_d_asd, CT_d_mdd, CT_d_ocd, CT_d_sz];
disorder_names = {'Bipolar Disorder'; 'ADHD'; 'ASD'; 'Depression'; 'OCD'; 'Schizophrenia'};
%% Compute Transdiagnostic Gradients (Figure 2)

% 1) derive cross-disorder correlation matrix
corr_dis = corr(disorders');

% 2) build/load gradients
gm = GradientMaps('kernel','na','approach','dm'); %in GradientMaps Script, set "sparsity" to 80 (due to low number of DK parcels)
gm = gm.fit(corr_dis);
psych_grad = gm.gradients{1};

%Figure 2A
handles = scree_plot(gm.lambda{1});set(handles.axes,'FontName','Gill Sans MT','FontSize',12,'XTick',[1 7],'YTick',[0 0.4],'YTickLabel',{'0','0.4'});
CT_psych_gradient1 = psych_grad(:,1);
CT_psych_gradient2 = psych_grad(:,2);

%alternatively, load pre-computed gradients from GitHub
%load('CT_psych_gradients');

%plot gradients to surface (without midbrain values) (*Figure 2B*)
plot_G1_no_midbrain = parcel_to_surface(CT_psych_gradient1);
plot_G1_no_midbrain(midbrain_verticesfsa5) = (max(plot_G1_no_midbrain)+min(plot_G1_no_midbrain))/2;
plot_G2_no_midbrain = parcel_to_surface(CT_psych_gradient2);
plot_G2_no_midbrain(midbrain_verticesfsa5) = (max(plot_G2_no_midbrain)+min(plot_G2_no_midbrain))/2;

figure, plot_cortical(plot_G1_no_midbrain);
figure, plot_cortical(plot_G2_no_midbrain); 
%% Compare Gradients to Valk 2020 (DOI: 10.1126/sciadv.abb3417) normative structural covariance gradient of cortical thickness
%Figure 2C
%re-parcel both normative CT covariance gradients to DK space (via fsa5 surface)
Valk_G1_fsa5 = parcel_to_surface(Valk2020_CT_structural_covariance_gradient(:,1),'schaefer_400_fsa5');
Valk_G2_fsa5 = parcel_to_surface(Valk2020_CT_structural_covariance_gradient(:,2),'schaefer_400_fsa5');
Valk_G1_DK = surface_to_parcel(Valk_G1_fsa5)';
Valk_G2_DK = surface_to_parcel(Valk_G2_fsa5)';
Valk_G1_DK(midbrain_parcels) = [];%remove midbrain parcels for spin test
Valk_G2_DK(midbrain_parcels) = [];
%test association with principal CT gradient via spin test
r_G1_CT = corr(Valk_G1_DK,CT_psych_gradient1);
r_G2_CT = corr(Valk_G1_DK,CT_psych_gradient2);
[p_spin_G1] = spin_test(Valk_G1_DK,CT_psych_gradient1, 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');
[p_spin_G2] = spin_test(Valk_G1_DK,CT_psych_gradient2, 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');
%test association with second CT gradient via spin test
r_G1_CT2 = corr(Valk_G2_DK,CT_psych_gradient1);
r_G2_CT2 = corr(Valk_G2_DK,CT_psych_gradient2);
[p_spin_G1_2] = spin_test(Valk_G2_DK,CT_psych_gradient1, 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');
[p_spin_G2_2] = spin_test(Valk_G2_DK,CT_psych_gradient2, 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');

%for visualization, color in midbrain with intermediate color (white)
Valk_plot_G1_no_midbrain = parcel_to_surface(Valk_G1_DK);
Valk_plot_G1_no_midbrain(midbrain_verticesfsa5)=(max(Valk_plot_G1_no_midbrain)+min(Valk_plot_G1_no_midbrain))/2;
figure; plot_cortical(Valk_plot_G1_no_midbrain)

%plot correlation between G1/G2 and CT gradient (*Figure 2C*)
Trans_Gradients = [CT_psych_gradient1 CT_psych_gradient2];
col = [0 0 0];
figure;
for j = 1:2
    axs = subplot(1, 2, j); hold on
    s   = scatter(Valk_G1_DK, Trans_Gradients(:,j), 38, col, 'filled');
    P1      = polyfit(Valk_G1_DK, Trans_Gradients(:,j), 1);
    yfit_1  = P1(1) * Valk_G1_DK + P1(2);
    hold on
    col_line = [0.66 0.13 0.11];%[53.3/255 0 0]
    plot(Valk_G1_DK, yfit_1, 'color',col_line, 'LineWidth', 3);
    ylim([min(Trans_Gradients(:,j))-0.01 max(Trans_Gradients(:,j))+0.01]);
    xlim([min(Valk_G1_DK)-0.01 max(Valk_G1_DK)+0.01]);
    set(gca,'xtick',[],'ytick',[]);
    axis square
end

%% Contextualize Gradients with cytoarchitectonic classes (von Economo Koskinas; see https://www.karger.com/Article/Abstract/103258 or https://enigma-toolbox.readthedocs.io/en/latest/pages/11.02.voneconomo/index.html)
%Figure 2D
%stratify gradients according to cytoarchitectonic classes 
%adapted 2025
ve = dlmread('economo_koskinas_fsa5.csv'); %original order: {'Agranular', 'Frontal', 'Parietal', 'Polar', 'Granular'}

% parcellate von economo classes to DK
label_vector     = dlmread(['aparc_fsa5' '.csv']); 
uparcel          = unique(label_vector); 
ve_dk      = zeros(length(uparcel),1); 

for ii = 1:length(uparcel) 
    thisparcel          = uparcel(ii); 
    ve_dk(ii)           = mode(ve(label_vector == thisparcel));
end
% remove midbrain
ve_dk(midbrain_parcels)= [];

% average gradient loadings per cyto class

class_meanG1 = nan(1,5);
class_meanG2 = nan(1,5);

for ii = 1:5
    id = find(ve_dk == ii);
    class_meanG1(ii) = mean(CT_psych_gradient1(id));
    class_meanG2(ii) = mean(CT_psych_gradient2(id));
end
classes_spider = [class_meanG1; class_meanG2];

%Reorder for better visualization
classes_spider_reordered(:,1:4) = classes_spider(:,2:5);
classes_spider_reordered(:,5) = classes_spider(:,1);
vek_reordered = {'Frontal','Parietal', 'Polar', 'Granular','Agranular'}; 
AxesLim = [- 0.1 -0.1; 0.15 0.15];

figure;
spider_plot(classes_spider_reordered, 'AxesLabels', vek_reordered,'AxesLabelsEdge', 'none','AxesLimits', [-0.26 -0.26 -0.26 -0.26 -0.26; 0.26 0.26 0.26 0.26 0.26],'AxesPrecision',2,...
    'AxesDisplay','one','AxesInterval',5,  'Color', [124, 30, 43;107 142 35; 1 0 1 ]/255,...
    'AxesFont', 'Gill Sans MT','LabelFont', 'Gill Sans MT','AxesFontSize', 12, 'LabelFontSize', 15);
legend_str = {'G1', 'G2'};legend(legend_str, 'Location', 'northeast'); legend boxoff

%% Identify genes whose expression pattern follows G1 / G2. 
% Identified genes were then fed into the CSEA tool for developmental enrichment analyses (http://genetics.wustl.edu/jdlab/csea-tool-2/)

%% Load AHBA data and test association with transdiagnostic gradients using null models tsting for spatial and gene specificity. Requires GAMBA Toolbox (Wei et al.)

%AHBA
ahba = fetch_ahba();
regions = ahba.label;
gene_symbols = ahba.Properties.VariableNames(2:end)';
expression = table2array(ahba(1:68,2:end));

% 1. assess spatial specificity via alexander-bloch spin test

for i = 1:length(expression)
    [p_spin_G1(i)] = spin_test_spin_map1_only(CT_psych_gradient1, expression(:,i), 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');
    [p_spin_G2(i)] = spin_test_spin_map1_only(CT_psych_gradient2, expression(:,i), 'surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000, 'type', 'pearson');
end
p_G1_sig = p_spin_G1(p_spin_G1<0.01);
p_G2_sig = p_spin_G1(p_spin_G2<0.01);

G1_gene_set_01 = gene_symbols(p_spin_G1<0.01);
G2_gene_set_01 = gene_symbols(p_spin_G2<0.01);

[r_G1, p_G1] = corr(CT_psych_gradient1,expression,'rows','complete'); %1:68 = only cortical parcellations (69:82 = subcortical); 2:end = all genes
[r_G2, p_G2] = corr(CT_psych_gradient2,expression,'rows','complete'); %1:68 only cortical parcellations (69:82 = subcortical); 2:end = all genes

r_G1_sig = r_G1(p_spin_G1<0.01);
r_G1_sig_positive = r_G1_sig>0;
G1_positive_genes = G1_gene_set_01(r_G1_sig_positive);
G1_negative_genes = G1_gene_set_01(~r_G1_sig_positive);

T_genes_G1 = table(G1_gene_set_01, r_G1_sig', p_G1_sig'); 
T_genes_G1.Properties.VariableNames = {'Gene','r','p_spin'};
writetable(T_genes_G1,'Genes_G1.xlsx','Sheet',1,'Range','A1')

% 2. assess gene specificity using GAMBA
%check whether previously identified gene set is significantly stronger
%associated with phenotype than random genes (that are also over-expressed
%in the brain!)

%res_G1_null_brain = permutation_null_brain(CT_psych_gradient1, G1_gene_set_05, expression, gene_symbols, 'brain')
res_G1_null_brain_p01 = permutation_null_brain(CT_psych_gradient1, G1_gene_set_01, expression, gene_symbols, 'brain');
res_G2_null_brain_p01 = permutation_null_brain(CT_psych_gradient2, G2_gene_set_01, expression, gene_symbols, 'brain');

% check whether imaging profiles associate with gene expression profiles, based
% on the null-coexpression model (where random genes with similar
% coexpression level is conserved).
res_G1_null_coexp_p01 = permutation_null_coexp(CT_psych_gradient1, G1_gene_set_01, expression, gene_symbols);
res_G2_null_coexp_p01 = permutation_null_coexp(CT_psych_gradient2, G2_gene_set_01, expression, gene_symbols);

%Gene set is then fed in to the CSEA Tool and presented in *Figure 2E*
%% Neurosynth functional decoding (Figure 2E)
%Functional Neurosynth decoding

cd /Users/m.hettwer/Documents/MATLAB/ENIGMA_Cross_Disorder_Gradients/Data/neurosynth_z_values;
BH =     {'action_association-test_z_lh.shape.gii'
    'action_association-test_z_rh.shape.gii'
    'affective_association-test_z_lh.shape.gii'
    'affective_association-test_z_rh.shape.gii'
    'attention_association-test_z_lh.shape.gii'
    'attention_association-test_z_rh.shape.gii'
    'auditory_association-test_z_lh.shape.gii'
    'auditory_association-test_z_rh.shape.gii'
    'autobiographical_memory_association-test_z_lh.shape.gii'
    'autobiographical_memory_association-test_z_rh.shape.gii'
    'cognitive_control_association-test_z_lh.shape.gii'
    'cognitive_control_association-test_z_rh.shape.gii'
    'emotion_association-test_z_lh.shape.gii'
    'emotion_association-test_z_rh.shape.gii'
    'episodic_memory_association-test_z_lh.shape.gii'
    'episodic_memory_association-test_z_rh.shape.gii'
    'eye_movement_association-test_z_lh.shape.gii'
    'eye_movement_association-test_z_rh.shape.gii'
    'face_association-test_z_lh.shape.gii'
    'face_association-test_z_rh.shape.gii'
    'inhibition_association-test_z_lh.shape.gii'
    'inhibition_association-test_z_rh.shape.gii'
    'language_association-test_z_lh.shape.gii'
    'language_association-test_z_rh.shape.gii'
    'motor_association-test_z_lh.shape.gii'
    'motor_association-test_z_rh.shape.gii'
    'multisensory_association-test_z_lh.shape.gii'
    'multisensory_association-test_z_rh.shape.gii'
    'pain_association-test_z_lh.shape.gii'
    'pain_association-test_z_rh.shape.gii'
    'reading_association-test_z_lh.shape.gii'
    'reading_association-test_z_rh.shape.gii'
    'reward_association-test_z_lh.shape.gii'
    'reward_association-test_z_rh.shape.gii'
    'semantics_association-test_z_lh.shape.gii'
    'semantics_association-test_z_rh.shape.gii'
    'social_cognition_association-test_z_lh.shape.gii'
    'social_cognition_association-test_z_rh.shape.gii'
    'verbal_association-test_z_lh.shape.gii'
    'verbal_association-test_z_rh.shape.gii'
    'visual_association-test_z_lh.shape.gii'
    'visual_association-test_z_rh.shape.gii'
    'visual_perception_association-test_z_lh.shape.gii'
    'visual_perception_association-test_z_rh.shape.gii'
    'visuospatial_association-test_z_lh.shape.gii'
    'visuospatial_association-test_z_rh.shape.gii'
    'working_memory_association-test_z_lh.shape.gii'
    'working_memory_association-test_z_rh.shape.gii'};
LH_ind = find(endsWith(BH, 'lh.shape.gii'));
LH = BH(LH_ind);
RH_ind = find(endsWith(BH, 'rh.shape.gii'));
RH = BH(RH_ind);
%%reparcel gradients and neurosynth data to same space (71 parcels), since
%%our Conte69 to DK converision csv file describes the 71 (instead of the
%%68) parcel version
Grad1_surf = parcel_to_surface(CT_psych_gradient1,'aparc_conte69');
Grad1 = surface_to_parcel(Grad1_surf, 'aparc_conte69');
Grad_surf2 = parcel_to_surface(CT_psych_gradient2,'aparc_conte69');
Grad2 = surface_to_parcel(Grad_surf2, 'aparc_conte69');

surf =load_conte69('surfaces');

%% bin neurosynth data
tmp = nan(24,20);
for j = 1:24
    L = gifti(strcat('/Users/m.hettwer/Documents/MATLAB/ENIGMA_Cross_Disorder_Gradients/Data/neurosynth_z_values/', LH{j}));
    R = gifti(strcat('/Users/m.hettwer/Documents/MATLAB/ENIGMA_Cross_Disorder_Gradients/Data/neurosynth_z_values/', RH{j}));
    funcmap= [L.cdata;R.cdata];
    %bins1   = [];
    bins1   = quantileranks(Grad1,20); %Grad1 binned into 1:20 bins and allocated to 71 DK parcels
    grad_fc = zeros(1,64984);
    for i = 1:71
        a = i-1; %because parcels are numbered from 0 to 70 and not 1 to 71 (due to previous Python indexing); bins from 1:20
        grad_fc(:,aparc_conte69==a) = bins1(i);
    end
    for i = 1:20
        tmp(j,i) = mean(funcmap(grad_fc==i)); %size: 24 x 20
    end
end

%zscore binned functional map
Grad1_binned = nan(1,20);
for i = 1:20
    Grad1_binned(i) = mean(Grad1(bins1==i));
end
tmpz2= zscore(tmp,[],2);
tmpz22 = tmpz2;
tmpz2(tmpz2<0.5) = nan;
weighted_mean = nanmean(tmpz2.*(Grad1_binned),2);
[~,b,c] = unique(weighted_mean);

tmp2 = nan(24,20);
grd_binned = nan(1,20);
for j = 1:24
    L = gifti(strcat('/Users/m.hettwer/Documents/MATLAB/ENIGMA_Cross_Disorder_Gradients/Data/neurosynth_z_values/', LH{j}));
    R = gifti(strcat('/Users/m.hettwer/Documents/MATLAB/ENIGMA_Cross_Disorder_Gradients/Data/neurosynth_z_values/', RH{j}));
    funcmap= [L.cdata;R.cdata];
%     bins1   = [];
    bins1   = quantileranks(Grad2,20);
    grad_fc = zeros(1,64984);
    for i = 1:71
        a = i-1;
        grad_fc(:,find(aparc_conte69==a)) = bins1(i);
    end
    for i = 1:20
        tmp2(j,i) = mean(funcmap(find(grad_fc==i)));
    end
end
for i = 1:20
    grd_binned(i) = mean(Grad2(bins1==i));
end
tmpz= zscore(tmp2,[],2);
tmpzz = tmpz;
tmpz(tmpz<0.5) = nan;
weighted_mean = nanmean(tmpz.*(grd_binned),2);
[~,d,e] = unique(weighted_mean);

j=0;
for i = 1:24
    j=j+1;
    words{j} = LH{j};
    worrs3{j} = words{j}(1:end-32);
    words_final{j} = strrep(worrs3{j},'_',' ');
    wrd{j} = num2str(j);
end
eucd= sqrt(c.^2+d.^2);

%% color coding
% col = [0 0 0];
BBgrad1 = CT_psych_gradient1;
BBgrad2 = CT_psych_gradient2;

g = load('/Users/m.hettwer/Documents/MATLAB/Toolbox/cbrewer/cbrewer/cbrewer/colorbrewer.mat');
colorbrewer = g.colorbrewer;

G1peak = max(BBgrad1);
G2min = min(BBgrad2);
G2max = max(BBgrad2);

colourness = nan(68,3);
for ii = 1:length(BBgrad1)
    
    if BBgrad1(ii) > 0
        colourness(ii,2) = 1 - (G1peak - BBgrad1(ii));
    else
        colourness(ii,2) = 1 - (G1peak + abs(BBgrad1(ii)));
    end
    if BBgrad2(ii) > 0
        colourness(ii,1) = 1 - (G2max- BBgrad2(ii));
        colourness(ii,3) = 1 - (abs(G2min) + BBgrad2(ii));
    else
        colourness(ii,1) = 1 - (G2max + abs(BBgrad2(ii)));
        colourness(ii,3) = 1 - (abs(G2min) - abs(BBgrad2(ii)));
    end
    
end
% rescale colorbar 
colourness_rescale = nan(68,3);
colourness_vertices = nan(64984,3);
for col = 1:3
    colourness_rescale(:,col) = rescale(colourness(:,col), 0, 1);
    colourness_vertices(:,col)=parcel_to_surface(colourness_rescale(:,col),'aparc_conte69');
end

%compute densities for density plot next to scatters
[x, y] = ksdensity(BBgrad1);
y2(1,:) = x; x2(1,:) = y;
[y2(2,:), x2(2,:)] = ksdensity(BBgrad2);

%% scatter (*Figure 2F*)
colourness_rescale_green =     colourness_rescale;
colourness_rescale_green(:,2) = colourness_rescale_green(:,2)-min(colourness_rescale(:,2));

f= figure;
a(1) = axes('position', [0.1 0.1 0.4 0.4]); %plot colored parcels along gradients
scatter(BBgrad2, BBgrad1, 60, colourness_rescale_green,'filled');
axis([min(BBgrad2)-0.01 max(BBgrad2)+0.01 min(BBgrad1)-0.01 max(BBgrad1)+0.01])
set(gca,'xtick',[],'ytick',[])

a(2) = axes('position', [0.5 0.1 0.07 0.4]); %add densities
patch(y2(1,:),x2(1,:), ones(1,length(y2(1,:)))); axis off;
colormap(a(2), [190, 190, 190]/250)%[0.5, 0.5, 0.5])
ylim([min(BBgrad1) max(BBgrad1)])

a(3) = axes('position', [0.1 0.5 0.4 0.07]);
patch(x2(2,:), y2(2,:), ones(1,length(y2(2,:)))); axis off;
colormap(a(3), [190, 190, 190]/250)
xlim([min(BBgrad2) max(BBgrad2)])

hold on
a(4) = axes('position', [0.1 0.1 0.4 0.4]); col = [0,0,0]; %add cognitive function terms
x = e; y = c; scatter(x,y,60,col,'filled');%flipped, so G1 is on y axis
txty = words_final;
dx = 0.85;
dy = 0.2; % displacement so the text does not overlay the data points
text(x+dx, y+dy, txty,'FontSize',13, 'FontName', 'Gill Sans MT');
set(gca,'XTick',[],'YTick',[],'color','none');xlabel('G2');ylabel('G1')
hold off

%% END OF SCRIPT





