function [MI, p_MI] = compute_MI_lme_interaction(data, demographics_tbl, contrast)

for roi = 1:size(data,1)
    % Extract data for this ROI
    this_roi = squeeze(data(roi,:,:))';
    this_roi(:,roi)=[]; %remove nans to fit GLM (diagonal in FC data)
    missing = find(sum(this_roi,1)==0);
    this_roi(:,missing)=[];
    roi %track progress
    parfor sroi = 1:size(this_roi,2)
        really_this_roi = this_roi(:,sroi);
        tbl = [array2table(really_this_roi,'VariableNames',{'mri'}) demographics_tbl];
        lme = fitlme(tbl,'mri ~ Resilience + age + Resilience*age + sex + site  + (1|subj)');

        if sum(strcmp(lme.Coefficients.Name, 'Resilience_resilient'))>0
            contrast = 'Resilience_resilient';
        else
            contrast = 'Resilience_vulnerable';
        end

        change_1(sroi,roi) = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age')) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*1;
        change_0(sroi,roi) = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age')) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*0;
        if sum(strcmp(lme.Coefficients.Name, 'sex_female'))>0 %which sex is coded first
            baseline_1(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,contrast))*1 ...
                + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_female'))*(1/2) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*1*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'site_site1'))*(1/3) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'site_site2'))*(1/3);
            baseline_0(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,contrast))*0 ...
                + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_female'))*(1/2) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*0*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'site_site1'))*(1/3) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'site_site2'))*(1/3);

        else
            this_site = lme.Coefficients.Name(find(contains(lme.Coefficients.Name,'site')));
            baseline_1(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,contrast))*1 ...
                + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_male'))*(1/2) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*1*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,this_site{1}))*(1/3) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,this_site{2}))*(1/3);
            baseline_0(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,contrast))*0 ...
                + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_male'))*(1/2) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,['age:' contrast]))*0*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,this_site{1}))*(1/3) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,this_site{2}))*(1/3);

        end
    end
    [MI.([contrast '_0'])(roi,:) ,p_MI.([contrast '_0'])(roi,:)] = corr(baseline_0(:,roi), change_0(:,roi), 'type', 'Spearman', 'rows', 'complete');
    [MI.([contrast '_1'])(roi,:) ,p_MI.([contrast '_1'])(roi,:)] = corr(baseline_1(:,roi), change_1(:,roi), 'type', 'Spearman', 'rows', 'complete');
end








