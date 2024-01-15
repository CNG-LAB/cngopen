function [MI, p_MI, ratio] = compute_MI_lme(data, demographics_tbl)

MI = nan(size(data,1),1);
p_MI = nan(size(data,1),1);

for roi = 1:size(data,1)
    %tic
    % Extract data for this ROI
    this_roi = squeeze(data(roi,:,:))';
    this_roi(:,roi)=[]; %remove nans to fit GLM (diagonal in FC data)
    missing = find(sum(this_roi,1)==0);
    this_roi(:,missing)=[];
    %roi %track progress
    parfor sroi = 1:size(this_roi,2)
        really_this_roi = this_roi(:,sroi);
        tbl = [array2table(really_this_roi,'VariableNames',{'mri'}) demographics_tbl];
        lme = fitlme(tbl,'mri ~ age  + sex + site + (1|subj)');

        change(sroi,roi) = lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'));
        % is male or female coded? depends which one comes first
        % site can be site1 site2 site3, use site_index

       site_idx =  find(contains(lme.Coefficients.Name,'site'));
        if sum(strcmp(lme.Coefficients.Name, 'sex_female'))>0
            baseline(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_female'))*(1/2) + lme.Coefficients.Estimate(site_idx(1))*(1/3) + lme.Coefficients.Estimate(site_idx(2))*(1/3);
        else
            baseline(sroi,roi) = lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'age'))*14 + lme.Coefficients.Estimate(strcmp(lme.Coefficients.Name,'sex_male'))*(1/2) + lme.Coefficients.Estimate(site_idx(1))*(1/3) + lme.Coefficients.Estimate(site_idx(2))*(1/3);
        end
    end
        [MI(roi,:) ,p_MI(roi,:)] = corr(baseline(:,roi), change(:,roi), 'type', 'Spearman', 'rows', 'complete');
        ratio(roi,:) = sum(change(:,roi)>0) / size(this_roi,2); 
%toc
end