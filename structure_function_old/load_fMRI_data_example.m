%%  load data

% Dependencies used:
% colormap: cbrewer: https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
% brainspace: download here: https://brainspace.readthedocs.io/en/latest/pages/install.html
% surfstat: https://mica-mni.github.io/surfstat/
% scientific colourmaps: https://www.fabiocrameri.ch/colourmaps/
% raincloudplots: https://github.com/RainCloudPlots/RainCloudPlots


% Prepare the MPC matrices according to 
% https://github.com/MICA-MNI/micaopen/tree/master/MPC

% Prepare the resting-state matrices like this:
% dummy code of my  approach: 

for load_fc = 1
    HCP400_fc1 = zeros(length(ID),1200,400);
    HCP400_fc2 = HCP400_fc1;
    HCP400_fc3 = HCP400_fc1;
    HCP400_fc4 = HCP400_fc1;
    
    matric_HCP400 = zeros(1206,400,400);
    for i = 1:length(ID)
        i
        try
            mycifti1 = ft_read_cifti([ID '*rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'])
            mycifti2 = ft_read_cifti([ID '*rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'])
            mycifti3 = ft_read_cifti([ID '*rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'])
            mycifti4 = ft_read_cifti([ID '*rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'])
            
            ts1  = mycifti1.dtseries(mycifti1.brainstructure <= 2,:);
            ts2  = mycifti2.dtseries(mycifti2.brainstructure <= 2,:);
            ts3  = mycifti3.dtseries(mycifti3.brainstructure <= 2,:);
            ts4  = mycifti4.dtseries(mycifti4.brainstructure <= 2,:);
            
            
            HCP400_fc1(i,:,:) = labelmean(ts1',HCP400.cdata,'ignorewarning');
            HCP400_fc2(i,:,:) = labelmean(ts2',HCP400.cdata,'ignorewarning');
            HCP400_fc3(i,:,:) = labelmean(ts3',HCP400.cdata,'ignorewarning');
            HCP400_fc4(i,:,:) = labelmean(ts4',HCP400.cdata,'ignorewarning');
            
            fc400(i,:,:) = (corr(squeeze(HCP400_fc1(i,:,:))) + corr(squeeze(HCP400_fc2(i,:,:))) + corr(squeeze(HCP400_fc3(i,:,:))) + ...
                corr(squeeze(HCP400_fc4(i,:,:))))./4;
            %Fischer z transform
            fc400z(i,:,:) = 0.5*log((1+ fc400(i,:,:))./(1- fc400(i,:,:)));
        catch
        end
    end
end

% define the IDs to keep
keep_mpc = find(squeeze(mean(MPC(:,1,1:400),3))>0);
keep_fc = find(squeeze(mean(fc400z(:,1,1:400),3))>0);

keep = intersect(keep_mpc,keep_fc);

hr_fc(eye(400)==1)=0; % zeros in the diagonal
hr_mpc(eye(400)==1)=0; % zeros in the diagonal

% make heritability map:
% use the pedigree structure of HCP in combination with solar-eclipse
% scripts: https://www.nitrc.org/projects/se_linux  
        
        