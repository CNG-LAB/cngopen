
%% Made by Sofie Louise Valk between 2014 and 2021
%% s.valk@fz-juelich.de / valk@cbs.mpg.de

% Dependencies used:
% brainspace: download here: https://brainspace.readthedocs.io/en/latest/pages/install.html
% surfstat: https://mica-mni.github.io/surfstat/
% scientific colourmaps: https://www.fabiocrameri.ch/colourmaps/
% raincloudplots: https://github.com/RainCloudPlots/RainCloudPlots
% canlab-core: https://github.com/canlab/CanlabCore for prediction

% load brain surface
% SN = brain surface to project results on (fsaverage 5)
% mask = no regions in mid-surface cut
% keepallrest = %movement data + bad quality

% load the gradients
% [closed_data] unfortunately data cannot be shared, please contact Sofie Valk 
% for details.
%
% Basically, vertex-level functional connectivity matrices are created for each
% individual at each time-point, using functional data that is projected on
% fsaverage5 surfaces. Following, gradients are constructed using
% brainspace (https://brainspace.readthedocs.io/en/latest/) and aligned
% using procrustes alignment to a template functional connectivity matrix
% based on the HCP S1200 young adult sample (van Essen, 2013).
% Following the next analysis are performed, similar to this work: 
% https://elifesciences.org/articles/64694. We reran analysis using a
% parcellation framework (e.g. summarized the data within the Schaefer 400
% parcels (Schaefer, 2018) which gave comparable results- but made
% computations faster- and also reran analysis using the BCT toolbox to
% generate alternative measures of integration and segregation: clustering coefficient
% and path length (https://sites.google.com/site/bctnet/).

% vertex level eccentricity for participants across timepoints
for i = 1:992
    GGG(i,:) = sqrt((G1(i,:).^2)+(G2(i,:).^2)+(G3(i,:).^2));
end

% make difference score for eccentricity
for difference_score = 1
        keep1 = find(~isnan(GGG(:,2))); % remove those without scanning data
        keep3 = mintersect(keep1, find(keepallrest2==1)); % remove outliers with 0.3 mm/degree movement
        keep_bv = keep3;
        keeptest = zeros(length(id),1);
        keeptest(keep3) = 1;
   
        eZ = zeros(992,20484);
        eZ(:,mask==1) = GGG(1:992,mask==1);
        eZ(isnan(eZ)) = 0;
        %remove outliers based on correlation with mean eccentrcity (motivated by Hong, 2020)
        r = corr(mean(eZ)',eZ');     
        outliers_ggg = (r>0.5)';
       
        keep = find(keeptest==1);
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        eChange_bl  = zeros(size(eZ))-666; % -666 operation is chosen to exclude 0 data otherwise unnoticed by the substraction procedure 
        eChange_t1  = zeros(size(eZ))-666;
        eChange_t2  = zeros(size(eZ))-666;
        
        %     who have data there?
        for i = 1:length(id) %id is participant id per timepoint
            if keeptest(i,:).*outliers_ggg(i,:) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeptest(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeptest(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeptest(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeptest(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeptest(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeptest(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeptest(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeptest(inter3);
                end
                if keeptest(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    eChange_bl(i,:,:)    = (eZ(i,:,:) - eZ(inter0,:,:));
                end
                if keeptest(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    eChange_t1(i,:,:)    = (eZ(i,:,:) - eZ(inter1,:,:));
                end
                if keeptest(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    eChange_t2(i,:,:)    = (eZ(i,:,:) - eZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(GGG));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = eChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = eChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = eChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        Z_last = ZChange_last;
end


% make difference score for G1
for difference_score = 1
        keep1 = find(~isnan(GGG(:,2)));
        keep3 = mintersect(keep1, find(keepallrest2==1));
        keep_bv = keep3;
        keeptest = zeros(length(id),1);
        keeptest(keep3) = 1;
   
        eZ = zeros(992,20484);
        eZ(:,mask==1) = G1all(1:992,mask==1);
        eZ(isnan(eZ)) = 0;
        %remove outliers
        r = corr(mean(eZ)',eZ');     
        outliers_g1 = (r>0.5)';
   
        keep = find(keeptest==1);
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        eChange_bl  = zeros(size(eZ))-666;
        eChange_t1  = zeros(size(eZ))-666;
        eChange_t2  = zeros(size(eZ))-666;
        
        %     who have data there?
        for i = 1:length(id)
            if keeptest(i,:).*outliers_g1(i,:) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeptest(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeptest(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeptest(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeptest(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeptest(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeptest(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeptest(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeptest(inter3);
                end
                if keeptest(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    eChange_bl(i,:,:)    = (eZ(i,:,:) - eZ(inter0,:,:));
                end
                if keeptest(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    eChange_t1(i,:,:)    = (eZ(i,:,:) - eZ(inter1,:,:));
                end
                if keeptest(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    eChange_t2(i,:,:)    = (eZ(i,:,:) - eZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(GGG));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = eChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = eChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = eChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        G1_last = ZChange_last;
end

% make difference score for G2
for difference_score = 1
        keep1 = find(~isnan(GGG(:,2)));
        keep3 = mintersect(keep1, find(keepallrest2==1));
        keep_bv = keep3;
        keeptest = zeros(length(id),1);
        keeptest(keep3) = 1;
   
        eZ = zeros(992,20484);
        eZ(:,mask==1) = G2all(1:992,mask==1);
        eZ(isnan(eZ)) = 0;
        r = corr(mean(eZ)',eZ');     
        outliers_g2 = (r>0.5)';
   
       
        keep = find(keeptest==1);
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        eChange_bl  = zeros(size(eZ))-666;
        eChange_t1  = zeros(size(eZ))-666;
        eChange_t2  = zeros(size(eZ))-666;
        
        %     who have data there?
        for i = 1:length(id)
            if keeptest(i,:).*outliers_g2(i,:) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeptest(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeptest(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeptest(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeptest(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeptest(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeptest(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeptest(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeptest(inter3);
                end
                if keeptest(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    eChange_bl(i,:,:)    = (eZ(i,:,:) - eZ(inter0,:,:));
                end
                if keeptest(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    eChange_t1(i,:,:)    = (eZ(i,:,:) - eZ(inter1,:,:));
                end
                if keeptest(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    eChange_t2(i,:,:)    = (eZ(i,:,:) - eZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(GGG));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = eChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = eChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = eChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        G2_last = ZChange_last;
end

% make difference score for G3
for difference_score = 1
        keep1 = find(~isnan(GGG(:,2)));
        keep3 = mintersect(keep1, find(keepallrest2==1));
        keep_bv = keep3;
        keeptest = zeros(length(id),1);
        keeptest(keep3) = 1;
   
        eZ = zeros(992,20484);
        eZ(:,mask==1) = G3all(1:992,mask==1);
        eZ(isnan(eZ)) = 0;
        r = corr(mean(eZ)',eZ');     
        outliers_g3 = (r>0.5)';
   

        keep = find(keeptest==1);
        minit = zeros(size(id));
        tp_t0_there = minit;
        tp_t1_there = minit;
        tp_t2_there = minit;
        tp_t3_there = minit;
        daysfrombl  = minit;
        daysfromt1  = minit;
        daysfromt2  = minit;
        eChange_bl  = zeros(size(eZ))-666;
        eChange_t1  = zeros(size(eZ))-666;
        eChange_t2  = zeros(size(eZ))-666;
        
        %     who have data there?
        for i = 1:length(id)
            if keeptest(i,:).*outliers_g3(i,:) == 1
                i
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* keeptest(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* keeptest(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* keeptest(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* keeptest(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* keeptest(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* keeptest(inter2);
                end
                tp_t3_there(i,1)    = ~isempty(inter3) .* keeptest(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* keeptest(inter3);
                end
                if keeptest(inter0)~=0
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    eChange_bl(i,:,:)    = (eZ(i,:,:) - eZ(inter0,:,:));
                end
                if keeptest(inter1)~=0
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    eChange_t1(i,:,:)    = (eZ(i,:,:) - eZ(inter1,:,:));
                end
                if keeptest(inter2)~=0
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    eChange_t2(i,:,:)    = (eZ(i,:,:) - eZ(inter2,:,:));
                end
            else
            end
        end
        
        ZChange_last = zeros(size(GGG));
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:,:)  = eChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:)  = eChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:)  = eChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:,:);
        
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:) = eChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:) = eChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:,:);
        ZChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:) = eChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:,:);
        
        
        G3_last = ZChange_last;
end
