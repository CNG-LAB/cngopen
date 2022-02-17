
% Fig 5. Principal gradient of heritability and individual variability. 
% A) Principal gradient of MPC and rsFC, left: gradients based on mean data 
% on MPC and rsFC and right gradients of heritable data alone, lower left 
% panel: mean versus heritable MPC G1, as well as heritably along the 
% principal mean gradient in MPC; lower right panel: mean versus heritable
% rsFC G1, as well as heritably along the principal mean gradient in rsFC; 
% B) Principal gradient of individual variation (std) in MPC and rsFC, 
% lower panel left: correlation between mean and std MPC G1, and std along
% the mean gradient of MPC; lower panel right: correlation between mean and 
% std rsFC G1, and std along the mean gradient of FC. Source data are 
% provided as a Source Data file.


% Gradient construction
% for more details see also brainspace.readthedocs.io

% gradient in humans and their differences: 

for heritable_gradients_mpc = 1
    % Gradient diffusion map embedding with procrustes alignment
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({MPCm, hr_mpc});
    
    mpc_g1 = gm.aligned{1}
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =mpc_g1(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = mpc_g1(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(batlow)
  
     mpc_h = rescale(gm.aligned{2}(:,1))
     
     heri_ct = zeros(1,20484);
     for i = 1:200
         heri_ct(:,find(parcels400==i+1)) = mpc_h(i);
     end
     for i = 1:200
         heri_ct(:,find(parcels400==i+1001)) = mpc_h(i+200);
     end
     
     f = figure,
     BoSurfStatViewData(-heri_ct,SN,'')
     colormap(bilbao)
     
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1),'filled','k'),lsline
    xlim([-0.15 0.21])
    
    % heritability along principal gradient
    bins1   = [];
    bins1   = quantileranks(-gm.aligned{1}(:,1),10);
  
    for j=1:10
        for i = 1:10
            mean_herit_m(j,i) = mean(mean(hr_mpc(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_herit_m)
    colormap((cbrewer('seq','Reds',99)))
    colorbar
  
end

for heritable_gradients_fc = 1
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({fc400m, hr_fc});
    
    fc_g1 = gm.aligned{1}
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =fc_g1(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = fc_g1(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(batlow)
  
     fc_h = rescale(gm.aligned{2}(:,1))
     
     heri_ct = zeros(1,20484);
     for i = 1:200
         heri_ct(:,find(parcels400==i+1)) =fc_h(i);
     end
     for i = 1:200
         heri_ct(:,find(parcels400==i+1001)) = fc_h(i+200);
     end
     
     f = figure,
     BoSurfStatViewData(-heri_ct,SN,'')
     colormap(bilbao)
     
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1),'filled','k'),lsline
    xlim([-0.15 0.21])
    
    % heritability along principal gradient
    bins1   = [];
    bins1   = quantileranks(-gm.aligned{1}(:,1),10);
  
    for j=1:10
        for i = 1:10
            mean_herit_m(j,i) = mean(mean(hr_fc(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_herit_m)
    colormap((cbrewer('seq','Reds',99)))
    colorbar
end

for std_gradients = 1
    MPCstd = squeeze(std(MPC(keep,:,:)));
    MPCstd(eye(size(MPCstd))==1) = 0;
    
    % gradients pure and std and heritable
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({meanMPC, MPCstd});
    
    mpc_stg1 = gm.aligned{2}
       
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =mpc_stg1(i,1);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = mpc_stg1(i+200,1);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(roma)
    
    bins1   = [];
    bins1   = quantileranks(-mpc_g1(:,1),10);
    
    for j=1:10
        for i = 1:10
            mean_stdmpc_m(j,i) = mean(mean(MPCstd(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_stdmpc_m,[0.1 0.8])
    colormap((cbrewer('seq','Reds',99)))
    colorbar
      
    scatter(gm.aligned{1}(:,1),gm.aligned{2}(:,1));
    
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1), 'filled','k'),lsline
    xlim([-0.15 0.21])
   
    
    fc400std = squeeze(std(fc400z(keep,:,:)));
    fc400std(eye(size(fc400std))==1) = 0;
    
    gm = GradientMaps('kernel','na','approach','dm','align','pa');
    gm = gm.fit({fc400m,fc400std});
    
    fc_stg1 = (gm.aligned{2})
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) =fc_stg1(i,1);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = fc_stg1(i+200,1);
    end
    
    f = figure,
    BoSurfStatViewData(-heri_ct,SN,'')
    colormap(roma)
    
    bins1   = [];
    bins1   = quantileranks(fc_g1(:,1),10);
    
    for j=1:10
        for i = 1:10
            mean_stdfc_m(j,i) = mean(mean(fc400std(find(bins1==j),find(bins1==i))))
        end
    end
    
    f = figure,
    imagesc(mean_stdfc_m,[0.1 0.2])
    colormap((cbrewer('seq','Reds',99)))
    colorbar
    
    corr(gm.aligned{1}(:,1),gm.aligned{2}(:,1));
    
    f = figure,
    scatter(-gm.aligned{1}(:,1),-gm.aligned{2}(:,1), 'filled','k'),lsline
    xlim([-0.14 0.15])
  
end










