 

% compute the correlation between microstructure profile covariance and
% functional connectivity

% make the mean maps of both metrics
fc400m = squeeze(mean(fc400z(keep,:,:),1));
MPCm = squeeze(mean(MPC(keep,:,:),1));

% correlate both metrics
for i = 1:400
    [r p] = corr(fc400m(i,:)',MPCm(i,:)')
    fc_mpc(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = fc_mpc(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = fc_mpc(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))


% similar approach for the heritable components
for i = 1:400
    [r p] = corr(fc400m(i,:)',hr_fc(i,:)')
    fc_fch(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = fc_fch(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = fc_fch(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))

     
% similar approach for the heritable components
for i = 1:400
    [r p] = corr(MPCm(i,:)',hr_mpc(i,:)')
    mpc_mpch(i) = r;
end

fsa5_map = zeros(1,20484);
for i = 1:200
    fsa5_map(:,find(parcels400==i+1)) = mpc_mpch(i);
end
for i = 1:200
    fsa5_map(:,find(parcels400==i+1001)) = mpc_mpch(i+200);
end

f = figure,
BoSurfStatViewData(fsa5_map,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))

