function p_perm_fdr = get_perm_fdr_p(original_t, permutation_t_values,n_parcels)
p_prop_pos = nan(n_parcels,1);p_prop_neg = nan(n_parcels,1);
for i = 1:n_parcels
    a = original_t(i); %true t-values
    b = permutation_t_values(i,:); %t-values from 10000 (or n_parcels) permutations
    p_prop_pos(i) = sum(b<a) / size(permutation_t_values(i,:),2);
    p_prop_neg(i) = sum(b>a) / size(permutation_t_values(i,:),2);
end
p_perm_fdr = fdr_bh(p_prop_pos,0.025) + fdr_bh(p_prop_neg,0.025); %two tailed testing, FDR level 0.025 on each side
