function winsorized_data = winsorize_by_sd(data,n_sd)

if nargin < 2
    n_sd = 3;  % Default is to winsorize data that exceeds 3 standart deviations
end

winsorized_data = data;
sd = std(data,'omitnan');
m = mean(data,'omitnan');
winsorized_data(winsorized_data < (m - (n_sd*sd)) ) = m - (n_sd*sd);
winsorized_data(winsorized_data > (m + (n_sd*sd)) ) = m + (n_sd*sd);

end
