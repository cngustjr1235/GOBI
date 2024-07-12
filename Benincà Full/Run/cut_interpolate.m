function [y_fit, t_fit, num_data, residual] = cut_interpolate(y, window_size, overlapping_ratio, method)
%% cut

window_move = ceil(window_size * (1-overlapping_ratio));

y_cut = {};
start = 1;

length_timeseries = length(y(:,1));
num_component = length(y(1,:));

while(1)
    if start + window_size > length_timeseries
        start = length_timeseries - window_size + 1;
    end
    y_tmp = y(start:start+window_size-1, :);
    for component = 1:num_component
        y_tmp(y_tmp(:, component) < 0, component) = 0;
        % normalize data
        y_tmp(:, component) = y_tmp(:, component) - min(y_tmp(:, component));
        if max(y_tmp(:, component)) ~= 0
            y_tmp(:, component) = y_tmp(:, component) / max(y_tmp(:, component));
        end
    end
    y_cut{end+1} = y_tmp;
    if start+window_size > length_timeseries
        break;
    end
    start = start + window_move;
end

num_data = length(y_cut);

%% interpolate data (TO BE)

y_fit = y_cut;
t_fit = linspace(0, 1, window_size).';

%% calculate residual
residual_total = zeros(num_data, num_component);
for data = 1:num_data
    y_cut_tmp = y_cut{data};
    y_fit_tmp = y_fit{data};
    for component = 1:num_component
        residual_tmp = mean((y_fit_tmp(:, component)-y_cut_tmp(:, component)).^2);
        residual_total(data, component) = residual_tmp;
    end
end
residual = mean(mean(residual_total));

end

