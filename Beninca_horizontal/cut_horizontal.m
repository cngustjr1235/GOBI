function [y_cut, t_cut, class_pair_list, class_count_list] = cut_horizontal(y, t, boundary)

num_component = length(boundary);
num_boundary_list = zeros(1, num_component);
for idx_comp = 1:num_component
    num_boundary_list(idx_comp) = length(boundary{idx_comp});
end
num_class_list = num_boundary_list + 1;

num_class_pair = prod(num_class_list);
class_pair_list = zeros(num_class_pair, num_component+1);

num_elem_repeat = 1;
num_total_repeat = num_class_pair;

for idx_comp = num_component:-1:1
    num_class = num_class_list(idx_comp);
    num_total_repeat = num_total_repeat / num_class;
    class_pair_list(:,idx_comp) = repmat(repelem(1:num_class, num_elem_repeat), 1, num_total_repeat).';
    num_elem_repeat = num_elem_repeat * num_class;
end
    
length_timeseries = length(t);

y_category = zeros(length_timeseries, num_component);

for idx_comp = 1:num_component
    y_comp = y(:,idx_comp);
    boundary_comp = boundary{idx_comp};
    category = 1;

    for bound = boundary_comp
        y_category(y_comp >= bound, idx_comp) = category;
        y_comp(y_comp >= bound) = NaN;
        category = category + 1;
    end
    y_category(~isnan(y_comp), idx_comp) = category;
end

y_cut = cell(1, num_class_pair);
t_cut = cell(1, num_class_pair);

[~, idx_pair_list] = ismember(y_category, class_pair_list(:,1:num_component), "rows");

for idx_pair = 1:num_class_pair
    t_cut{idx_pair} = t(idx_pair_list == idx_pair);
    y_cut{idx_pair} = y(idx_pair_list == idx_pair, :);
end

class_count_list = {};
for idx_comp = 1:num_component
    class_count = zeros(1, num_class_list(idx_comp));
    for class = 1:num_class_list(idx_comp)
        class_indices = find(class_pair_list(:,idx_comp) == class);
        class_count(class) = sum(arrayfun(@(idx) length(y_cut{idx}), class_indices));
    end
    class_count_list{idx_comp} = class_count;
end

for idx_pair = 1:num_class_pair
    class_pair_list(idx_pair, num_component+1) = length(t_cut{idx_pair});
end

end



