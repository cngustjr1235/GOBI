function [] = plot_heatmap(dimension, component_list, score)
    labels = ["Cyclopoids", "Calanoids", "Rotifers", "Protozoa", "Nanophytoplankton", "Picophytoplankton", "FilamentousDiatoms", "Ostracods", "Harpacticoids", "Bacteria", "Nitrogen", "Phosphorus"];
    
    c_nan = [0.9 0.9 0.9];
    font_s = 14;
    tmp = 1 - linspace(0, 1, 501)';
    cmap_score = [[tmp],[tmp],[ones(501,1)]];
    
    x_var = [];
    for i = 0:(2^dimension-1)
        x_var_tmp = [];
        for j = 1:dimension
            if rem(fix(i/(2^(dimension-j))), 2) == 0
                x_var_tmp = [x_var_tmp, "+"];
            else
                x_var_tmp = [x_var_tmp, "-"];
            end
        end
        x_var = [x_var, strjoin(x_var_tmp, " ")];
    end
    
    y_var = [];
    for i = 1:length(score)
        y_label = [];
        for j = 1:dimension+1
            idx = component_list(i, j);
            y_label = [y_label, labels(idx)];
        end
        y_var = [y_var, strjoin(y_label, ' -> ')];
    end
    
    y_label = [];
    for i = 1:dimension
        y_label = [y_label, sprintf("C%d", i)];
    end
    y_label = strjoin([y_label, "T"], ", ");
    
    h = heatmap(score);
    
    score_range = [0,1];
    h.FontName = 'Arial';
    h.Colormap = cmap_score;
    h.ColorLimits = score_range;
    h.FontSize = font_s;
    h.MissingDataColor = c_nan;
    h.XData = x_var;
    h.YData = y_var;
    h.YLabel = y_label;
    s = struct(h);
    cbh = s.Colorbar;
    set(cbh, 'YTick', score_range)
end