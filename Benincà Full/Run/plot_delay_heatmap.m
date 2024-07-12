function [] = plot_delay_heatmap(score)
    
    c_nan = [0.9 0.9 0.9];
    font_s = 14;
    tmp = 1 - linspace(0, 1, 501)';

    % positive only
    cmap_score = [[tmp], [tmp], [ones(501,1)]];
    score_range = [0,1];

    % negative only
    % cmap_score = [[ones(501,1)],[1-tmp],[1-tmp]];    
    % score_range = [-1,0];

    % positive and negative
    % cmap_score = [[ones(500,1);tmp],[1-tmp(1:end-1);tmp],[1-tmp; ones(500,1)]];
    % score_range = [-1,1];

    
    x_var = 0:10;
    y_var = 0:10;

    x_label = "Delay of Rotifers(day)";    
    y_label = "Delay of Protozoa(day)";

    
    h = heatmap(score);
    
    h.FontName = 'Arial';
    h.Colormap = cmap_score;
    h.ColorLimits = score_range;
    h.FontSize = font_s;
    h.MissingDataColor = c_nan;
    h.XData = x_var;
    h.YData = y_var;
    h.XLabel = x_label;
    h.YLabel = y_label;
    s = struct(h);
    s.Axes.XAxisLocation = "top";
    cbh = s.Colorbar;
    set(cbh, 'YTick', score_range)
end