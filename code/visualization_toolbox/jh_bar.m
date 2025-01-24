function [fig_bar, fig_err, fig_dot] = jh_bar(data, varargin)

% input parsing
parser = inputParser;
addOptional(parser, 'err', [], @(x) isvector(x)&isnumeric(x) )
addOptional(parser, 'sample', [], @iscell)
addParameter(parser, 'Color',[.6 .6 .6])
addParameter(parser, 'FaceAlpha', 1)

addParameter(parser, 'DrawError',true)
addParameter(parser, 'ColorError',[])
addParameter(parser, 'LineWidthError',1.2)

addParameter(parser, 'DrawMean',false)
addParameter(parser, 'DrawMeanLine',false)

addParameter(parser, 'DrawPoint',true)
addParameter(parser, 'DrawPointLine',false)
addParameter(parser, 'DotDensity',0.2)
addParameter(parser, 'PointSize',15)
addParameter(parser, 'PointAlpha',0.5)

addParameter(parser, 'DrawStats',false)
addParameter(parser, 'StatType','param')
addParameter(parser, 'Stats',[])

parse(parser, varargin{:})


dot_density = parser.Results.DotDensity;
color = parser.Results.Color;
face_alpha = parser.Results.FaceAlpha;
color_error = parser.Results.ColorError; 
line_width_error = parser.Results.LineWidthError;
is_draw_error = parser.Results.DrawError;
is_draw_mean = parser.Results.DrawMean;
is_draw_mean_line = parser.Results.DrawMeanLine;
is_draw_point = parser.Results.DrawPoint;
is_draw_point_line = parser.Results.DrawPointLine;
point_size = parser.Results.PointSize;
point_alpha = parser.Results.PointAlpha;
is_draw_stat = parser.Results.DrawStats;
stat_type = parser.Results.StatType;
stats = parser.Results.Stats;


if isempty(color_error)
    color_error = [1 1 1];
end

%%
if iscell(data)
    [avg,err] = jh_mean_err(data);
    sample = data;
elseif isvector(data)
    avg = data;
    err = parser.Results.err;
    sample = parser.Results.sample;
end
%%
hold on

% bar
fig_bar = {};
fig_err = {};
for graph_i = 1:length(avg)
    if isnan(avg(graph_i)); continue; end
    % bar
    fig = bar(graph_i, avg(graph_i));
    fig.EdgeColor = [0 0 0];
    fig.FaceAlpha = face_alpha;
    fig.LineWidth = .75; fig.EdgeColor =  [.2 .2 .2]; [.35 .35 .35]; 
    if iscell(color)
        fig.FaceColor = color{graph_i};
    else
        fig.FaceColor = color;
    end
    fig_bar{graph_i} = fig;

end

% error 
fig_err = errorbar(1:length(avg), avg, err,'marker','square');
fig_err.Color = [.25 .25 .25]; fig_err.LineStyle = 'none'; fig_err.LineWidth = line_width_error;  
fig_err.MarkerFaceColor = 'none'; fig_err.MarkerSize = 5; fig_err.MarkerEdgeColor = [.25 .25 .25]; 
if ~is_draw_error; fig_err.Color = 'none'; end
if ~is_draw_mean; fig_err.MarkerEdgeColor = 'none'; end

% sample dots
if is_draw_point
    fig_dot = {};
    x_draw_all = {};
    y_draw_all = {};
    for graph_i = 1:length(sample)
        if isempty(sample{graph_i}); continue; end

        y_draw = sort(sample{graph_i});
        y_draw(isnan(y_draw)) = [];
        y_draw = round(y_draw,5);
        epsilon = 0.000001;
        y_unique = uniquetol(y_draw, epsilon); % machine epsilon issue
        num_repeated = arrayfun(@(x) sum((x-epsilon)<y_draw & (x+epsilon)>y_draw),y_unique); % machine epsilon issue
    
        x_draw = arrayfun(@(x) linspace(-fig_bar{graph_i}.BarWidth/2*dot_density, fig_bar{graph_i}.BarWidth/2*dot_density, x - mod(x,2)),num_repeated,'UniformOutput',false);
        x_draw(mod(num_repeated,2)==1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2)==1),'UniformOutput',false);
        x_draw = cell2mat(reshape(x_draw,1,[])) + graph_i + .17;
    
        fig = scatter(x_draw, y_draw, point_size, 'k','filled');
        fig.MarkerFaceAlpha = point_alpha;
        fig_dot{graph_i} = fig;

        x_draw_all{graph_i} = x_draw;
        y_draw_all{graph_i} = y_draw;
    end

    if is_draw_point_line
        fig_line = plot(cell2mat(x_draw_all'),cell2mat(y_draw_all'),'color',[.25 .25 .25 point_alpha*0.5], 'linewidth',1);   
    end
end

% stats
if is_draw_stat
    range = ylim;
    range = diff(range);

    for graph_i = 1:length(sample)
        temp = sample{graph_i};
        if ~isempty(stats)
            p = stats(graph_i);
        elseif strcmp(stat_type, 'param')
            [~,p] = ttest(temp);
        elseif strcmp(stat_type, 'nonparam')
            p = signrank(temp);
        end

        x = graph_i;
        
        if avg(graph_i) > 0
            y = avg(graph_i) + 1.05 * err(graph_i) + .05 * range;
        else
            y = avg(graph_i) - 1.05 * err(graph_i) - .05 * range;
        end
        mark = [];
        if p < .1; mark = '+'; end
        if p < .05; mark = '*'; end
        if p < .01; mark = '**'; end
        if p < .001; mark = '***'; end
        if p < .1
            text(x, y, mark, 'FontSize',15, 'HorizontalAlignment','center','VerticalAlignment','middle')
        end

    end
end

