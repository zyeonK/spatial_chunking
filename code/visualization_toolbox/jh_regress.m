function [r,p,fig_scatter] = jh_regress(xData, yData, varargin)

% argument parsing
switch nargin
    case 1
    case 2
        draw_param = 'on';
    otherwise
        draw_param = varargin{1};
        varargin(1) = [];
end

% remove nan
xData = reshape(xData,[],1);
yData = reshape(yData,[],1);
ind_remove = isnan(xData) | isnan(yData);
xData(ind_remove) = [];
yData(ind_remove) = [];

% precision issue
xData = round(xData,10);
yData = round(yData,10);

% caculate correlation
[r,p] = corr(reshape(xData,[],1),reshape(yData,[],1), varargin{:});

% plot
if strcmp(draw_param, 'dot')
    fig_scatter = scatter(xData,yData,'filled');
    fig_scatter.MarkerFaceColor = [0.5 0.5 0.5];
    fig_scatter.MarkerFaceAlpha = .6;

elseif strcmp(draw_param, 'line')
    fig_line = lsline;     
    fig_line.LineWidth = 2;

elseif strcmp(draw_param, 'ci')

    mdl = fitlm(xData, yData);
    x = linspace(min(xData), max(xData), 1000);
    [y,ci] = predict(mdl,x');

    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    line(x, y, 'color', [.75 .75 .75], 'linewidth',2);

elseif strcmp(draw_param, 'on')
    fig_scatter = scatter(xData,yData,'filled');
    fig_scatter.MarkerFaceColor = [0.5 0.5 0.5];
    fig_scatter.MarkerFaceAlpha = .6;

    mdl = fitlm(xData, yData);
    x = linspace(min(xData), max(xData), 1000);
    [y,ci] = predict(mdl,x');

    x_lim = xlim; y_lim = ylim;

    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 

    if p<0.05
        line(x, y, 'color', 'r', 'linewidth',2);
    else
        line(x, y, 'color', [.5 .5 .5], 'linewidth',2);
    end


    xlim(x_lim); ylim(y_lim);

end

if ~strcmp(draw_param, 'off')
    set(gca,'LineWidth',1.5);
    set(gca,'FontName','Helvetica','FontSize',13, 'FontWeight','bold')
    box off
end
