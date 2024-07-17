function createStackedHistogram(data, xlim_values, ylim_values, x_axis_title, scenarios, colors, n_bars, rotate90, hide_axes)
    % createStackedHistogram - Function to create a stacked histogram.
    %
    % Inputs:
    %   data          - Data for the histogram.
    %   xlim_values   - Limits for the x-axis [min, max]. If empty, automatic scaling is used.
    %   ylim_values   - Limits for the y-axis [mi, max] (as probabilities in %). If empty, automatic scaling is used.
    %   x_axis_title  - Title for the x-axis.
    %   scenarios     - Scenario indices corresponding to the data.
    %   colors        - Colors for each scenario.
    %   n_bars        - Number of bins for the histogram.
    %   rotate90      - Boolean indicating whether to rotate the plot 90 degrees.
    %   hide_axes     - Boolean indicating whether to hide both axes.

    % Number of scenarios
    uniqueScenarios = unique(scenarios);
    numScenarios = length(uniqueScenarios);

    % Define histogram bins
    if isempty(xlim_values)
        figure
        s = scatter(data, data);
        xlim_values = get(gca, 'XLim');
        close(gcf)
    end

    edges = linspace(xlim_values(1), xlim_values(2), n_bars + 1);

    % Initialize matrix to store counts
    counts = zeros(length(edges)-1, numScenarios);

    % Calculate histogram counts for each scenario
    for i = 1:numScenarios
        idx = scenarios == uniqueScenarios(i);
        counts(:, i) = histcounts(data(idx), edges);
    end

    % Convert counts to probabilities (percentage)
    probabilities = (counts / length(data)) * 100;

    % Create the stacked histogram
    if rotate90
        % For rotated histogram
        h = barh(edges(1:end-1) + diff(edges)/2, probabilities, 'stacked', 'BarWidth', 1);
        ylabel(x_axis_title);
        xlabel('Probability');
        xtickformat('percentage');
    else
        % For normal histogram
        h = bar(edges(1:end-1) + diff(edges)/2, probabilities, 'stacked', 'BarWidth', 1);
        xlabel(x_axis_title);
        ylabel('Probability');
        ytickformat('percentage');
    end

    % Apply colors to bars
    for i = 1:numScenarios
        h(i).FaceColor = colors(i, :);
    end

    % Set x-axis and y-axis limits
    if rotate90
        ylim(xlim_values);
        if ~isempty(ylim_values)
            xlim(ylim_values);
        end        
    else
        xlim(xlim_values);
        if ~isempty(ylim_values)
            ylim(ylim_values);
        end        
    end

    % Hide axes if specified
    if hide_axes
        if rotate90
            ylabel('');
            xlabel('');
            set(gca, 'XColor', 'w', 'xtick', [], 'ytick', []);
        else
            xlabel('');
            ylabel('');
            set(gca, 'YColor', 'w', 'xtick', [], 'ytick', []);
        end
    end

end
