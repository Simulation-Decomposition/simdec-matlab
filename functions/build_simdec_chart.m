
function b = build_simdec_chart(output,scenario,scenario_legend,main_colors,axistitle,var_names_dec)
% significance calculates how much variability of the output is explained by inputs 
% 
% Uses function hex2rgb.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% REQUIRED INPUTS:
%
%   output       - an array of Y values
%   scenario     - an array of associated scenario indices
%               
% OPTIONAL INPUTS (specify empty matrix [] if not using):
%
%   main_colors      - hex values of the main colors, should be as many as the
%                      states for the first in decomposition (or the most
%                      influential input variable). For example,
%                      main_colors = {'#00B0F0'; '#1FDF4D'; '#1FDF4D'}
%   scenario_legend  - a scenario table that shows which states of which
%                      variables compose different sceanrios
%   axistitle        - a title for the x axis.
%   var_name_dec     - cell array with names for variables for decomposition
%                      (scenario legend)
%
% OUTPUT
%   the SimDec graph and the scenario table.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova & ChatGPT, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



% separating result by scenario
result_dec=cell(1,max(scenario));

for i=1:max(scenario)
    s=size(output(scenario==i),1);
    result_dec(i)=mat2cell(output(scenario==i),s,1);
end

% defining edges of bins
l=min(output);
h=max(output);
bins=100;
edges=[l:((h-l)/bins):h];

% define frequency of each scenario NPV for each bin
f=zeros(max(scenario),100);
for i=1:max(scenario)
    f(i,:)=histcounts(cell2mat(result_dec(i)),edges);
end

% building stacked histogram
figure
b=bar(1:bins,permute(f,[2 1]),1,'stack');

% axis titles
xlabel(axistitle);
ylabel('Probability');

% colours
    
    % main colors 
    if isempty(main_colors)
        main_colors = {'#00B0F0'; '#E7D819'; '#1FDF4D'};
    end 

    N_main_colors = max(scenario_legend(:,2));
    N_shades = scenario_legend(end,1) / N_main_colors;
    color = NaN(scenario_legend(end,1),3);

        % shades
        sc = 1;
        if N_shades >= 3
            range = 0.8; 
        else
            range = 0.5;
        end
        step = range * 2 / (N_shades - 1);
        beta = [-range : step : range];
        
        for m = 1 : N_main_colors
            current_main_color = hex2rgb(main_colors(m));
            
            for sh = 1 : N_shades
                color(sc,:) = brighten (current_main_color, beta(sh));
                sc = sc+1;
            end
        end

        %filling

        for i=1:max(scenario)
            b(i).FaceColor = color(i,:);
        end
        
        
        for i=1:max(scenario)
            b(i).EdgeColor = [0.5 0.5 0.5];
        end
        hold on

% ticks
% real axis ticks
figure
h=histogram(output,bins);
xticks=get(gca,'xtick');
% to fix matlab R2022b update issue, when the last X tick is not
%   necesseraly at the very end
if max(h.BinEdges) > max(xticks)
    xticks = [xticks, xticks(end)+xticks(2)]; 
end
yticks=get(gca,'ytick');
close
% X axis

% getting corresponding min and max on the artificial axis
min_art=(xticks(1)-min(output))*bins/(max(output)-min(output));
max_art=(xticks(end)-max(output))*bins/(max(output)-min(output))+bins;
distance_art=max_art-min_art;

% corresponding ticks on the artificial axis
xticks_art=(xticks-xticks(1)).*distance_art./(xticks(end)-xticks(1))+min_art;

% Y axis
total=size(output,1);
yticks_art=yticks/total;
a=[cellstr(num2str(yticks_art'*100))]; % converting values into percentage
pct = char(ones(size(a,1),1)*'%'); % creating vector of % signs
new_yticks = [char(a),pct]; % add the '%' signs after the percentage values

% getting new labels to the graph
set(gca,'XLim',[min(xticks_art), max(xticks_art)],'XTick',xticks_art,'XTickLabel',xticks,'YLim',[0 yticks(end)],'YTickLabel',new_yticks);
hold off


%% showing the legend

    % rename states
    
var1 = {'Color'};
var2 = var_names_dec;
var3 = {'Ymin','Ymean','Ymax', 'Probability'};
var = [var1  var2  var3];


t = array2table(scenario_legend,'VariableNames',var);
fig = uifigure;
uit = uitable(fig,'Data', t);

    for i = 1 : scenario_legend(end,1)
        s = uistyle('BackgroundColor', color(i,:));
        addStyle(uit,s,'cell',[i,1]);
    
    end


end
