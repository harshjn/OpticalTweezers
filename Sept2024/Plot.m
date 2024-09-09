T=300         %T #Kelvin is 25degreeCelsius
r=10e-6        %Radius of the outer circle
periodR = 2*pi*r;
k_b=1.38e-23;
kBT = k_b*T;
%%
Alpha(:,1)=10*ones(20,1);
Alpha(:,2)=100*ones(20,1);
Alpha(:,3)=1000*ones(20,1);
F0Mat = Alpha;

F0Mat(:,1) = linspace(1e-15,10e-15,20);
F0Mat(:,2) = linspace(1e-15,10e-15,20)*10;
F0Mat(:,3) = linspace(1e-15,10e-15,20)*100;

jMat = F0Mat;

%% After inputting j Values from the three output files.
LambdaMat = F0Mat*periodR/kBT;

for k = 1:3
    figure
    plot( LambdaMat(:,k)./Alpha(:,k), jMat(:,k)./Alpha(:,k))

end
%%

% Assuming LambdaMat, jMat, and Alpha are already calculated and available

% Create figure with high resolution
figure('Position', [100, 100, 800, 600]);

% Color scheme
colors = {'#1f77b4', '#ff7f0e', '#2ca02c'}; % Blue, Orange, Green
labels = {'\alpha=10', '\alpha=100', '\alpha=1000'};

% Plot data
for k = 1:2
    plot(LambdaMat(:,k)./Alpha(:,k), jMat(:,k)./Alpha(:,k),'o-', ...
        'Color', colors{k}, 'LineWidth', 8);
    hold on;
end
for k = 3
    plot(LambdaMat(:,k)./Alpha(:,k), jMat(:,k)./Alpha(:,k),'o-', ...
        'Color', colors{k}, 'LineWidth', 4);
    hold on;
end

% Customize the plot
%set(gca, 'YScale', 'log');
xlabel('\lambda / \alpha', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('j / \alpha', 'FontSize', 34, 'Interpreter', 'tex');
title('Custom Plot', 'FontSize', 36);

% Improve tick labels
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.TickLength = [0.02 0.02];
ax.FontSize = 12;

% Add grid
grid on;
ax.GridAlpha = 0.2;
ax.MinorGridAlpha = 0.1;

% Add legend
legend(labels, 'FontSize', 12, 'Location', 'best');
legend box off;

% Add text for alpha
text(0.05, 0.95, '\alpha', 'Units', 'normalized', 'FontSize', 14, ...
    'VerticalAlignment', 'top', 'Interpreter', 'tex');

% Adjust layout
set(gcf, 'Color', 'w');
box on;

% Save the figure (optional)
print('-dpng', '-r300', 'high_quality_plot.png');


