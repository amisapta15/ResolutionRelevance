
%function subplots_main( RESULTS_H_ALL, RESULTS_H_CDM_ALL, RESULTS_H_CDM_ALL_UP, RESULTS_H_CDM_ALL_DOWN, RESULTS_H_ind_ALL, RESULTS_H_pair_ALL, sample_nums_ca, sample_nums_sub, folderPath, save_result_name, NAMES)
%[ca_ind_data_S, sub_ind_data_S, ca_pair_data_S, sub_pair_data_S, ca_ind_data_S_CDM, sub_ind_data_S_CDM, ca_pair_data_S_CDM, sub_pair_data_S_CDM] = 


sample_nums_ca = [4, 6, 8, 10, 12, 14, 16]; %
sample_nums_sub = [4, 6, 8, 10, 12, 14, 16, 19]; %4, 6, 8, 10, 12, 14, 16, 19


COLORS = [ 1, 0, 0; 0, 0, 1];
REGIONS = ['CA1' , 'SUB'];

  
count_CA = sum(contains(NAMES, 'CA'));
count_SUB = sum(contains(NAMES, 'SUB'));

color_CA = [linspace(0, 1, count_CA)', zeros(count_CA, 1), zeros(count_CA, 1)];
color_sub = [zeros(count_SUB, 1), zeros(count_SUB, 1), linspace(0, 1, count_SUB)'];



figure('Position', [100, 100, 800, 600]);

subplot(3, 1, 1);
ca=0;
sub=0;
for i = 1:numel(NAMES)
            
    
    if contains(NAMES(i), 'CA')
        ca = ca+1;
        color = color_CA(ca,:);
        
    end
    
    if contains(NAMES(i), 'SUB')
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end
      
    PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, RESULTS_H_ind_ALL{i}, 'UniformOutput', false)),  '-o', 'LineWidth',0.5, 'MarkerSize', 8, 'Color', color, 'DisplayName', [displayName,', ','H_{ind}']); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
   
    
    
    PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false)),  '-^', 'LineWidth', 0.5, 'MarkerSize', 8, 'Color', color, 'DisplayName', [displayName,', ','H','-','CDM']); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;


    PLOT = errorbar(sample_nums, cell2mat(cellfun(@mean, RESULTS_H_ind_ALL{i}, 'UniformOutput', false)), cell2mat(cellfun(@std, RESULTS_H_ind_ALL{i}, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');

    
    mean_UP = cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    grid on;


    
    
end
xlabel('Population Size');
%ylabel('$H_{ind} - H_{data}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside', 'Orientation', 'vertical');

%%%%%%%%%%%%%%%%%%%%%%%%%



subplot(3, 1, 2);
ca=0;
sub=0;
for i = 1:numel(NAMES)

    
    
    
    %Sind_min_S = cellfun(@(a, b) a - b, RESULTS_H_ind_ALL{i}, RESULTS_H_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM = cellfun(@minus, RESULTS_H_ind_ALL{i}, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM_UP = cellfun(@minus, RESULTS_H_ind_ALL{i}, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false);
    Sind_min_S_CDM_DOWN = cellfun(@minus, RESULTS_H_ind_ALL{i},  RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false);

    
    %Sind_min_S_CDM_UP = cellfun(@mean, Sind_min_S_CDM_UP{i});
    %Sind_min_S_CDM_DOWN = cellfun(@mean, Sind_min_S_CDM_DOWN{i});

    if contains(NAMES(i), 'CA')
        ca = ca+1;
        color = color_CA(ca,:);
        
    end
    
    if contains(NAMES(i), 'SUB')
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end
        

    %PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, Sind_min_S, 'UniformOutput', false)),  '-o', 'LineWidth', 1, 'MarkerSize', 8, 'MarkerEdgeColor', 'auto', 'MarkerFaceColor', 'auto', 'Color', color, 'DisplayName', displayName); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    

    
    PLOT = plot( sample_nums,cell2mat(cellfun(@mean, Sind_min_S_CDM, 'UniformOutput', false)), 'LineStyle',  '-', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, 'Color', color, 'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;



    %errorbar(sample_nums, cell2mat(cellfun(@mean, Sind_min_S, 'UniformOutput', false)), cell2mat(cellfun(@std, Sind_min_S, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
    grid on;
    hold on;
    %errorbar(sample_nums, cell2mat(cellfun(@mean, Sind_min_S_CDM, 'UniformOutput', false)), cell2mat(cellfun(@std, Sind_min_S_CDM, 'UniformOutput', false)), 'o-', 'Color', color, 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
    grid on;
    hold on;

    
    mean_UP = cell2mat(cellfun(@mean, Sind_min_S_CDM_UP, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, Sind_min_S_CDM_DOWN, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;

    
    
end
xlabel('Population Size');
ylabel('$H_{ind} - H_{data}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside');




subplot(3, 1, 3);
ca=0;
sub=0;
for i = 1:numel(NAMES)


    %norm_Sind_min_S = cellfun(@(a, b) (a - b)/b, RESULTS_H_ind_ALL{i}, RESULTS_H_ALL{i}, 'UniformOutput', false);
    
    norm_Sind_min_S_CDM = cellfun(@(a, b) (a - b)./b, RESULTS_H_ind_ALL{i}, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM_UP = cellfun(@(a, b, c) (a-c)./b, RESULTS_H_ind_ALL{i}, RESULTS_H_CDM_ALL{i}, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false);
    Sind_min_S_CDM_DOWN = cellfun(@(a, b,c) (a-c)./b, RESULTS_H_ind_ALL{i}, RESULTS_H_CDM_ALL{i}, RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false);




    if contains(NAMES(i), 'CA')
        %ca_ind_data_S = norm_Sind_min_S;
        ca_ind_data_S_CDM = norm_Sind_min_S_CDM;
        ca = ca+1;
        color = color_CA(ca,:);
    end
    
    if contains(NAMES(i), 'SUB')
        %sub_ind_data_S = norm_Sind_min_S;
        sub_ind_data_S_CDM = norm_Sind_min_S_CDM;
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end
        
    
    %PLOT = plot( sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S, 'UniformOutput', false)),  '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', color, 'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    
    PLOT = plot( sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S_CDM, 'UniformOutput', false)),  '-o', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color, 'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;

    %errorbar(sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S, 'UniformOutput', false)), cell2mat(cellfun(@std, norm_Sind_min_S, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    %errorbar(sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S_CDM, 'UniformOutput', false)), cell2mat(cellfun(@std, norm_Sind_min_S_CDM, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    
    grid on;
    hold on;

    mean_UP = cell2mat(cellfun(@mean, Sind_min_S_CDM_UP, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, Sind_min_S_CDM_DOWN, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;


    
    
end
xlabel('Population Size');
ylabel('$\frac{H_{ind} - H_{data}}{H_{data}}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside');


spacing = 0.03;
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
set(subplot(3, 1, 2), 'Position', get(subplot(3, 1, 2), 'Position') - [0, 0, 0, spacing]);
set(subplot(3, 1, 3), 'Position', get(subplot(3, 1, 3), 'Position') - [0, spacing, 0, 0]);




hold off;
saveas(gcf, [folderPath, save_result_name, '_ind'], 'jpg');
saveas(gcf, [folderPath , save_result_name, '_ind'], 'fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i = 1:numel(NAMES)
%    RESULTS_H_pair_ALL{i} = cellfun(@(x) cell2mat(x), RESULTS_H_pair_ALL{i}, 'UniformOutput', false)
%end


figure('Position', [100, 100, 800, 600]);
subplot(3, 1, 1);
ca=0;
sub=0;
for i = 1:numel(NAMES)

    
    if contains(NAMES(i), 'CA')
        ca = ca+1;
        color = color_CA(ca,:);
        
    end
    
    if contains(NAMES(i), 'SUB')
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end

    PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, RESULTS_H_pair_ALL{i}, 'UniformOutput', false)),  '-o', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color, 'DisplayName', [displayName,', ','H_{pair}']); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    
    
    %PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, RESULTS_H_ALL{i}, 'UniformOutput', false)),  '-s', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color, 'DisplayName', [displayName,', ','H']); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    
    
    PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false)),  '-^', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color, 'DisplayName', [displayName,', ','H','-','CDM']); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;

  
    %errorbar(sample_nums, cell2mat(cellfun(@mean, RESULTS_H_pair_ALL{i}, 'UniformOutput', false)), cell2mat(cellfun(@std, RESULTS_H_pair_ALL{i}, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    %errorbar(sample_nums, cell2mat(cellfun(@mean, RESULTS_H_ALL{i}, 'UniformOutput', false)), cell2mat(cellfun(@std, RESULTS_H_ALL{i}, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    %errorbar(sample_nums, cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false)), cell2mat(cellfun(@std, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');

    grid on;
    hold on;

    mean_UP = cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none',  'DisplayName', displayName);
    hold on;


    
    
end
xlabel('Population Size');
%ylabel('$H_{pair} - H_{data}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside', 'Orientation', 'vertical');

%%%%%%%%%%%%%%%%%%%%%%%%%



subplot(3, 1, 2);
ca=0;
sub=0;
for i = 1:numel(NAMES)

    
    
    
    %Sind_min_S = cellfun(@(a, b) a - b, RESULTS_H_pair_ALL{i}, RESULTS_H_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM = cellfun(@minus, RESULTS_H_pair_ALL{i}, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM_UP = cellfun(@minus, RESULTS_H_pair_ALL{i}, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false);
    Sind_min_S_CDM_DOWN = cellfun(@minus, RESULTS_H_pair_ALL{i},  RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false);

    
    %Sind_min_S_CDM_UP = cellfun(@mean, Sind_min_S_CDM_UP{i});
    %Sind_min_S_CDM_DOWN = cellfun(@mean, Sind_min_S_CDM_DOWN{i});

    if contains(NAMES(i), 'CA')
        ca = ca+1;
        color = color_CA(ca,:);
        
    end
    
    if contains(NAMES(i), 'SUB')
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end
        

    %PLOT = plot( sample_nums,  cell2mat(cellfun(@mean, Sind_min_S, 'UniformOutput', false)),  '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', color, 'DisplayName', displayName); 
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    

    
    PLOT = plot( sample_nums,cell2mat(cellfun(@mean, Sind_min_S_CDM, 'UniformOutput', false)),  '-o', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color,  'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    

    %errorbar(sample_nums, cell2mat(cellfun(@mean, Sind_min_S, 'UniformOutput', false)), cell2mat(cellfun(@std, Sind_min_S, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    %errorbar(sample_nums, cell2mat(cellfun(@mean, Sind_min_S_CDM, 'UniformOutput', false)), cell2mat(cellfun(@std, Sind_min_S_CDM, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    grid on;
    hold on;


    mean_UP = cell2mat(cellfun(@mean, Sind_min_S_CDM_UP, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, Sind_min_S_CDM_DOWN, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;


    
    
end
xlabel('Population Size');
ylabel('$H_{pair} - H_{data}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside');




subplot(3, 1, 3);
ca=0;
sub=0;
for i = 1:numel(NAMES)


    %norm_Sind_min_S = cellfun(@(a, b) (a - b)/b, RESULTS_H_pair_ALL{i}, RESULTS_H_ALL{i}, 'UniformOutput', false);
    norm_Sind_min_S_CDM = cellfun(@(a, b) (a - b)./b, RESULTS_H_pair_ALL{i}, RESULTS_H_CDM_ALL{i}, 'UniformOutput', false);
    Sind_min_S_CDM_UP = cellfun(@(a, b, c) (a-c)./b, RESULTS_H_pair_ALL{i}, RESULTS_H_CDM_ALL{i}, RESULTS_H_CDM_ALL_UP{i}, 'UniformOutput', false);
    Sind_min_S_CDM_DOWN = cellfun(@(a, b,c) (a-c)./b, RESULTS_H_pair_ALL{i}, RESULTS_H_CDM_ALL{i}, RESULTS_H_CDM_ALL_DOWN{i}, 'UniformOutput', false);




    if contains(NAMES(i), 'CA')
        %ca_pair_data_S = norm_Sind_min_S;
        ca_pair_data_S_CDM = norm_Sind_min_S_CDM;
        ca = ca+1;
        color = color_CA(ca,:);
    end
    
    if contains(NAMES(i), 'SUB')
        %sub_pair_data_S = norm_Sind_min_S;
        sub_pair_data_S_CDM = norm_Sind_min_S_CDM;
        sub = sub+1;
        color = color_sub(sub,:);
    end
    displayNameParts = strsplit(char(NAMES(i)), '_');
    displayName = strjoin(displayNameParts, ',');
    displayName = strrep(displayName, '.mat', '');
    displayName = strrep(displayName, '_', '\_');
    
    if contains(NAMES(i), 'CA')
        sample_nums = sample_nums_ca;
    end
    if contains(NAMES(i), 'SUB')
        sample_nums = sample_nums_sub;
    end
        
    
    %PLOT = plot( sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S, 'UniformOutput', false)),  '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', color, 'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;
    
    PLOT = plot( sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S_CDM, 'UniformOutput', false)),  '-o', 'LineWidth', 1, 'MarkerSize', 8, 'Color', color,  'DisplayName', displayName);
    set(gca, 'XTick', 1:max(sample_nums)); % Set the x-axis ticks to integers
    grid on;
    hold on;

    %errorbar(sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S, 'UniformOutput', false)), cell2mat(cellfun(@std, norm_Sind_min_S, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    %errorbar(sample_nums, cell2mat(cellfun(@mean, norm_Sind_min_S_CDM, 'UniformOutput', false)), cell2mat(cellfun(@std, norm_Sind_min_S_CDM, 'UniformOutput', false)),  'o-', 'Color', color, 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
    grid on;
    hold on;


    mean_UP = cell2mat(cellfun(@mean, Sind_min_S_CDM_UP, 'UniformOutput', false));
    mean_DOWN = cell2mat(cellfun(@mean, Sind_min_S_CDM_DOWN, 'UniformOutput', false));
    fill([sample_nums, fliplr(sample_nums)], [mean_DOWN,  fliplr(mean_UP)], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;


    
    
end
xlabel('Population Size');
ylabel('$\frac{H_{pair} - H_{data}}{H_{data}}$', 'Interpreter', 'latex');
legend('Location', 'northeastoutside');


spacing = 0.03;
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
set(subplot(3, 1, 2), 'Position', get(subplot(3, 1, 2), 'Position') - [0, 0, 0, spacing]);
set(subplot(3, 1, 3), 'Position', get(subplot(3, 1, 3), 'Position') - [0, spacing, 0, 0]);




hold off;
saveas(gcf, [folderPath , save_result_name, '_pair'], 'jpg');
saveas(gcf, [folderPath , save_result_name, '_pair'], 'fig');



%end




