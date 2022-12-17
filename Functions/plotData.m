function plotData(cell_data, plot_type, title_fig)    
    % Mass plots
    n_plots = size(cell_data,1);
    leg = cell(n_plots,1);
    figure
    t = tiledlayout(2,2);

    for i=1:n_plots
        leg{i,1} = cell_data{i,3};
        X1 = []; X2 = []; X3 = []; X4 = []; T = [];
    
        if plot_type == "mass"
            max_m = max(max(cell_data{i,2}(:,:,1)/1000));  % Max mass to ylim
            X1 = cell_data{i,2}(:,1,1)/1000; % Mass tank 3 in Kg
            X2 = cell_data{i,2}(:,2,1)/1000; % Mass tank 3 in Kg
            X3 = cell_data{i,2}(:,3,1)/1000; % Mass tank 3 in Kg
            X4 = cell_data{i,2}(:,4,1)/1000; % Mass tank 3 in Kg
            T = cell_data{i,1}(:,:,1)/60;        % Time in min
            y_label = "Mass (Kg)";
        end
        
        if plot_type == "height"
            max_m = max(max(cell_data{i,2}(:,:)));  % Max mass to ylim
            X1 = cell_data{i,2}(:,1,1); % Mass tank 3 in Kg
            X2 = cell_data{i,2}(:,2,1); % Mass tank 3 in Kg
            X3 = cell_data{i,2}(:,3,1); % Mass tank 3 in Kg
            X4 = cell_data{i,2}(:,4,1); % Mass tank 3 in Kg
            T = cell_data{i,1}(:,:,1)/60;        % Time in min
            y_label = "Height (cm)";
        end

        if plot_type == "flow"
            max_m = max(max(cell_data{i,2}(:,:)));  % Max mass to ylim
            X1 = cell_data{i,2}(:,1,1); % Mass tank 3 in Kg
            X2 = cell_data{i,2}(:,2,1); % Mass tank 3 in Kg
            X3 = cell_data{i,2}(:,3,1); % Mass tank 3 in Kg
            X4 = cell_data{i,2}(:,4,1); % Mass tank 3 in Kg
            T = cell_data{i,1}(:,:,1)/60;        % Time in min
            y_label = "q (cm3/s)";
        end

        % Tank 3
        nexttile(1)
        plot(T, X3)
        xlabel("Time (min)", 'Interpreter','latex')
        ylabel(y_label, 'Interpreter','latex')
        title("Tank 3", 'Interpreter','latex')
        ylim([0 max_m])
        grid on
        hold on
        
        % Tank 4
        nexttile(2)
        plot(T, X4)
        xlabel("Time (min)", 'Interpreter','latex')
        ylabel(y_label, 'Interpreter','latex')
        title("Tank 4", 'Interpreter','latex')
        ylim([0 max_m])
        grid on
        hold on
        
        % Tank 1 
        nexttile(3)
        plot(T, X1)
        xlabel("Time (min)", 'Interpreter','latex')
        ylabel(y_label, 'Interpreter','latex')
        title("Tank 1", 'Interpreter','latex')
        ylim([0 max_m])
        grid on
        hold on        
        
        % Tank 2
        nexttile(4)
        plot(T, X2)
        xlabel("Time (min)", 'Interpreter','latex')
        ylabel(y_label, 'Interpreter','latex')
        title("Tank 2", 'Interpreter','latex')
        ylim([0 max_m])
        grid on
    
        hold on
        
    end

    % Legends
    hl = legend(leg, 'Interpreter','latex');
    hl.Layout.Tile = 'south';
    title(t,title_fig, 'Interpreter','latex');
end