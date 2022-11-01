function plotData(cell_data, plot_type, title_fig)    
    % Mass plots
    n_plots = size(cell_data,1)
    leg = cell(n_plots,1);
    figure

    for i=1:n_plots
        steps = size(cell_data{i,1},3);
        for j=1:steps
            leg{i,1} = cell_data{i,3};
            X1 = []; X2 = []; X3 = []; X4 = []; T = [];
        
            if plot_type == "mass"
                max_m = max(max(cell_data{i,2}(:,:,j)/1000));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j)/1000; % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j)/1000; % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j)/1000; % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j)/1000; % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                y_label = "Mass (Kg)";
            end
            
            if plot_type == "height"
                max_m = max(max(cell_data{i,2}(:,:,j)));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j); % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j); % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j); % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j); % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                y_label = "Height (cm)";
            end
    
            if plot_type == "flow"
                max_m = max(max(cell_data{i,2}(:,:,j)));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j); % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j); % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j); % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j); % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                y_label = "q (cm3/s)";
            end
    
            % Tank 3
            subplot(2,2,1)
            plot(T, X3)
            xlabel("Time (min)")
            ylabel(y_label)
            title("Tank 3")
            ylim([0 max_m])
            grid on
            hold on
            
            % Tank 4
            subplot(2,2,2)
            plot(T, X4)
            xlabel("Time (min)")
            ylabel(y_label)
            title("Tank 4")
            ylim([0 max_m])
            grid on
            hold on
            
            % Tank 1 
            subplot(2,2,3)
            plot(T, X1)
            xlabel("Time (min)")
            ylabel(y_label)
            title("Tank 1")
            ylim([0 max_m])
            grid on
            hold on        
            
            % Tank 2
            subplot(2,2,4)
            plot(T, X2)
            xlabel("Time (min)")
            ylabel(y_label)
            title("Tank 2")
            ylim([0 max_m])
            grid on
        
            hold on
        end
    end

    % Legends
    subplot(2,2,1)
    legend(leg)
    subplot(2,2,2)
    legend(leg)
    subplot(2,2,3)
    legend(leg)
    subplot(2,2,4)
    legend(leg)

    sgtitle(title_fig)
end