function plotSteps(cell_data, plot_type, steps)    
    % Mass plots
    n_plots = size(cell_data,1);
    leg = cell(length(steps),1);

    for i=1:length(steps)
        leg{i} = strcat("Step ",num2str(steps(i)));
    end

    for i=1:n_plots
        figure
        t = tiledlayout(2,2)

        title_fig = strcat(cell_data{i,3}," Steps");

        for j = 1:size(cell_data{i,2},3)
            X1 = []; X2 = []; X3 = []; X4 = []; T = [];
        
            if plot_type == "mass"
                max_m = max(max(cell_data{i,2}(:,:,j)/1000));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j)/1000; % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j)/1000; % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j)/1000; % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j)/1000; % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                x_lim = max(T);
                y_label = "Mass (Kg)";
            end
            
            if plot_type == "height"
                max_m = max(max(cell_data{i,2}(:,:)));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j); % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j); % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j); % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j); % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                x_lim = max(T);
                y_label = "Height (cm)";
            end
    
            if plot_type == "flow"
                max_m = max(max(cell_data{i,2}(:,:)));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j); % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j); % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j); % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j); % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                x_lim = max(T);
                y_label = "q (cm3/s)";
            end

            if plot_type == "height"
                max_m = max(max(cell_data{i,2}(:,:)));  % Max mass to ylim
                X1 = cell_data{i,2}(:,1,j); % Mass tank 3 in Kg
                X2 = cell_data{i,2}(:,2,j); % Mass tank 3 in Kg
                X3 = cell_data{i,2}(:,3,j); % Mass tank 3 in Kg
                X4 = cell_data{i,2}(:,4,j); % Mass tank 3 in Kg
                T = cell_data{i,1}(:,:,j)/60;        % Time in min
                x_lim = max(T);
                y_label = " ";
            end
    
            % Tank 3
            nexttile(1)
            plot(T, X3)
            title("Tank 3", 'Interpreter','latex')
            xlabel("Time (min)", 'Interpreter','latex')
            ylabel(y_label, 'Interpreter','latex')
            ylim([0 max(X3)])
            xlim([T(1) x_lim])
            grid on
            hold on
            
            % Tank 4
            nexttile(2)
            plot(T, X4)
            title("Tank 4", 'Interpreter','latex')
            xlabel("Time (min)", 'Interpreter','latex')
            ylabel(y_label, 'Interpreter','latex')
            ylim([0 max(X4)])
            xlim([T(1) x_lim])
            grid on
            hold on
            
            % Tank 1 
            nexttile(3)
            plot(T, X1)
            title("Tank 1", 'Interpreter','latex')
            xlabel("Time (min)", 'Interpreter','latex')
            ylabel(y_label, 'Interpreter','latex')
            ylim([0 max(X1)])
            xlim([T(1) x_lim])
            grid on
            hold on        
            
            % Tank 2
            nexttile(4)
            plot(T, X2)
            xlabel("Time (min)", 'Interpreter','latex')
            ylabel(y_label, 'Interpreter','latex')
            title("Tank 2", 'Interpreter','latex')
            ylim([0 max(X2)])
            xlim([T(1) x_lim])
            grid on
        
            hold on
        end  
        
        hl = legend(leg, 'Interpreter','latex');
        hl.Layout.Tile = 'south';
        title(t,title_fig, 'Interpreter','latex');
    end



end