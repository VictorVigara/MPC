function evalTF(y_ss, time_constant, steps, y_discrete_steps_normalized,T_discrete_steps_init, F)
figure
title_fig = strcat("Normalised steps vs Transfer functions(y_i/",F,")" );
t = tiledlayout(2,2)
order_plots = [3; 4; 1; 2]; 

for i=1:size(y_ss,2)
    for j=1:length(steps)
        y_ss_ij = y_ss(j,i); 
        tc_ij = time_constant(j,i); 

        G_ij = tf([y_ss_ij],[tc_ij 1]); 
        [y_step, t_step] = step(G_ij, size(T_discrete_steps_init,1)); 
        y_norm = y_discrete_steps_normalized(:,i,j); 
        t_norm = T_discrete_steps_init(:,1,j)-T_discrete_steps_init(1,1,j);

        T_max = max(max(t_norm/60),max(t_step/60)); 
        Y_max = max(max(y_norm), max(y_step));
        
        % Get y values from normalised steps for same t_stamps than step
        y_n = interp1(t_norm, y_norm, t_step); 
        if isnan(y_n(end))
            y_n = y_n(1:end-1); 
            y_s = y_step(1:end-1);
        end
        % Compare outputs
        fit = goodnessOfFit(y_s, y_n, 'NRMSE'); 
        fit = (1-fit)*100;
        
        nexttile(order_plots(i))
        plot(t_step/60, y_step,'r')
        hold on
        plot(t_norm/60, y_norm,'b')
        hold on
        title_cat = strcat("Tank ", string((i)));
        title(title_cat, 'Interpreter','latex')
        xlabel("Time (min)", 'Interpreter','latex')
        ylim([0 Y_max])
        xlim([0 T_max])
        fit_text = strcat("Fit step ", string(steps(j)), ": ", string(fit),"%")
        y_med = y_step(end)/2; 
        y_size = y_step(end)/10;
        y_pos = y_med-y_size*(j-1);
        text(t_step(end)/60/2, y_pos, fit_text)
        grid on
        hold on
    end
end
leg = ["Transfer function","Normalised step"];
hl = legend(leg, 'Interpreter','latex');
hl.Layout.Tile = 'south';
title(t,title_fig, 'Interpreter','latex');
end