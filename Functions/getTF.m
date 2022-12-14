function [y_ss, time_constant] = getTF(norm_step,steps,time)
fprintf("\nTRANSFER FUNCTIONS FROM NORMALISED STEP\n")
y_ss = zeros(length(steps), size(norm_step,2)); 
time_constant = zeros(length(steps),size(norm_step,2)); 

for i=1:size(norm_step,2)
    for j=1:length(steps)
        y_step = norm_step(:,i,j);  % Get normalised step response from yi step j
        time_step = time;
        % Remove repeated values  in y_step that are not
        % valid for latter interp1 function
        if length(y_step) ~= length(unique(y_step)) 
            [y_step, idx] = unique(y_step); 
            time_step = time(idx); 
        end
        y_ss(j,i) = max(y_step);            % Get steady state gain
        y_time_constant = 0.632*y_ss(j,i);  % Get gain at 63.2% to determine time constant
        % Get time constant in seconds by interpolating previous value  
        time_constant(j,i) = interp1(y_step, time_step, y_time_constant) - time_step(1);
        fprintf("Transfer function from F1 to Y%d - Step = %f", i, steps(j))
        G = tf(y_ss(j,i), [time_constant(j,i) 1])
    end
end

end