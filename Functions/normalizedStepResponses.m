function [normSteps_F1, normSteps_F2] = normalizedStepResponses(y,y_ss,u,u_ss,steps)
    for i=1:length(y_ss)
        for j=1:length(u_ss)
            for k=1:length(steps)
                normSteps(:,i,j,k) = (y(:,i,k) - y_ss(i))/(u(j)*steps(k) - u_ss(j));
            end
        end
    end
    normSteps_F1 = squeeze(normSteps(:,:,1,:)); 
    normSteps_F2 = squeeze(normSteps(:,:,2,:));
end