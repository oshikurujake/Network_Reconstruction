warning('off')
clear all
% Configeration of the resistor grid
hn = 5;
vn = hn;
R = 100;
rn = numel(ones([vn,hn-1]))+numel(ones([vn-1,hn]));
all_r = 0.5*R*rand([1,rn]);
average = mean(all_r,'all');


HR = [];
for i = 1:vn
    a = all_r(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
    HR = [HR;a];
end
VR = [];
for i = 1:(vn-1)
    b = all_r(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
    VR = [VR;b];
end

% Forward solution using kirchhoff method
[U_meas,boundary_ratio,bn,CM] = kirchhoff_grid(hn,vn,HR,VR);
% Hyperparameter
hyper = 0.03%linspace(2e-3,1e-1,50);
accu_error = [];
% Set a small pentrubation for each resistance 0.1 ohm
delta = mean(all_r,'all')/10000;
for lambda = hyper
    % Count how many interations
    count = 1;
    % Run Resistance_matrix
    initial = R*ones([1,rn]);
    HR_sim = [];
    for i = 1:vn
        a = initial(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
        HR_sim = [HR_sim;a];
    end
    VR_sim = [];
    for i = 1:(vn-1)
        b = initial(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
        VR_sim = [VR_sim;b];
    end
    %% Modify afterwards  
    Fwd_u = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
    u_diff = [];
    
    % Generate jacobian
    for i = 1:length(initial)
        new = initial;
        new(i) = initial(i) + delta;
        HR_sim = [];
        for i = 1:vn
            a = new(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
            HR_sim = [HR_sim;a];
        end
        VR_sim = [];
        for i = 1:(vn-1)
            b = new(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
            VR_sim = [VR_sim;b];
        end
        data1 = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
        % Obtain the voltage after small pentrubation delta
        u_diff = [u_diff,data1-Fwd_u];
    end
    J = u_diff/delta;

    % Obtain the resistance updates.
    del_r = round(inv(J.'*J+lambda^2*eye(rn,rn))*J.'*(U_meas-Fwd_u),8);
    
%     % L2 Norm
%     Euc_norm = norm(U_meas-Fwd_u);
%     s_norm = [Euc_norm];
%     
%     % Relative measurememnts error
%     sum_resi = sqrt(sum((Fwd_u./U_meas -1).^2,'all')/length(data1));
%     s_resi = [sum_resi];
%     
%     % Standard deviation
%     stand = std(U_meas-Fwd_u);
%     s_std = [stand];
    
    % Average/Macro error for the resistor.
    error = sqrt(sum(abs(initial./all_r-1),'all')/rn);
    ave_error = [error];
    result = initial;
    
    % Individual/Microscopic error for each element
    i_error = abs(initial./all_r-1);
    indi_error = [i_error];
    while error > 0.05
        
        % Complete the update
        initial = initial+del_r';
        result = [result;initial];
        
        % Transfer to a grid
        HR_sim = [];
        for i = 1:vn
            a = initial(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
            HR_sim = [HR_sim;a];
        end
        VR_sim = [];
        for i = 1:(vn-1)
            b = initial(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
            VR_sim = [VR_sim;b];
        end
 
        Fwd_u = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
%         Euc_norm = norm(U_meas-Fwd_u);
        
        J = [];
        u_diff = [];
        % Generate jacobian
        for i = 1:length(initial)
%             data1 = NaN;
            new = initial;
            delta_i = delta;
                
                new(i) = initial(i) + delta_i;
                HR_sim = [];
                for i = 1:vn
                    a = new(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
                    HR_sim = [HR_sim;a];
                end
                VR_sim = [];
                for i = 1:(vn-1)
                    b = new(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
                    VR_sim = [VR_sim;b];
                end
                [data1,ratio,bound_nodes,CM] = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
                while isnan(data1) ==  1
                    disp('small pentrubation cause NAN');
                    delta_i = delta_i+delta;
                    new = initial;
    %                 delta_i = delta;

                    new(i) = initial(i) + delta_i;
                    HR_sim = [];
                    for i = 1:vn
                        a = new(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
                        HR_sim = [HR_sim;a];
                    end
                    VR_sim = [];
                    for i = 1:(vn-1)
                        b = new(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
                        VR_sim = [VR_sim;b];
                    end
                    [data1,ratio,bound_nodes,CM] = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
    %                 delta_i = 2* delta_i;
                end
            % Obtain the voltage after small pentrubation delta
            u_diff = [u_diff,data1-Fwd_u];
        end
        % Jacobian
        J = u_diff/delta;
        
        if isnan(det(J.'*J)) ==  1
            disp('The inverse of Jacobian matrix is NAN');
            
            break
        end
        % With regulation
        del_r = round(inv(J.'*J+lambda^2*eye(rn,rn))*J.'*(U_meas-Fwd_u),8);
        % without regulation
%         del_r = inv(J.'*J)*J.'*(U_meas-Fwd_u);
%         sum_resi = sqrt(sum((data1./U_meas -1).^2,'all')/length(data1));
%         Euc_norm = norm(U_meas-data1);
%         stand = std(U_meas-data1);
        error = sqrt(sum(abs(initial./all_r-1).^2,'all')/rn);
        i_error = abs(initial./all_r-1);
        indi_error = [indi_error;i_error];
        ave_error = [ave_error,error];
%         s_norm = [s_norm,Euc_norm];
%         s_resi = [s_resi,sum_resi];
%         s_std = [s_std,stand];
        if error > 100
            break
        end
        
        count = count +1;
%         if count > 100
%             disp('Results do not coverge');
%             break
%         end
    end
    
    accu_error = [accu_error, error];
end
% plot(hyper,accu_error)
% ylim([0,0.5]);
% for lambda = hyper(find(accu_error <=min(accu_error)))
%     % Count how many interations
%     count = 1;
%     % Run Resistance_matrix
%     initial = R*ones([1,rn]);
%     HR_sim = [];
%     for i = 1:vn
%         a = initial(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
%         HR_sim = [HR_sim;a];
%     end
%     VR_sim = [];
%     for i = 1:(vn-1)
%         b = initial(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
%         VR_sim = [VR_sim;b];
%     end
%     %% Modify afterwards  
%     Fwd_u = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
%     u_diff = [];
%     
%     % Generate jacobian
%     for i = 1:length(initial)
%         new = initial;
%         new(i) = initial(i) + delta;
%         HR_sim = [];
%         for i = 1:vn
%             a = new(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
%             HR_sim = [HR_sim;a];
%         end
%         VR_sim = [];
%         for i = 1:(vn-1)
%             b = new(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
%             VR_sim = [VR_sim;b];
%         end
%         data1 = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
%         % Obtain the voltage after small pentrubation delta
%         u_diff = [u_diff,data1-Fwd_u];
%     end
%     J = u_diff/delta;
% 
%     % Obtain the resistance updates.
%     del_r = inv(J.'*J+lambda^2*eye(rn,rn))*J.'*(U_meas-Fwd_u);
%     
% %     % L2 Norm
% %     Euc_norm = norm(U_meas-Fwd_u);
% %     s_norm = [Euc_norm];
% %     
% %     % Relative measurememnts error
% %     sum_resi = sqrt(sum((Fwd_u./U_meas -1).^2,'all')/length(data1));
% %     s_resi = [sum_resi];
% %     
% %     % Standard deviation
% %     stand = std(U_meas-Fwd_u);
% %     s_std = [stand];
%     
%     % Average/Macro error for the resistor.
%     error = sqrt(sum(abs(initial./all_r-1),'all')/rn);
%     ave_error = [error];
%     result = initial;
%     
%     % Individual/Microscopic error for each element
%     i_error = abs(initial./all_r-1);
%     indi_error = [i_error];
%     while error > 0.05
%         
%         % Complete the update
%         initial = initial+del_r';
%         result = [result;initial];
%         
%         % Transfer to a grid
%         HR_sim = [];
%         for i = 1:vn
%             a = initial(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
%             HR_sim = [HR_sim;a];
%         end
%         VR_sim = [];
%         for i = 1:(vn-1)
%             b = initial(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
%             VR_sim = [VR_sim;b];
%         end
%  
%         Fwd_u = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
% %         Euc_norm = norm(U_meas-Fwd_u);
%         
%         J = [];
%         u_diff = [];
%         % Generate jacobian
%         for i = 1:length(initial)
%             new = initial;
%             new(i) = initial(i) + delta;
%             HR_sim = [];
%             for i = 1:vn
%                 a = new(((i-1)*(2*hn-1)+1):((i-1)*(2*hn-1)+hn-1));
%                 HR_sim = [HR_sim;a];
%             end
%             VR_sim = [];
%             for i = 1:(vn-1)
%                 b = new(((i-1)*(2*hn-1)+hn):((i-1)*(2*hn-1)+hn+hn-1));
%                 VR_sim = [VR_sim;b];
%             end
%             data1 = kirchhoff_grid(hn,vn,HR_sim,VR_sim);
%             % Obtain the voltage after small pentrubation delta
%             u_diff = [u_diff,data1-Fwd_u];
%         end
%         % Jacobian
%         J = u_diff/delta;
%         
%         if isnan(det(J.'*J)) ==  1
%             disp('The inverse of Jacobian matrix is NAN');
%             
%             break
%         end
%         % With regulation
%         del_r = inv(J.'*J+lambda^2*eye(rn,rn))*J.'*(U_meas-Fwd_u);
%         % without regulation
% %         del_r = inv(J.'*J)*J.'*(U_meas-Fwd_u);
% %         sum_resi = sqrt(sum((data1./U_meas -1).^2,'all')/length(data1));
% %         Euc_norm = norm(U_meas-data1);
% %         stand = std(U_meas-data1);
%         error = sqrt(sum(abs(initial./all_r-1).^2,'all')/rn);
%         i_error = abs(initial./all_r-1);
%         indi_error = [indi_error;i_error];
%         ave_error = [ave_error,error];
% %         s_norm = [s_norm,Euc_norm];
% %         s_resi = [s_resi,sum_resi];
% %         s_std = [s_std,stand];
%         if error > 100
%             break
%         end
%         
%         count = count +1;
%         if count > 50
%             disp('Results do not coverge');
%             break
%         end
%     end
%     
%     accu_error = [accu_error, error];
% end