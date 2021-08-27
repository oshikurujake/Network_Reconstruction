h = 4;% Horizontav vength/evements number
v = 4;% Verticav vength/evements number
hn = h+1; % Horizontav nodes
vn = v+1; % Verticav nodes
nod = hn*vn;
R = 100;
% all_r = R*rand([1,40]);
% Hyperparameter




hyper = linspace(1e-2,1e-1,50);
accu_error = [];
for lambda = hyper
    HR = [all_r(1:4);all_r(10:13);all_r(19:22);all_r(28:31);all_r(37:40)];
    VR = [all_r(5:9);all_r(14:18);all_r(23:27);all_r(32:36)];
    count = 1;
    HR_o = HR;
    VR_o = VR;

    run Resistor_network_EIT_measuremnets
    % run Resistance_matrix
    U_meas = data1;
    first_iter = R*ones([1,40]);
    HR = [first_iter(1:4);first_iter(10:13);first_iter(19:22);first_iter(28:31);first_iter(37:40)];
    VR = [first_iter(5:9);first_iter(14:18);first_iter(23:27);first_iter(32:36)];
    HR_1 = HR;
    VR_1 = VR;
    run Resistor_network_EIT_measuremnets
    % run Resistance_matrix

    Fwd_u = data1;

    % Set a small pentrubation for each resistance 0.1 ohm
    %delta = mean(all_r,'all')/500;
    delta = 1;
    u_diff = [];


    clear J del_r
    for i = 1:length(first_iter)
        new = first_iter;
        new(i) = first_iter(i) + delta;
        HR = [new(1:4);new(10:13);new(19:22);new(28:31);new(37:40)];
        VR = [new(5:9);new(14:18);new(23:27);new(32:36)];
        run Resistor_network_EIT_measuremnets
    %     run Resistance_matrix
        u_diff = [u_diff,data1-Fwd_u];
    end
    J = u_diff/delta;
    % [r,j] = size(u_diff);

    % Generate jacobian

    del_r = inv(J.'*J+lambda^2*eye(40,40))*J.'*(U_meas-Fwd_u);

    Euc_norm = norm(U_meas-Fwd_u);
    s_norm = [Euc_norm];
    sum_resi = sqrt(sum((Fwd_u./U_meas -1).^2,'all')/208);
    s_resi = [sum_resi];
    stand = std(U_meas-Fwd_u);
    s_std = [stand];
    count = 1;
    error = sqrt(sum(abs(first_iter./all_r-1),'all')/40);
    ave_error = [error];
    result = first_iter;
    i_error = abs(first_iter./all_r-1);
    indi_error = [i_error];
    while error > 0.05
    %     del_r = meas_diff\J;
        first_iter = first_iter+del_r';
        result = [result;first_iter];
        HR = [first_iter(1:4);first_iter(10:13);first_iter(19:22);first_iter(28:31);first_iter(37:40)];
        VR = [first_iter(5:9);first_iter(14:18);first_iter(23:27);first_iter(32:36)];
        run Resistor_network_EIT_measuremnets
    %     run Resistance_matrix
        Fwd_u = data1;
        Euc_norm = norm(U_meas-Fwd_u);
        J = [];
        u_diff = [];
        for i = 1:length(first_iter)
            new = first_iter;
            new(i) = first_iter(i) + delta;
            HR = [new(1:4);new(10:13);new(19:22);new(28:31);new(37:40)];
            VR = [new(5:9);new(14:18);new(23:27);new(32:36)];
    %         run Resistance_matrix
            run Resistor_network_EIT_measuremnets
            u_diff = [u_diff,data1-Fwd_u];
        end
        J = u_diff/delta;
        % [r,j] = size(u_diff);

        % Calculate delta R
    %     del_r = (U_meas-Fwd_u)\J;
        if isnan(det(J.'*J)) ==  1
            disp('The inverse of Jacobian matrix is NAN');
            
            break
        end

        del_r = inv(J.'*J+lambda^2*eye(40,40))*J.'*(U_meas-Fwd_u);

        sum_resi = sqrt(sum((data1./U_meas -1).^2,'all')/208);
        Euc_norm = norm(U_meas-data1);
        stand = std(U_meas-data1);
        error = sqrt(sum((abs(first_iter./all_r)-1).^2,'all')/40);
        i_error = abs(first_iter./all_r-1);
        indi_error = [indi_error;i_error];
        ave_error = [ave_error,error];
        s_norm = [s_norm,Euc_norm];
        s_resi = [s_resi,sum_resi];
        s_std = [s_std,stand];
        if error > 100
        break
        end
        count = count+1;
        if count > 50
        break
        end
    end
    
    accu_error = [accu_error, error];
end

% 
% % consider first element in HR
% HR(1,1) = HR(1,1)+delta;
% 
% run Resistor_network_EIT_measuremnets
% 
% pentru = data1;
% % Calculate first column of Jacobian matrix
% Jcob = (pentru - Fwd_u)/delta;
