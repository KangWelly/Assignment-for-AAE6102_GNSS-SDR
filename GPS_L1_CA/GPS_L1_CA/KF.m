x_est=[-2.416972450612282e+06;5.385212144527982e+06;2.408103440105958e+06;0;0];
n_states=5;
dtr     = pi/180;
P_est = diag([1e6, 1e6, 1e6, 1e2, 1e1]); % covariance matrix
Q = diag([1e-6, 1e-6, 1e-6, 1e-4*dt, 1e-6*dt]);
n_sat=5;
R = diag(ones(n_sat,1) * 5^2);
A = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 dt;
    0 0 0 0 1];
R = diag(ones(n_sat,1) * 5^2);
result={};
for epoch=1:size(navSolutions.PRN,2)
    result{epoch}=x_est;
    %n_sat=size(navSolutions.PRN(epoch),1);
    rho_true = zeros(n_sat, 1);
    X=navSolutions.satllitePosition{epoch};

    %---STEP 1
    x_pred = A * x_est;
    pos=x_pred;
    %---STEP 2
    P_pred = A * P_est * A' + Q;

    rho_meas=[];  H=[];
    for i=1:n_sat
        if epoch == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2;
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;

            %--- Correct satellite position (do to earth rotation) --------
            % Convert SV position at signal transmitting time to position
            % at signal receiving time. ECEF always changes with time as
            % earth rotates.
            Rot_X = e_r_corr(traveltime, X(:, i));

            %--- Find the elevation angel of the satellite ----------------
            [az(i), el(i), ~] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));

            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                    0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
            weight(i)=sin(el(i))^2;
        end % if iter == 1 ... ... else

        H(i, :) =  [ (-(Rot_X(1) - pos(1))) / norm(Rot_X - pos(1:3), 'fro') ...
            (-(Rot_X(2) - pos(2))) / norm(Rot_X - pos(1:3), 'fro') ...
            (-(Rot_X(3) - pos(3))) / norm(Rot_X - pos(1:3), 'fro') ...
            1 0];
        rho_meas(i,1)=navSolutions.correctedP(i,epoch);
    end
    H=H';
    %---STEP 3
    K = P_pred * H' / (H * P_pred * H' + R(1:n_sat,1:n_sat));
    z_pred = H * x_pred;
    %---STEP 4
    x_est = x_pred + K * (rho_meas - z_pred);
    %---STEP 5
    P_est = (eye(n_states) - K * H) * P_pred;
   
end
