function [H,max_rel_drift] = check_angular_momentum(X,t,mass_bodies)

% Verify Angular Momentum Calculation
n = numel(mass_bodies);

% Preallocate
H = zeros(length(t),3);

% Loop through calculated position and velocities of the bodies
for k = 1:length(t)
    
    % Preallcoate H
    H_k = zeros(3,1);
    
    for i = 1:n
        
        r_i = X(k,3*(i-1)+1 : 3*i)'; % position at every "ith" body
        v_i = X(k,3*n + 3*(i-1)+1 : 3*n + 3*i)'; % velocity at every "ith" body
        
        H_k = H_k + mass_bodies(i) * cross(r_i, v_i); % H = mi * (ri x vi)
    end
    
    H(k,:) = H_k';
end



% Relative Drift
H0 = H(1,:); % initial angular momentum 

% norm of the difference between the inital and actual angular momentum
delta_H = vecnorm(H - H0, 2, 2); 
rel_drift = delta_H / norm(H0); % see how much angular momentum is changing

max_rel_drift = max(rel_drift); % find maximum drift

