function [T, U, E] = check_energy(t, X, n, mass_bodies, G)

% Preallocate Kinetic and Potential Energy 
T = zeros(length(t),1);
U = zeros(length(t),1);

for k = 1:length(t)
    
    T_k = 0;
    U_k = 0;
    
    % Calculate Kinetic Energy 
    for i = 1:n
        v_i = X(k,3*n + 3*(i-1)+1 : 3*n + 3*i); % velocity 
        T_k = T_k + 0.5 * mass_bodies(i) * dot(v_i, v_i); % 0.5*mv^2
    end
    
    % Calculate Potential Energy 
    for i = 1:n
        r_i = X(k,3*(i-1)+1 : 3*i); % position of "ith" body
        
        for j = i+1:n
            r_j = X(k,3*(j-1)+1 : 3*j); % position of the other 3 bodies
            r_ij = norm(r_i - r_j); % difference in position
            
            % U = G*mi*mj/rij
            U_k = U_k + G * mass_bodies(i)*mass_bodies(j) / r_ij; 
        end
    end
    
    T(k) = T_k;
    U(k) = U_k;
end

E = T - U; % Total energy

end