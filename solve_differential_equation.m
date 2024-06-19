function concentrations_over_time = solve_differential_equation(master, alpha, rho, C0, dt, steps)
    C = C0;
    concentrations_over_time = zeros(length(C0), steps+1);
    concentrations_over_time(:, 1) = C0;
    for k = 2:steps+1
        dC = dCdt(C, master, alpha, rho);
        C = C + dt * dC;
        concentrations_over_time(:, k) = C;
    end
end
