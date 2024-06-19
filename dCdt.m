function dCdt_result = dCdt(C, master, alpha, rho)
    n = length(C);
    dCdt_result = zeros(size(C));
    for i = 1:n
        sum1 = alpha * C(i) * (1 - C(i));
        sum2 = 0;
        for j = 1:n
            sum2 = sum2 + master(i, j) * C(i) - master(i, j) * C(j);
        end
        dCdt_result(i) = sum1 - rho * sum2;
    end
end
