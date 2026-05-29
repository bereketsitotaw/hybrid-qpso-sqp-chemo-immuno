function solution = unpack_solution(z_opt, N)
%UNPACK_SOLUTION Extract state/control trajectories from NLP decision vector.

    idx = 1:(N + 1);
    solution.x = z_opt(idx);
    solution.y = z_opt(idx + (N + 1));
    solution.z = z_opt(idx + 2 * (N + 1));
    solution.u = z_opt(idx + 3 * (N + 1));
    solution.v = z_opt(idx + 4 * (N + 1));
    solution.T = z_opt(end);
    solution.time = linspace(0, solution.T, N + 1);
    solution.J = NaN;
end
