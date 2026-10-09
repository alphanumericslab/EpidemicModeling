function coef = nonnegative_least_squares(x, y, max_iter, tolerance)
    % NONNEGATIVE_LEAST_SQUARES Solve NNLS by deterministic cyclic coordinate descent.
    % coef = nonnegative_least_squares(x,y,max_iter,tolerance)
    % x: observations-by-features; y: observation vector. Defaults: 10000, 1e-10.
    % Output: nonnegative feature column vector. Identical coordinate order,
    % initialization and stopping rule are used in Python. Zero columns are skipped.
    % Author: Reza Sameni | Emory University

    if nargin < 3
        max_iter = 10000;
    end

    if nargin < 4
        tolerance = 1e-10;
    end

    y = y(:);

    assert(size(x, 1) == numel(y) && all(isfinite(x(:))) && all(isfinite(y)), ...
        'Invalid regression data.');
    coef = zeros(size(x, 2), 1);
    residual = y;
    norms = sum(x.^2, 1);

    for k = 1:max_iter
        old = coef;

        for j = 1:numel(coef)

            if norms(j) == 0
                continue
            end

            delta = max(0, coef(j) + x(:, j)' * residual / norms(j)) - coef(j);
            coef(j) = coef(j) + delta;
            residual = residual - delta * x(:, j);
        end

        if max(abs(coef - old), [], 'all') <= tolerance * (1 + max(abs(coef), [], 'all'))
            break
        end
    end
end
