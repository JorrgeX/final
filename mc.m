dimensions = [5, 10, 15, 20];
N = 1e6;

% Method One: Sampling in the Unit Cube C^d
disp('Method One: Sampling in the Unit Cube C^d');
for d = dimensions
    X = rand(N, d) - 0.5;
    sum_sq = sum(X.^2, 2);
    inside_ball = sum_sq <= 1;
    volume_estimate = sum(inside_ball) / N;
    fprintf('Estimated Volume for dimension %d: %f\n', d, volume_estimate);
end

% Method Two: Sampling in the Unit Ball B^d
disp(' ');
disp('Method Two: Sampling in the Unit Ball B^d');
for d = dimensions
    U = randn(N, d);
    norm_U = sqrt(sum(U.^2, 2));
    X = U ./ norm_U;
    r = rand(N, 1).^(1/d);
    X = X .* r;
    inside_cube = max(abs(X), [], 2) <= 0.5;
    vol_ball = pi^(d/2) / gamma(d/2 + 1);
    volume_estimate = sum(inside_cube) / N * vol_ball;
    fprintf('Estimated Volume for dimension %d: %f\n', d, volume_estimate);
end