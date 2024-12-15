% dimensions = [5, 10, 15, 20];
% N = 1e6;

% % Method One: Sampling in the Unit Cube C^d
% disp('Method One: Sampling in the Unit Cube C^d');
% for d = dimensions
%     X = rand(N, d) - 0.5;
%     sum_sq = sum(X.^2, 2);
%     inside_ball = sum_sq <= 1;
%     volume_estimate = sum(inside_ball) / N;
%     fprintf('Estimated Volume for dimension %d: %f\n', d, volume_estimate);
% end

% % Method Two: Sampling in the Unit Ball B^d
% disp(' ');
% disp('Method Two: Sampling in the Unit Ball B^d');
% for d = dimensions
%     U = randn(N, d);
%     norm_U = sqrt(sum(U.^2, 2));
%     X = U ./ norm_U;
%     r = rand(N, 1).^(1/d);
%     X = X .* r;
%     inside_cube = max(abs(X), [], 2) <= 0.5;
%     vol_ball = pi^(d/2) / gamma(d/2 + 1);
%     volume_estimate = sum(inside_cube) / N * vol_ball;
%     fprintf('Estimated Volume for dimension %d: %f\n', d, volume_estimate);
% end


ds = [5, 10, 15, 20];

for i = 1:length(ds)
    fprintf("Setting d = %d\n",ds(i))
    % part a
    vol = estimateVolumeCubeIntersectionBall(ds(i),1e6);
    fprintf("a) vol = %d\n",vol)
    

    % part b
    vol = estimateVolumeBallIntersectionCube(ds(i),1e6);
    fprintf("b) vol = %d\n",vol)
end


function vol = estimateVolumeBallIntersectionCube(d, n)
    % Generate n random points inside a d-dimensional unit ball
    xi = randn(n, d);
    radii = sqrt(sum(xi.^2, 2))*ones(1,d);
    % Scale points to be uniformly inside the ball
    scale = rand(n, 1).^(1/d) * ones(1,d);
    xi = (xi .* scale) ./ (radii);

    % Check if points are inside the unit cube
    ind  = find(max(abs(xi),[],2) <= 0.5);
    Na = length(ind);
    
    % Calculate the fraction of points inside the cube
    fractionInside = Na / n;

    % Volume of the d-dimensional unit ball
    ballVolume = pi^(d / 2) / (d / 2  * gamma(d / 2 ));
   
    % Approximate volume of the intersection
    vol = fractionInside * ballVolume;
end

function vol = estimateVolumeCubeIntersectionBall(d, numSamples)
    xi = rand(numSamples, d) - 0.5;
    
    ind = find(sum(xi.^2, 2) <= 1); %see if it's in the B^d
    Na = length(ind); 

    vol = Na / numSamples; % Volume of B^d \cap C^d
end