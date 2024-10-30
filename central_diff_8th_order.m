function derivative = central_diff_8th_order(x, y)
    % Check if x and y have the same size
    if length(x) ~= length(y)
        error('x and y must have the same size');
    end

    % Ensure there are enough points for eighth-order central difference
    if length(x) < 9
        error('x and y must contain at least nine elements for eighth-order central differentiation');
    end

    % Initialize derivative vector
    n = length(y);
    derivative = zeros(size(y));

    % Uniform step size
    h = x(2) - x(1);

    % Eighth-order central differences for interior points
    for i = 5:(n-4)
        derivative(i) = (y(i-4) * 1/280 - y(i-3) * 4/105 + y(i-2) * 1/5 - y(i-1) * 4/5 + y(i+1) * 4/5 - y(i+2) * 1/5 + y(i+3) * 4/105 - y(i+4) * 1/280) / h;
    end

    % First-order differences at the very boundaries
    derivative(1) = (y(2) - y(1)) / h; % Forward difference for the first point
    derivative(2) = (y(3) - y(2)) / h; % Forward difference for the second point
    derivative(3) = (y(4) - y(3)) / h; % Forward difference for the third point
    derivative(4) = (y(5) - y(4)) / h; % Forward difference for the fourth point
    derivative(n-3) = (y(n-3) - y(n-4)) / h; % Backward difference for the fourth-to-last point
    derivative(n-2) = (y(n-2) - y(n-3)) / h; % Backward difference for the third-to-last point
    derivative(n-1) = (y(n-1) - y(n-2)) / h; % Backward difference for the second-to-last point
    derivative(n) = (y(n) - y(n-1)) / h; % Backward difference for the last point
end
