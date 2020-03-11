% 2D semi-analytical solution, using predictor-corrector method for 2
% subdomains split along the line gamma

% X: a vector of (nearby) particle positions
% x0: the x-position of the particle at which the solution is centered
% r: vector of distances between the particle at x0 and particles at the
%     positions in X
% D1: diffusion coefficient for subdomain 1 (on the left)
% D2: diffusion coefficient for subdomain 2 (on the right)
% gamma: subdomain boundary
% dt: length of the time step over which mass-transfers are conducted

% Xc: x-correction for transfers to opposite subdomain, calculated so as to
%     conserve mass

function solution = PrCo_2D_SBS(X, x0, r, D1, D2, gamma, dt)

    if x0(1) <= gamma

%         solution for x < gamma
        c1 = @(r) 1 / (4 * pi * D1 * dt) .* exp(-(r'.^2) ./ (4 * D1 * dt));

%         solution for x > gamma
        c2 = @(r) 1 / (4 * pi * D2 * dt) .* exp(-(r'.^2) ./ (4 * D2 * dt));

%         Eq. (27)
        Xc = x0(1) - sqrt(D2 / D1) * (x0(1) - gamma);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions, Eq. (32) for i = k, r
        negIndex = find(X(:, 1) <= gamma);
        posIndex = find(X(:, 1) > Xc);

        solution1 = zeros(size(X, 1), 1);
        solution2 = zeros(size(X, 1), 1);

%         Evaluate Eq. (32) for i = k, r, at the particles in their support
        solution1(negIndex) = c1(r(negIndex));
        solution2(posIndex) = c2(r(posIndex));

%         This is Eq. (26)
        solution = solution1 + solution2;

    else

%         solution for x < b
        c1 = @(r) 1 / (4 * pi * D1 * dt) .* exp(-(r'.^2) ./ (4 * D1 * dt));

%         solution for x > b
        c2 = @(r) 1 / (4 * pi * D2 * dt) .* exp(-(r'.^2) ./ (4 * D2 * dt));

%         Eq. (27)
        Xc = x0(1) - sqrt(D1 / D2) * (x0(1) - gamma);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions, Eq. (32) for i = k, r
        negIndex = find(X(:, 1) < Xc);
        posIndex = find(X(:, 1) >= gamma);

        solution1 = zeros(size(X, 1), 1);
        solution2 = zeros(size(X, 1), 1);

%         Evaluate Eq. (32) for i = k, r, at the particles in their support
        solution1(negIndex) = c1(r(negIndex));
        solution2(posIndex) = c2(r(posIndex));

%         This is Eq. (31)
        solution = solution1 + solution2;

    end

end

