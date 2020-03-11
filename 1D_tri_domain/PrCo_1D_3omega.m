% 1D semi-analytical solution, using predictor-corrector method

% X: a vector of (nearby) particle positions
% x0: the position of the particle at which the solution is centered
% D1: diffusion coefficient for subdomain 1 (on the left)
% D2: diffusion coefficient for subdomain 2 (on the right)
% gamma: subdomain boundary
% dt: length of the time step over which mass-transfers are conducted

% Xc: x-correction for transfers to opposite subdomain, calculated so as to
%     conserve mass

function solution = PrCo_1D_3omega(X, x0, D1, D2, gamma, dt)

%     Solution for source on left side of domain--i.e., Eq. (26)
    if x0 <= gamma

%         solution for x < gamma
        c1 = @(x) 1 / sqrt(4 * pi * D1 * dt) * exp(-(x - x0).^2/(4 * D1 * dt));

%         solution for x > Xc
        c2 = @(x) 1 / sqrt(4 * pi * D2 * dt) * exp(-(x - x0).^2/(4 * D2 * dt));

%         Eq. (27)
        Xc = x0 - sqrt(D2 / D1) * (x0 - gamma);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions, Eq's (21) and (22)
        negIndex = find(X <= gamma);
        posIndex = find(X > Xc);

        solution1 = zeros(size(X));
        solution2 = zeros(size(X));

%         Evaluate Eq's (21) and (22) at the particles in their support
        solution1(negIndex) = c1(X(negIndex));
        solution2(posIndex) = c2(X(posIndex));

%         This is Eq. (26)
        solution = solution1 + solution2;

%     Solution for source on right side of domain--i.e., Eq. (20)
    else

%         solution for x < Xc
        c1 = @(x) 1 / sqrt(4 * pi * D1 * dt) * exp(-(x - x0).^2/(4 * D1 * dt));

%         solution for x > gamma
        c2 = @(x) 1 / sqrt(4 * pi * D2 * dt) * exp(-(x - x0).^2/(4 * D2 * dt));

%         Eq. (25)
        Xc = x0 - sqrt(D1 / D2) * (x0 - gamma);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions, Eq's (21) and (22)
        negIndex = find(X < Xc);
        posIndex = find(X >= gamma);

        solution1 = zeros(size(X));
        solution2 = zeros(size(X));

%         Evaluate Eq's (21) and (22) at the particles in their support
        solution1(negIndex) = c1(X(negIndex));
        solution2(posIndex) = c2(X(posIndex));

%         This is Eq. (20)
        solution = solution1 + solution2;

    end
end

