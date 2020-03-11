% 2D semi-analytical solution, using predictor-corrector method for 2
% subdomains split along the line gamma

% X: a vector of (nearby) particle positions
% [x0, y0]: the position of the particle at which the solution is centered
% r: vector of distances between the particle at x0 and particles at the
%     positions in X
% D1-D4: diffusion coefficient for subdomains 1-4 (correspond to standard
%     numbers assigned to quadrant in the cartesian plane)
% gammaX/Y: subdomain boundaries
% dt: length of the time step over which mass-transfers are conducted

% Xc: x-correction for transfers to opposite subdomain, calculated so as to
%     conserve mass

function solution = PrCo_2D_2x2(X, x0, y0, r, D1, D2, D3, D4, gammaX, gammaY, dt)

%     gaussian solution, as a function of separation distance and DD, the
%     diffusion coefficient, Eq. (32)
    c = @(r, DD) 1 / (4 * pi * DD * dt) .* exp(-(r'.^2) ./ (4 * DD * dt));

    if (x0 >= gammaX) && (y0 >= gammaY) % Quadrant 1

%         this corresponds to Eq. (35)
        Xc2 = x0 - sqrt(D2 / D1) * (x0 - gammaX);
        Xc3 = x0 - sqrt(D3 / D1) * (x0 - gammaX);
        Yc3 = y0 - sqrt(D3 / D1) * (y0 - gammaY);
        Yc4 = y0 - sqrt(D4 / D1) * (y0 - gammaY);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions. this corresponds to the pieces of Eq. (34)
        idx1 = find((X(:, 1) >= gammaX) & (X(:, 2) >= gammaY));
        idx2 = find((X(:, 1) < Xc2)    & (X(:, 2) >= gammaY));
        idx3 = find((X(:, 1) < Xc3)    & (X(:, 2) < Yc3));
        idx4 = find((X(:, 1) >= gammaX) & (X(:, 2) < Yc4));

        sol = zeros(size(X, 1), 4);

%         evaluate the pieces of Eq. (34)
        sol(idx1, 1) = c(r(idx1), D1);
        sol(idx2, 2) = c(r(idx2), D2);
        sol(idx3, 3) = c(r(idx3), D3);
        sol(idx4, 4) = c(r(idx4), D4);

%         This is Eq. (34)
        solution = sum(sol, 2);

    elseif (x0 < gammaX) && (y0 >= gammaY) % Quadrant 2

%         this corresponds to Eq. (35)
        Xc1 = x0 - sqrt(D1 / D2) * (x0 - gammaX);
        Yc3 = y0 - sqrt(D3 / D2) * (y0 - gammaY);
        Xc4 = x0 - sqrt(D4 / D2) * (x0 - gammaX);
        Yc4 = y0 - sqrt(D4 / D2) * (y0 - gammaY);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions. this corresponds to the pieces of Eq. (34)
        idx1 = find((X(:, 1) >= Xc1)  & (X(:, 2) >= gammaY));
        idx2 = find((X(:, 1) < gammaX) & (X(:, 2) >= gammaY));
        idx3 = find((X(:, 1) < gammaX) & (X(:, 2) < Yc3));
        idx4 = find((X(:, 1) >= Xc4)  & (X(:, 2) < Yc4));

        sol = zeros(size(X, 1), 4);

%         evaluate the pieces of Eq. (34)
        sol(idx1, 1) = c(r(idx1), D1);
        sol(idx2, 2) = c(r(idx2), D2);
        sol(idx3, 3) = c(r(idx3), D3);
        sol(idx4, 4) = c(r(idx4), D4);

%         This is Eq. (34)
        solution = sum(sol, 2);

    elseif (x0 < gammaX) && (y0 < gammaY) % Quadrant 3

%         this corresponds to Eq. (35)
        Xc1 = x0 - sqrt(D1 / D3) * (x0 - gammaX);
        Yc1 = y0 - sqrt(D1 / D3) * (y0 - gammaY);
        Yc2 = y0 - sqrt(D2 / D3) * (y0 - gammaY);
        Xc4 = x0 - sqrt(D4 / D3) * (x0 - gammaX);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions. this corresponds to the pieces of Eq. (34)
        idx1 = find((X(:, 1) >= Xc1)  & (X(:, 2) >= Yc1));
        idx2 = find((X(:, 1) < gammaX) & (X(:, 2) >= Yc2));
        idx3 = find((X(:, 1) < gammaX) & (X(:, 2) < gammaY));
        idx4 = find((X(:, 1) >= Xc4)  & (X(:, 2) < gammaY));

        sol = zeros(size(X, 1), 4);

%         evaluate the pieces of Eq. (34)
        sol(idx1, 1) = c(r(idx1), D1);
        sol(idx2, 2) = c(r(idx2), D2);
        sol(idx3, 3) = c(r(idx3), D3);
        sol(idx4, 4) = c(r(idx4), D4);

%         This is Eq. (34)
        solution = sum(sol, 2);

    elseif (x0 >= gammaX) && (y0 < gammaY) % Quadrant 4

%         this corresponds to Eq. (35)
        Yc1 = y0 - sqrt(D1 / D4) * (y0 - gammaY);
        Xc2 = x0 - sqrt(D2 / D4) * (x0 - gammaX);
        Yc2 = y0 - sqrt(D2 / D4) * (y0 - gammaY);
        Xc3 = x0 - sqrt(D3 / D4) * (x0 - gammaX);

%         find the particles in the support of the "keep" and "redistribute"
%         solutions. this corresponds to the pieces of Eq. (34)
        idx1 = find((X(:, 1) >= gammaX) & (X(:, 2) >= Yc1));
        idx2 = find((X(:, 1) < Xc2)    & (X(:, 2) >= Yc2));
        idx3 = find((X(:, 1) < Xc3)    & (X(:, 2) < gammaY));
        idx4 = find((X(:, 1) >= gammaX) & (X(:, 2) < gammaY));

        sol = zeros(size(X, 1), 4);

%         evaluate the pieces of Eq. (34)
        sol(idx1, 1) = c(r(idx1), D1);
        sol(idx2, 2) = c(r(idx2), D2);
        sol(idx3, 3) = c(r(idx3), D3);
        sol(idx4, 4) = c(r(idx4), D4);

%         This is Eq. (34)
        solution = sum(sol, 2);

    else
        fprintf('***Error***: your cases are not mutually exclusive\n')
    end

end

