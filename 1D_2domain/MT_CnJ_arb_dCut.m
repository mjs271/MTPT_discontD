% analytical solution for composite diffusion, adapted from Carslaw and
% Jaeger, "Conduction of Heat in Solids", p. 363

% Inputs/Assumptions:
%     1. X, the 1D position vector, must be in ascending order and contain zero
%     2. x0 is the x-location of the Dirac pulse source
%     3. D1 is diffusion coefficient for x < gamma
%     4. D2 is diffusion coefficient for x > gamma
%     5. maxT is the final simulation time of interest
%     6. the initial condition is a pulse of mass 1, centered at x = x0
%     7. gamma is the location of the interface between the subdomains with
%     distinct diffusion coefficients

function solution = MT_CnJ_arb_dCut(X, x0, D1, D2, gamma, maxT)

% Eq. (12)
if x0 >= gamma

%     solution for x < 0
    c1 = @(x) (D2 * D1 * D2^(-1 / 2)) /...
              ((D2 * sqrt(D1) + D1 * sqrt(D2)) * sqrt(pi * D1 * maxT)) *...
              exp(-(x - gamma - (x0 - gamma) * sqrt(D1 / D2)).^2 / (4 * D1 * maxT));

%     solution for x > 0
    c2 = @(x) 1 / (2 * sqrt(pi * D2 * maxT)) *...
              exp(-(x - x0).^2 / (4 * D2 * maxT)) + ...
              (D2 * sqrt(D1) - D1 * sqrt(D2)) /...
              (2 * (D2 * sqrt(D1) + D1 * sqrt(D2)) * sqrt(pi * D2 * maxT)) *...
              exp(-(x + x0 - 2 * gamma).^2 / (4 * D2 * maxT));

    negIndex = find(X <= gamma);
    posIndex = find(X > gamma);      

    solution = zeros(size(X));

    solution(negIndex) = c1(X(negIndex));
    solution(posIndex) = c2(X(posIndex));

% Eq. (16)
else
    
    temp = D1;
    D1 = D2;
    D2 = temp;
    
%     solution for x > 0
    c1 = @(x) (D2 * D1 * D2^(-1 / 2)) /...
              ((D2 * sqrt(D1) + D1 * sqrt(D2)) * sqrt(pi * D1 * maxT)) *...
              exp(-(x - gamma - (x0 - gamma) * sqrt(D1 / D2)).^2 / (4 * D1 * maxT));

%     solution for x < 0
    c2 = @(x) 1 / (2 * sqrt(pi * D2 * maxT)) *...
              exp(-(x - x0).^2 / (4 * D2 * maxT)) + ...
              (D2 * sqrt(D1) - D1 * sqrt(D2)) /...
              (2 * (D2 * sqrt(D1) + D1 * sqrt(D2)) * sqrt(pi * D2 * maxT)) *...
              exp(-(x + x0 - 2 * gamma).^2 / (4 * D2 * maxT));

    negIndex = find(X <= gamma);
    posIndex = find(X > gamma);

    solution = zeros(size(X));

    solution(negIndex) = c2(X(negIndex));
    solution(posIndex) = c1(X(posIndex));

end

