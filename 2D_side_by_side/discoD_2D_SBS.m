%==========================================================================
% This script generates MTPT and RWPT solutions for the 2D discontinuous
% D(x) problem with 2 side-by-side subdomains.
% This script will generate Figures 6 and 7 in:
%     "A Mass-transfer Particle-tracking Method for Simulating Transport
%     with Discontinuous Diffusion Coefficients,"
%     Advances in Water Resources, 2020.
%==========================================================================

% close all
clear variables

% length of domain (x-direction)
L = 5e1;
% height of domain (y-direction)
H = 5e1;
% number of particles for MT simulation
Nx = 2e2 + 1;
Ny = 2e2 + 1;
N = Nx * Ny;
% number of particles for random-walk simulation
nRW = 1e7;
% total simulation time
maxT = 6;

% number of iterations to normalize using the Sinkhorn-Knopp Algorithm
nNorm = 1e3;

% number of members in ensemble
numEns = 3;

% position of point-mass source
x0 = -2;
y0 = 0;

% subdomain boundary locations
gamma = 2;

% time step length and calculated number of time steps to take
dt = 1e-1;
nSteps = ceil(maxT / dt);

% diffusion coefficients (D2vec is length numEns)
D1 = 5;
D2vec = [2.5 1.0 0.5];

% position vectors for stationary MT particles
X = linspace(-L / 2, L / 2, Nx)';
Y = linspace(-H / 2, H / 2, Ny)';
[mgX, mgY] = meshgrid(X, Y);
mtX = [reshape(mgX, N, 1) reshape(mgY, N, 1)];

% dirac IC
mass = zeros(N, 1);
idx = dsearchn(mtX, [x0 y0]);
mass(idx, :) = 1 / ((L * H) / N);
mass0 = mass;

%% MT Simulation

MTsolnVec = zeros(N, 3);

% D2 ensemble loop
for ens = 1 : numEns

%     set the current D2
    D2 = D2vec(ens);

%     calculate the max interaction distance for rangesearch()
    dist = 3 * sqrt(4 * max(D1, D2) * dt);

%     conduct the rangesearch to find nearby particles
    [idx, r] = rangesearch(mtX, mtX, dist, 'BucketSize', ceil(1e-2 * N));

%     determine how many particles are nearby and preallocate the vectors
%     to build the sparse weight matrix
    Nclose = sum(cellfun('length', idx));
    row = zeros(Nclose, 1);
    col = zeros(Nclose, 1);
    val = zeros(Nclose, 1);

%     calculate the entries of the weight matrix
    start = 1;
    for i = 1 : N
        num = length(idx{i});
        finish = start - 1 + num;
        row(start : finish) = i;
        col(start : finish) = idx{i};
        val(start : finish) = PrCo_2D_SBS(mtX(idx{i}, :), mtX(i, :), r{i}, D1, D2, gamma, dt);
        start = finish + 1;
    end

%     when N is super large the sparse methods and clearing large arrays saves memory
    clear idx r

%     create the sparse weight matrix
    Wmat = sparse(row, col, val);
    clear row col val

%     normalize via SK, with nNorm iterations, ending with a column
%     normalization
    Wmat = sinkhornKnoppCol(Wmat, 'MaxIter', nNorm);

%     build the transfer matrix according to Algorithm 2
    Tmat = Wmat;

    clear Wmat

%     initialize the mass vector
    mass = mass0;

%     conduct the mass transfers
    for i = 1 : nSteps
        mass = Tmat * mass;
    end

    clear Tmat

%     save the final mass vector for plotting
    MTsolnVec(:, ens) = mass;

end

%% print out final mass and percent mass change

sumM0 = sum(mass0);

fprintf('MT Initial Mass = %f\n', sumM0)

for i = 1 : numEns
    ensMsum = sum(MTsolnVec(:, i));
    fprintf('MT Final Mass (Ens = %i) = %f\n', i, ensMsum)
    fprintf('      Percent Mass Change = %.2f\n', 100 * abs(sumM0 - ensMsum) / sumM0)
end

%% RW Simulation

% start all the particles at the point-source location
RWsolnVec = zeros(nRW, 2, 3);
RWsolnVec(:, 1, :) = x0;
RWsolnVec(:, 2, :) = y0;

% D2 ensemble loop
for ens = 1 : 3

    D2 = D2vec(ens);

    for i = 1 : nSteps
%         D value from where particles currently are
        Dpre = D1 .* (RWsolnVec(:, 1, ens) <= gamma) + ...
               D2 .* (RWsolnVec(:, 1, ens) > gamma);
        randy = randn(nRW, 2); % must use same random number for both steps
%         predictor step
        dX0 = sqrt(2 .* Dpre .* dt) .* randy;
        Xtemp = RWsolnVec(:, :, ens) + dX0;
%         D value from where particles land after predictor step
        Dcor = D1 .* (Xtemp(:, 1) <= gamma) + ...
               D2 .* (Xtemp(:, 1) > gamma);
%         corrector step
        dX = sqrt(2 .* Dcor .* dt) .* randy;
        RWsolnVec(:, :, ens) = RWsolnVec(:, :, ens) + dX;
    end

end

%%

nBinsX = 80;
nBinsY = 80;
binsX = linspace(X(1), X(end), nBinsX);
binsY = linspace(Y(1), Y(end), nBinsY);

RWPlot = zeros(nBinsX - 1, nBinsY - 1, 3);

% bin the particles to get concentrations
[RWPlot(:, :, 1), ~, ~] = histcounts2(RWsolnVec(:, 1, 1), RWsolnVec(:, 2, 1), binsX, binsY, 'Normalization', 'pdf');
[RWPlot(:, :, 2), ~, ~] = histcounts2(RWsolnVec(:, 1, 2), RWsolnVec(:, 2, 2), binsX, binsY, 'Normalization', 'pdf');
[RWPlot(:, :, 3), ~, ~] = histcounts2(RWsolnVec(:, 1, 3), RWsolnVec(:, 2, 3), binsX, binsY, 'Normalization', 'pdf');

% generate a meshgrid at the midpoints for plotting
binsX = binsX(2 : end) - (binsX(2) - binsX(1)) / 2;
binsY = binsY(2 : end) - (binsY(2) - binsY(1)) / 2;
[BX, BY] = meshgrid(binsX, binsY);

%%  Plots

color = lines(7);

% for positioning the subplots to look nice
posCell = cell(6, 1);
posCell{1} = [0.03, 0.53, 0.30, 0.42];
posCell{2} = [0.34, 0.53, 0.30, 0.42];
posCell{3} = [0.67, 0.53, 0.30, 0.42];
posCell{4} = [0.03, 0.05, 0.30, 0.42];
posCell{5} = [0.34, 0.05, 0.30, 0.42];
posCell{6} = [0.67, 0.05, 0.30, 0.42];
%     [left bottom width height]

%% heatmap plot for final time step

figure(2)
clf
set(gcf, 'Position', [0, 100, 1400, 900])
for i = 1 : 3
    subplot('Position', posCell{i})
    pcolor(BX, BY, RWPlot(:, :, i)');
    hold on;
    shading flat;
    axis equal; axis tight;
    contour(BX, BY, RWPlot(:, :, i)', 3, 'w', 'linewidth', 2);
    p1 = plot(x0, y0, 'bx', 'linewidth', 2, 'MarkerSize', 20);
    p2 = plot([gamma gamma], [Y(1) Y(end)], '-r', 'linewidth', 2);
    xlim([X(1) X(end)]); ylim([Y(1) Y(end)]);
    legend([p1 p2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                     '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'},...
                     'Interpreter', 'latex', 'FontSize', 22, 'Location', 'northwest')
    text(X(1) + 1, 0 + H / 5, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1, '%1.1f')],...
                               'FontSize', 22, 'color', 'w', 'HorizontalAlignment',...
                               'left', 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end) - 1, 0 + H / 5, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i), '%1.1f')],...
                                 'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'right',...
                                 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    if i == 1
        ylabel('\textbf{\emph{y}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    end
    if i == 2
        title('\textbf{Random-Walk}', 'Interpreter', 'latex', 'FontSize', 22)
    end
    if i == 3
        cb = colorbar;
        ylabel(cb, '\textbf{Concentration [mol L{\boldmath$^{-2}$}]}',...
                   'Interpreter', 'latex', 'FontSize', 22)
    end
end
for i = 4 : 6
    MTPlot = reshape(MTsolnVec(:, i - 3), Ny, Nx);
    subplot('Position', posCell{i})
    imagesc(mtX(:, 1), mtX(:, 2), MTPlot);
    set(gca, 'YDir', 'normal');
    hold on;
    shading flat;
    axis equal; axis tight;
    contour(mgX, mgY, MTPlot, 3, 'w', 'linewidth', 2);
    p1 = plot(x0, y0, 'bx', 'linewidth', 2, 'MarkerSize', 20, 'linewidth', 2.5);
    p2 = plot([gamma gamma], [Y(1) Y(end)], '-r', 'linewidth', 2.5);
    xlim([X(1) X(end)]); ylim([Y(1) Y(end)]);
    [~, hobj, ~, ~] = legend([p1 p2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                                       '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'},...
                                       'Interpreter', 'latex', 'FontSize', 22,...
                                       'Location', 'northwest');
    set(hobj, 'linewidth', 2.5);
    text(X(1) + 1, 0 + H / 5, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1, '%1.1f')],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end) - 1, 0 + H / 5, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i - 3), '%1.1f')],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    xlabel('\textbf{\emph{x}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    if i == 4
        ylabel('\textbf{\emph{y}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    end
    if i == 5
        title('\textbf{Mass-Transfer}', 'Interpreter', 'latex', 'FontSize', 22)
    end
    if i == 6
        cb = colorbar;
        ylabel(cb, '\textbf{Concentration [mol L{\boldmath$^{-2}$}]}',...
               'Interpreter', 'latex', 'FontSize', 22)
    end
end

%% compare level sets for RW and MT simulations

posCell_1x3 = cell(6, 1);
posCell_1x3{1} = [0.03, 0.1, 0.30, 0.88];
posCell_1x3{2} = [0.36, 0.1, 0.30, 0.88];
posCell_1x3{3} = [0.69, 0.1, 0.30, 0.88];

figure(3)
clf
set(gcf, 'Position', [0, 100, 1400, 450])
for i = 1 : 3
    MTPlot = reshape(MTsolnVec(:, i), Ny, Nx);
    subplot('Position', posCell_1x3{i})
    contour(mgX, mgY, MTPlot, 3, 'EdgeColor', color(5, :), 'linewidth', 2); hold on
    contour(BX, BY, RWPlot(:, :, i)', 3, '--', 'EdgeColor', color(4, :), 'linewidth', 2);
    grid on
    hold on
    c1 = plot(NaN, '--', 'Color', color(4, :));
    c2 = plot(NaN, '-', 'Color', color(5, :));
    shading flat;
    axis equal; axis tight;
    p1 = plot(x0, y0, 'bx', 'linewidth', 2, 'MarkerSize', 20, 'linewidth', 2.5);
    p2 = plot([gamma gamma], [Y(1) Y(end)], '-r', 'linewidth', 2.5);
    xlim([X(1) X(end)]); ylim([Y(1) Y(end)]);
    [~, hobj, ~, ~] = legend([p1 p2 c1 c2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                           '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}',...
                           '\textbf{Random-Walk}', '\textbf{Mass-Transfer}'}, 'Interpreter',...
                           'latex', 'FontSize', 22, 'Location', 'northwest');
    set(hobj, 'linewidth', 2.5);
    text(X(1) + 1, 0 - H / 5, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1, '%1.1f')],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end) - 1, 0 - H / 5, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i), '%1.1f')],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    xlabel('\textbf{\emph{x}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    if i == 1
        ylabel('\textbf{\emph{y}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    end
end
