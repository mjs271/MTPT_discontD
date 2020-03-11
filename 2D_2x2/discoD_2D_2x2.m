%==========================================================================
% This script generates MTPT and RWPT solutions for the 2D discontinuous
% D(x) problem with 4 subdomains in a 2x2 tiling.
% This script will generate Figures 8 and 9 in:
%     "A Mass-transfer Particle-tracking Method for Simulating Transport
%     with Discontinuous Diffusion Coefficients,"
%     Advances in Water Resources, 2020.
%==========================================================================

% close all
clear variables

% length of domain (x-direction)
L = 8;
% height of domain (y-direction)
H = 8;
% number of particles for MT simulation
Nx = 2e2 + 1;
Ny = 2e2 + 1;
N = Nx * Ny;
% number of particles for random-walk simulation
nRW = 1e7;
% total simulation time
maxT = 3;

% number of iterations to normalize using the Sinkhorn-Knopp Algorithm
nNorm = 1e3;

% number of members in ensemble
numEns = 3;

% subdomain boundary locations
gammaX = 0;
gammaY = 0;

% position vectors for point-mass source
x0vec = [-1 -1 -1] .* 0.2;
y0vec = [-1 -1 -1] .* 0.2;

% time step length and calculated number of time steps to take
dt = 1e-1;
nSteps = ceil(maxT / dt);

% diffusion coefficients (vectors are length numEns)
DD = 0.3;
D1vec = [0.1 1.0 1.0] * DD;
D2vec = [0.5 0.1 1.0] * DD;
D3vec = [1.0 1.0 0.1] * DD;
D4vec = [0.3 1.0 1.0] * DD;

% position vectors for stationary MT particles
X = zeros(Nx, numEns);
Y = zeros(Ny, numEns);
mtX = zeros(N, 2, 3);
for i = 1 : numEns
    X(:, i) = linspace(x0vec(i) - L / 2, x0vec(i) + L / 2, Nx)';
    Y(:, i) = linspace(y0vec(i) - H / 2, y0vec(i) + H / 2, Ny)';
    [mgX, mgY] = meshgrid(X(:, i), Y(:, i));
    mtX(:, :, i) = [reshape(mgX, N, 1) reshape(mgY, N, 1)];
end
clear mgX mgY


%% MT Simulation

MTsolnVec = zeros(N, numEns);

% start the timer
tic

% ensemble loop
for ens = 1 : numEns

%     set current D's
    D1 = D1vec(ens);
    D2 = D2vec(ens);
    D3 = D3vec(ens);
    D4 = D4vec(ens);

%     set current point-source position and initialize mass vector
    x0 = x0vec(ens);
    y0 = y0vec(ens);
    mass = zeros(N, 1);
    idx0 = dsearchn(mtX(:, :, ens), [x0 y0]);
    mass(idx0, :) = 1 / ((L * H) / N);
    mass0 = mass;

%     calculate the max interaction distance for rangesearch()
    dist = 3 * sqrt(4 * max([D1 D2 D3 D4]) * dt);

%     conduct the rangesearch to find nearby particles
    [idx, r] = rangesearch(mtX(:, :, ens), mtX(:, :, ens), dist, 'BucketSize', ceil(1e-2 * N));

%     determine how many particles are nearby and preallocate the vectors
%     to build the sparse weight matrix
    Nclose = sum(cellfun('length', idx));
    row = zeros(Nclose, 1);
    col = zeros(Nclose, 1);
    val = zeros(Nclose, 1);

%     calculate the entries of the weight matrix
    start = 1;
    for i = 1 : N
        finish = start - 1 + length(idx{i});
        row(start : finish) = i;
        col(start : finish) = idx{i};
        val(start : finish) = PrCo_2D_2x2(mtX(idx{i}, :, ens), mtX(i, 1, ens),...
                                          mtX(i, 2, ens), r{i}, D1, D2, D3, D4,...
                                          gammaX, gammaY, dt);
        start = finish + 1;
    end

%     calculate the entries of the weight matrix
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

% print the run time
MT_time = toc

%% print out final mass and percent mass change

sumM0 = sum(mass0);

fprintf('MT Initial Mass = %f\n', sumM0)

for i = 1 : numEns
    ensMsum = sum(MTsolnVec(:, i));
    fprintf('MT Final Mass (Ens = %i) = %f\n', i, ensMsum)
    fprintf('    Percent Mass Change = %.2f\n', 100 * abs(sumM0 - ensMsum) / sumM0)
end

%% RW Simulation

RWsolnVec = zeros(nRW, 2, numEns);

% use a smaller time step to get less noisy results
dt = 1e-2;
nSteps = ceil(maxT / dt);

% start the timer
tic

for ens = 1 : numEns

%     set current point-source position
    D1 = D1vec(ens);
    D2 = D2vec(ens);
    D3 = D3vec(ens);
    D4 = D4vec(ens);

%     set current point-source position
    x0 = x0vec(ens);
    y0 = y0vec(ens);

% start all the particles at the point-source location
    RWsolnVec(:, 1, ens) = x0;
    RWsolnVec(:, 2, ens) = y0;
    RWVec = RWsolnVec(:, :, ens);

    for i = 1 : nSteps
%         D value from where particles currently are
        Dpre = D1 .* (RWVec(:, 1) >= gammaX) .* (RWVec(:, 2) >= gammaY) + ... % Upper Right (Quadrant 1)
               D2 .* (RWVec(:, 1) < gammaX)  .* (RWVec(:, 2) >= gammaY) + ... % Upper Left (Quadrant 2)
               D3 .* (RWVec(:, 1) < gammaX)  .* (RWVec(:, 2) < gammaY) + ... % Lower Left (Quadrant 3)
               D4 .* (RWVec(:, 1) >= gammaX) .* (RWVec(:, 2) < gammaY); % Lower Right (Quadrant 4)
        randy = randn(nRW, 2); % must use same random number for both steps
%         predictor step
        dX0 = sqrt(2 .* Dpre .* dt) .* randy;
        Xtemp = RWVec(:, :) + dX0;
%         D value from where particles land after predictor step
        Dcor = D1 .* (Xtemp(:, 1) >= gammaX) .* (Xtemp(:, 2) >= gammaY) + ... % Upper Right
               D2 .* (Xtemp(:, 1) < gammaX)  .* (Xtemp(:, 2) >= gammaY) + ... % Upper Left
               D3 .* (Xtemp(:, 1) < gammaX)  .* (Xtemp(:, 2) < gammaY) + ... % Lower Left
               D4 .* (Xtemp(:, 1) >= gammaX) .* (Xtemp(:, 2) < gammaY); % Lower Right
%         corrector step
           dX = sqrt(2 .* Dcor .* dt) .* randy;
        RWVec(:, :) = RWVec(:, :) + dX;
    end
    RWsolnVec(:, :, ens) = RWVec;
end

% print the run time
RW_runTime = toc

%%

% start the timer
tic

nBinsX = 101;
nBinsY = 101;
RWPlot = zeros(nBinsX - 1, nBinsY - 1, numEns);

BX = zeros(nBinsX - 1, nBinsY - 1, 3);
BY = zeros(nBinsX - 1, nBinsY - 1, 3);

% bin the particles to get concentrations
for i = 1 : numEns
    binsX = linspace(X(1, i), X(end, i), nBinsX);
    binsY = linspace(Y(1, i), Y(end, i), nBinsY);
    [RWPlot(:, :, i), ~, ~] = histcounts2(RWsolnVec(:, 1, i), RWsolnVec(:, 2, i),...
                                          binsX, binsY, 'Normalization', 'pdf');
    binsX = binsX(2 : end) - (binsX(2) - binsX(1)) / 2;
    binsY = binsY(2 : end) - (binsY(2) - binsY(1)) / 2;
    [BX(:, :, i), BY(:, :, i)] = meshgrid(binsX, binsY);
end

% print the run time
RW_binTime = toc

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

fig = 4;
figure(fig)
clf
for i = 1 : numEns
    subplot('Position', posCell{i})
    pcolor(BX(:, :, i), BY(:, :, i), RWPlot(:, :, i)');
    hold on;
    shading flat;
    axis equal; axis tight;
    contour(BX(:, :, i), BY(:, :, i), RWPlot(:, :, i)', 3, 'w', 'linewidth', 2);
    p1 = plot(x0vec(i), y0vec(i), 'bx', 'linewidth', 2, 'MarkerSize', 20, 'linewidth', 2.5);
    p2 = plot([gammaX gammaX], [Y(1, i) Y(end, i)], '-r', 'linewidth', 2.5);
    plot([X(1, i) X(end, i)], [gammaY gammaY], '-r', 'linewidth', 2);
    xlim([X(1, i) X(end, i)]); ylim([Y(1, i) Y(end, i)]);
    [~, hobj, ~, ~] = legend([p1 p2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                             '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'},...
                             'Interpreter', 'latex', 'FontSize', 22, 'Location', 'northwest');
    set(hobj, 'linewidth', 2.5);
    text(X(end, i) - 0.25, 0 + H / 5, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1vec(i))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontWeight', 'bold');
    text(X(1, i) + 0.25, 0 + H / 5, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(1, i) + 0.25, 0 - H / 5, ['{\boldmath $D_{3}=\ $}\bf' num2str(D3vec(i))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end, i) - 0.25, 0 - H / 5, ['{\boldmath $D_{4}=\ $}\bf' num2str(D4vec(i))],...
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
for i = numEns + 1 : 2 * numEns
    MTPlot = reshape(MTsolnVec(:, i - numEns), Ny, Nx);
    [mgX, mgY] = meshgrid(X(:, i - numEns), Y(:, i - numEns));
    subplot('Position', posCell{i})
    imagesc(mtX(:, 1, i - numEns), mtX(:, 2, i - numEns), MTPlot);
    set(gca, 'YDir', 'normal');
    hold on;
    shading flat;
    axis equal; axis tight;
    contour(mgX, mgY, MTPlot, 3, 'w', 'linewidth', 2);
    p1 = plot(x0vec(i - numEns), y0vec(i - numEns), 'bx', 'linewidth', 2.5, 'MarkerSize', 20);
    p2 = plot([gammaX gammaX], [Y(1, i - numEns) Y(end, i - numEns)], '-r', 'linewidth', 2.5);
    plot([X(1, i - numEns) X(end, i - numEns)], [gammaY gammaY], '-r', 'linewidth', 2);
    xlim([X(1, i - numEns) X(end, i - numEns)]); ylim([Y(1, i - numEns) Y(end, i - numEns)]);
    [~, hobj, ~, ~] = legend([p1 p2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                             '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'},...
                             'Interpreter', 'latex', 'FontSize', 22, 'Location', 'northwest');
    set(hobj, 'linewidth', 2.5);
    text(X(end, i - numEns) - 0.25, 0 + H / 5, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1vec(i - numEns))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontWeight', 'bold');
    text(X(1, i - numEns) + 0.25, 0 + H / 5, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i - numEns))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(1, i - numEns) + 0.25, 0 - H / 5, ['{\boldmath $D_{3}=\ $}\bf' num2str(D3vec(i - numEns))],...
         'FontSize', 22, 'color', 'w', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end, i - numEns) - 0.25, 0 - H / 5, ['{\boldmath $D_{4}=\ $}\bf' num2str(D4vec(i - numEns))],...
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
figure(fig)
set(gcf, 'Position', [0, 100, 1400, 900])

%% compare level sets for RW and MT simulations

posCell_1x3 = cell(6, 1);
posCell_1x3{1} = [0.03, 0.1, 0.30, 0.88];
posCell_1x3{2} = [0.36, 0.1, 0.30, 0.88];
posCell_1x3{3} = [0.69, 0.1, 0.30, 0.88];

fig = 6;
figure(fig)
clf
for i = 1 : numEns
    MTPlot = reshape(MTsolnVec(:, i), Ny, Nx);
    [mgX, mgY] = meshgrid(X(:, i), Y(:, i));
    subplot('Position', posCell_1x3{i})
    contour(mgX, mgY, MTPlot, 3, 'EdgeColor', color(5, :), 'linewidth', 2);
    hold on
    grid on
    contour(BX(:, :, i), BY(:, :, i), RWPlot(:, :, i)', 3, '--',...
            'EdgeColor', color(4, :), 'linewidth', 2);
    c1 = plot(NaN, '--', 'Color', color(4, :));
    c2 = plot(NaN, '-', 'Color', color(5, :));
    shading flat;
    axis equal; axis tight;
    p1 = plot(x0vec(i), y0vec(i), 'bx', 'linewidth', 2, 'MarkerSize', 20, 'linewidth', 2.5);
    p2 = plot([gammaX gammaX], [Y(1, i) Y(end, i)], '-r', 'linewidth', 2.5);
    plot([X(1, i) X(end, i)], [gammaY gammaY], '-r', 'linewidth', 2);
    xlim([X(1, i) X(end, i)]); ylim([Y(1, i) Y(end, i)]);
    [~, hobj, ~, ~] = legend([p1 p2 c1 c2], {'\textbf{Source Location} {\boldmath$(x_0)$}',...
                           '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}',...
                           '\textbf{Random-Walk}', '\textbf{Mass-Transfer}'}, 'Interpreter',...
                           'latex', 'FontSize', 22, 'Location', 'northwest');
    set(hobj, 'linewidth', 2.5);
    text(X(end, i) - 0.25, 0 + H / 8, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1vec(i))],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(1, i) + 0.25, 0 + H / 8, ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i))],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(1, i) + 0.25, 0 - H / 5, ['{\boldmath $D_{3}=\ $}\bf' num2str(D3vec(i))],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'left',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    text(X(end, i) - 0.25, 0 - H / 5, ['{\boldmath $D_{4}=\ $}\bf' num2str(D4vec(i))],...
         'FontSize', 22, 'color', 'k', 'HorizontalAlignment', 'right',...
         'VerticalAlignment', 'middle', 'Interpreter', 'latex');
    xlabel('\textbf{\emph{x}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    if i == 1
        ylabel('\textbf{\emph{y}-position [L]}', 'Interpreter', 'latex', 'FontSize', 22)
    end
end
figure(fig)
set(gcf, 'Position', [0, 100, 1400, 450])
