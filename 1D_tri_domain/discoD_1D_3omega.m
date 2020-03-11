%==========================================================================
% This script generates MTPT and RWPT solutions for the 1D discontinuous
% D(x) problem with 3 subdomains.
% This script will generate Figure 5 in:
%     "A Mass-transfer Particle-tracking Method for Simulating Transport
%     with Discontinuous Diffusion Coefficients,"
%     Advances in Water Resources, 2020.
%==========================================================================

% close all
clear variables

% position of point-mass source
x0 = -6.5;
% length of 1D domain
L = 5e1;
% number of particles for mass-transfer (MT) simulation
N = 5e3;
% number of particles for random-walk simulation
nRW = 1e6;
% total simulation time
maxT = 6;

% number of iterations to normalize using the Sinkhorn-Knopp Algorithm
nNorm = 1e3;

% number of members in ensemble (i.e., the number of D2 values to consider)
numEns = 3;

% subdomain boundary locations
gamma = [-10 ; 2];

% time step length and calculated number of time steps to take
dt = 1e-2;
nSteps = ceil(maxT / dt);

% diffusion coefficients (D2vec is length numEns)
D1 = 5;
D2vec = [2.5 1 0.5];
D3 = 0.05;

% locations of stationary MT particles
X = linspace(x0 - L / 2, x0 + L / 2, N)';

% dirac IC
[~, sourceX] = min(abs(X - x0));
X(sourceX) = x0;
mass = zeros(N, 1);
mass(sourceX) = 1 / (L / N);
mass0 = mass;

%% MT simulation

MTsolnVec = zeros(N, 3);

% D2 ensemble loop
for ens = 1 : numEns

%     choose the current D2 and build a vector of D's
    D2 = D2vec(ens);
    DVec = [D1 D2 D3];

%     calculate the max interaction distance for rangesearch()
    dist = 3 * sqrt(8 * max(DVec) * dt);

%     conduct the rangesearch to find nearby particles
    [idx, r] = rangesearch(X, X, dist, 'BucketSize', ceil(1e-2 * N));

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

%         this builds the weight matrix with the predictor-corrector solution
%         according to line 3 in Algorithm 2
        row(start : finish) = i;
        col(start : finish) = idx{i};
        [~, cutDex] = min(abs(X(i) - gamma));
        val(start : finish) = PrCo_1D_3omega(X(idx{i}), X(i), DVec(cutDex), DVec(cutDex + 1), gamma(cutDex), dt);
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
    fprintf('    Percent Mass Change = %.2f\n', 100 * abs(sumM0 - ensMsum) / sumM0)
end

%% RW Simulation

% start all the particles at the point-source location
RWsolnVec = x0 * ones(nRW, 3);

% D2 ensemble loop
for ens = 1 : numEns

    D2 = D2vec(ens);

    for i = 1 : nSteps
%         D value from where particles currently are
        Dpre = D1 .* (RWsolnVec(:, ens) <= gamma(1)) +...
               D2 .* (RWsolnVec(:, ens) > gamma(1) & RWsolnVec(:, ens) <= gamma(2))+...
               D3 .* (RWsolnVec(:, ens) > gamma(2));
        randy = randn(nRW, 1); % must use same random number for both steps
%         predictor step
        dX0 = sqrt(2 .* Dpre .* dt) .* randy;
        Xtemp = RWsolnVec(:, ens) + dX0;
%         D value from where particles land after predictor step
        Dcor = D1 .* (Xtemp <= gamma(1)) +...
               D2 .* (Xtemp > gamma(1) & Xtemp <= gamma(2)) +...
               D3 .* (Xtemp > gamma(2));
%         corrector step
        dX = sqrt(2 .* Dcor .* dt) .* randy;
        RWsolnVec(:, ens) = RWsolnVec(:, ens) + dX;
    end
end

nBins = 100;
bins = linspace(X(1), X(end), nBins);

% bin the particles to get concentrations
[RWPlot1, ~] = histcounts(RWsolnVec(:, 1), bins, 'Normalization', 'pdf');
[RWPlot2, ~] = histcounts(RWsolnVec(:, 2), bins, 'Normalization', 'pdf');
[RWPlot3, ~] = histcounts(RWsolnVec(:, 3), bins, 'Normalization', 'pdf');

% save the midpoints for plotting
bins = bins(2 : end) - (bins(2) - bins(1)) / 2;


%%  Plots

% calculated height to make plots pretty
height = 1.1 * max(max(MTsolnVec));

% get color order for plotting
color = lines(7);

% spatial mass plot for final time step
fig = 2;
figure(fig)
clf
fill([X(1) gamma(1) gamma(1) X(1) X(1)],height.*[0 0 1 1 0],[1 0 0],'FaceAlpha',0.05); grid on; hold on;
fill([gamma(2) gamma(1) gamma(1) gamma(2) gamma(2)],height.*[0 0 1 1 0],[0 1 0],'FaceAlpha',0.05);
fill([X(end) gamma(2) gamma(2) X(end) X(end)],height.*[0 0 1 1 0],[0 0 1],'FaceAlpha',0.05);
p2 = plot(X, MTsolnVec(:, 1), '-.+', 'color', color(1, :), 'linewidth', 2,...
          'MarkerSize', 10);
p3 = plot(bins, RWPlot1, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'g',...
          'linewidth', 3, 'MarkerSize', 12);
plot(X, MTsolnVec(:, 2), '-.+', 'color', color(2, :), 'linewidth', 2,...
          'MarkerSize', 10);
plot(bins, RWPlot2, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'g',...
          'linewidth', 3, 'MarkerSize', 12);
plot(X, MTsolnVec(:, 3), '-.+', 'color', color(3, :), 'linewidth', 2,...
          'MarkerSize', 10);
plot(bins, RWPlot3, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'g',...
          'linewidth', 3, 'MarkerSize', 12);
p4 = line([x0 x0], [0 height], 'color', 'b', 'linewidth', 5);
p5 = line([gamma(1) gamma(1)], [0 height], 'color', 'r', 'linewidth', 5);
line([gamma(2) gamma(2)], [0 height], 'color', 'r', 'linewidth', 5);
xlim([X(1) X(end)]); ylim([0 height]);
[~, hobj, ~, ~] = legend([p2 p3 p4 p5], {'\textbf{Mass-Transfer}',...
                         '\textbf{Random-Walk}', '\textbf{Source Location} {\boldmath($x_0$)}',...
                         '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'}, 'Interpreter',...
                         'latex', 'FontSize', 32, 'Location', 'northwest');
set(hobj, 'linewidth', 4);
rectangle('Position', [X(1) + 1.9, 0.58 * height, 5.9, 0.008], 'FaceColor', 'w',...
    'EdgeColor', 'k', 'LineWidth', 1)
text(X(1) + 2, 0.6 * height, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1, '%1.2f')],...
     'FontSize', 32, 'color', 'k', 'HorizontalAlignment', 'left',...
     'VerticalAlignment', 'middle', 'Interpreter', 'latex');
rectangle('Position', [X(N / 2) - 7.2, 0.84 * height, 5.9, 0.025], 'FaceColor', 'w',...
    'EdgeColor', 'k', 'LineWidth', 1)
for i = 1 : numEns
    text(X(N / 2) - 1.4, 0.95 * height - 0.008 * (i - 1),...
         ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i), '%1.2f')],...
          'FontSize', 32, 'color', color(i, :), 'HorizontalAlignment',...
          'right', 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
end
rectangle('Position', [X(end) - 7.8, 0.58 * height, 5.9, 0.008], 'FaceColor', 'w',...
    'EdgeColor', 'k', 'LineWidth', 1)
text(X(end) - 2, 0.6 * height, ['{\boldmath $D_{3}=\ $}\bf' num2str(D3, '%1.2f')],...
     'FontSize', 32, 'color', 'k', 'HorizontalAlignment', 'right',...
     'VerticalAlignment', 'middle', 'Interpreter', 'latex');
xlabel('\textbf{\emph{x}-position [L]}', 'Interpreter', 'latex', 'FontSize', 32)
ylabel('\textbf{Concentration [mol L{\boldmath$^{-1}$}]}', 'Interpreter', 'latex', 'FontSize', 32)
figure(fig)
set(gcf, 'Position', [0, 100, 1800, 900])
