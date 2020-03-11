%==========================================================================
% This script conducts a convergence analysis in terms of particle number,
% N, for the discontinuous D(x) MTPT method employing a semi-analytical solution.
% This script, will generate Figure 11 in:
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
% total simulation time
maxT = 6;

% number of iterations to normalize using the Sinkhorn-Knopp Algorithm
nNorm = 1e3;

% number of members in ensemble (i.e., the number of D2 values to consider)
numDEns = 3;

% subdomain boundary location
gamma = -1;

% time step length and calculated number of time steps to take
dt = 1e-2;
nSteps = ceil(maxT / dt);

% number of refinements in N and smalles value of N to begin with
numN = 9;
NvecStart = 1e2;

% 2x refinements in N
Nvec = NvecStart * 2.^[0 : numN - 1]

% diffusion coefficients (D2vec is length numDEns)
D1 = 5;
D2vec = [2.5 0.5 0.05];

%% MT simulation

% dimension of array is: # of error metrics x N ensemble x D2 ensemble
errVec = zeros(4, numN, 3);

% particle number ensemble loop
for nEns = 1 : numN

    % number of particles
    N = Nvec(nEns);

    % locations of stationary MT particles
    X = linspace(x0 - L / 2, x0 + L / 2, N)';

    % calculate analytical solution, using Eqs. (12)-(16)
    analyticalSoln = zeros(N, 3);
    for i = 1 : 3
        analyticalSoln(:, i) = MT_CnJ_arb_dCut(X, x0, D1, D2vec(i), gamma, maxT);
    end

    % calculated height to make plots pretty
    height = 1.1 * max(max(analyticalSoln));

    % dirac IC
    [~, sourceX] = min(abs(X - x0));
    X(sourceX) = x0;
    mass = zeros(N, 1);
    mass(sourceX) = 1 / (L / N);
    mass0 = mass;

    MTsolnVec = zeros(N, 3);

%     D2 ensemble loop
    for Dens = 1 : numDEns

%         choose the current D2
        D2 = D2vec(Dens);

%         calculate the max interaction distance for rangesearch()
        dist = 3 * sqrt(4 * max([D1 D2]) * dt);

%         conduct the rangesearch to find nearby particles
        [idx, r] = rangesearch(X, X, dist, 'BucketSize', ceil(1e-2 * N));

%         determine how many particles are nearby and preallocate the vectors
%         to build the sparse weight matrix
        Nclose = sum(cellfun('length', idx));
        row = zeros(Nclose, 1);
        col = zeros(Nclose, 1);
        val = zeros(Nclose, 1);

%         calculate the entries of the weight matrix
        start = 1;
        for i = 1 : N
            finish = start - 1 + length(idx{i});

%             this builds the weight matrix with the predictor-corrector solution
%             according to line 3 in Algorithm 2
            row(start : finish) = i;
            col(start : finish) = idx{i};
            val(start : finish) = PrCo_1D_2omega(X(idx{i}), X(i), D1, D2, gamma, dt);
            start = finish + 1;
        end

%         when N is super large the sparse methods and clearing large arrays saves memory
        clear idx r

%         create the sparse weight matrix
        Wmat = sparse(row, col, val);
        clear row col val

%         normalize via SK, with nNorm iterations, ending with a column
%         normalization
        Wmat = sinkhornKnoppCol(Wmat, 'MaxIter', nNorm);

%         build the transfer matrix according to Algorithm 2
        Tmat = Wmat;

%         initialize the mass vector
        mass = mass0;

%         conduct the mass transfers
        for i = 1 : nSteps
            mass = Tmat * mass;
        end

%         save the final mass vector for plotting
        MTsolnVec(:, Dens) = mass;

%         calculate the error in RMSE, then inf, 2, and 1 norms
        errVec(1, nEns, Dens) = sqrt(mean((MTsolnVec(:, Dens) - analyticalSoln(:, Dens)).^2));
        errVec(2, nEns, Dens) = norm(MTsolnVec(:, Dens) - analyticalSoln(:, Dens), inf);
        errVec(3, nEns, Dens) = norm(MTsolnVec(:, Dens) - analyticalSoln(:, Dens), 2);
        errVec(4, nEns, Dens) = norm(MTsolnVec(:, Dens) - analyticalSoln(:, Dens), 1);

    end

end

%% print out final mass and percent mass change

sumM0 = sum(mass0);

fprintf('MT Initial Mass = %f\n', sumM0)

for i = 1 : numDEns
    ensMsum = sum(MTsolnVec(:, i));
    fprintf('MT Final Mass (Ens = %i) = %f\n', i, ensMsum)
    fprintf('    Percent Mass Change = %.2f\n', 100 * abs(sumM0 - ensMsum) / sumM0)
end


%%  Final time solution plots

% get color order for plotting
color = lines(7);

% spatial mass plot for final time step
fig = 2;
figure(fig)
clf
fill([X(1) gamma gamma X(1) X(1)],height.*[0 0 1 1 0],[0 0.25 1],'FaceAlpha',0.05); grid on; hold on;
fill([X(end) gamma gamma X(end) X(end)],height.*[0 0 1 1 0],[1 0.25 0],'FaceAlpha',0.05);

p1 = plot(X, analyticalSoln(:, 1), 'k', 'linewidth', 17);
p2 = plot(X(:), MTsolnVec(:, 1), '-.+', 'color', color(1, :), 'linewidth', 2,...
          'MarkerSize', 10);

plot(X, analyticalSoln(:, 2), 'k', 'linewidth', 17)
plot(X(:), MTsolnVec(:, 2), '-.+', 'color', color(2, :), 'linewidth', 2,...
          'MarkerSize', 10);

plot(X, analyticalSoln(:, 3), 'k', 'linewidth', 17)
plot(X(:), MTsolnVec(:, 3), '-.+', 'color', color(3, :), 'linewidth', 2,...
          'MarkerSize', 10);

p4 = line([x0 x0], [0 height], 'color', 'b', 'linewidth', 5);
p5 = line([gamma gamma], [0 height], 'color', 'r', 'linewidth', 5);
xlim([X(1) X(end)]); ylim([0 height]);
[~, hobj, ~, ~] = legend([p1 p2 p4 p5], {'\textbf{Analytic Solution}',...
                         '\textbf{Mass-Transfer}', '\textbf{Source Location} {\boldmath$(x_0)$}',...
                         '\textbf{Subdomain Boundary} {\boldmath($\gamma$)}'}, 'Interpreter',...
                         'latex', 'FontSize', 32, 'Location', 'northwest');
set(hobj, 'linewidth', 4);
rectangle('Position', [X(end) - 6.8, 0.48 * height, 5.9, 0.0116], 'FaceColor', 'w',...
    'EdgeColor', 'k', 'LineWidth', 1)
for i = 1 : 3
    text(X(end) - 1, 0.6 * height - 0.004 * (i - 1), ['{\boldmath $D_{2}=\ $}\bf' num2str(D2vec(i), '%1.2f')],...
         'FontSize', 32, 'color', color(i, :), 'HorizontalAlignment',...
         'right', 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
end
rectangle('Position', [X(1) + 0.9, 0.58 * height, 5.9, 0.0036], 'FaceColor', 'w',...
    'EdgeColor', 'k', 'LineWidth', 1)
text(X(1) + 1, 0.6 * height, ['{\boldmath $D_{1}=\ $}\bf' num2str(D1, '%1.2f')], 'FontSize', 32, 'color', 'k',...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
xlabel('\textbf{\emph{x}-position [L]}', 'Interpreter', 'latex', 'FontSize', 32)
ylabel('\textbf{Concentration [mol L{\boldmath$^{-1}$}]}', 'Interpreter', 'latex', 'FontSize', 32)
figure(fig)
set(gcf, 'Position', [0, 100, 1800, 900])

%% Convergence plots

fig = 27;
figure(fig)
clf
for i = 1 : Dens
    subplot(1, 3, i)
    loglog(Nvec, errVec(1, :, i), '-o', 'LineWidth', 5.5)
    hold on
    loglog(Nvec, errVec(2, :, i), '-o', 'LineWidth', 5.5)
    xlabel('\boldmath{$N$}','Interpreter','latex', 'FontSize', 22)
    ylabel('\textbf{Error}','Interpreter','latex', 'FontSize', 22)
    legend({'\textbf{RMSE}', '\boldmath{$\ell^\infty$} \textbf{Norm}'},...
            'Interpreter', 'latex', 'FontSize', 22, 'Location', 'northeast')
    title(['\boldmath{$D_2 = $} \bf', num2str(D2vec(i))], 'Interpreter', 'latex', 'FontSize', 22)
end
figure(fig)
set(gcf, 'Position', [0, 100, 1400, 450])
