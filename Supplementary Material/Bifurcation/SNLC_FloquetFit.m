% Fit Floquet multipliers near a saddle-node on limit cycle (SNLC)
% using the square-root scaling model:
%
%               lambda(mu) = 1 - A * sqrt(mu - mu_c)
%
% Compatible with MATLAB R2018a+. It preferentially uses `lsqcurvefit`
% (Optimization Toolbox); if unavailable, the fallback is `fminsearch`.

%% ========== 1) Load Data ==========

% Option A: load from CSV (two columns: mu and lambda). Header allowed.
% dataFile = 'floquet_data.csv';
% mu = []; lambdaVals = [];
%
% if exist(dataFile,'file')
%     T = readtable(dataFile);
%     if width(T) < 2
%         error('CSV must contain at least two columns: first mu, second lambda.');
%     end
%     mu = T{:,1};
%     lambdaVals = T{:,2};
% end

% Option B (recommended for quick use): directly assign vectors.
% Remove the comment and fill in the vectors:
mu = K_keep;
lambdaVals = mu_keep;

if isempty(mu)
    error(['No data found. Please supply the CSV file "' dataFile '" ', ...
           'or directly assign vectors mu and lambdaVals inside this script.']);
end

mu = mu(:); 
lambdaVals = lambdaVals(:);

% Sort by mu
[mu, idx] = sort(mu);
lambdaVals = lambdaVals(idx);

%% ========== 2) Specify the Known Critical Parameter mu_c (K_c) ==========

Kc = 105.9;   % <<<=== Replace with the known critical value K_c ===>>>

fprintf('Performing fit using fixed Kc = %.6g.\n', Kc);

%% ========== 3) Fit Parameter A Only ==========

% Model: lambda(mu) = 1 - A * sqrt(mu - Kc)
model = @(A, m) 1 - A * sqrt(max(m - Kc, 0));

% Initial guess for A
A0 = max((1 - min(lambdaVals)) / sqrt(max(max(mu) - Kc, 1e-6)), 1e-4);

% Bounds: A >= 0
lb = [0]; 
ub = [Inf];

% Attempt lsqcurvefit
opts = optimoptions('lsqcurvefit','Display','off');
A_fit = lsqcurvefit(@(A,m) model(A,m), A0, mu, lambdaVals, lb, ub, opts);

%% ========== 4) Plot Results (with continuation extension) ==========

figure('Color','w','Position',[100 100 800 600]);
hold on; box on;

% Raw data
scatter(mu, lambdaVals, 40, 'filled', 'DisplayName', 'Data');

% Extend plotting range by 10%
extend_ratio = 0.1;
muRange = max(mu) - min(mu);
muFine = linspace(min(mu) - extend_ratio*muRange, ...
                  max(mu) + extend_ratio*muRange, 400);

% Fitted curve
plot(muFine, model(A_fit, muFine), 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');
yline(1, '--k', 'DisplayName', '\lambda = 1');

xlabel('\mu');
ylabel('\lambda(\mu)');
title(sprintf('Square-root Fit with Fixed K_c = %.3g', Kc), 'FontSize', 12);
legend('Location','best');

% Annotated box with fitted parameter
txt = sprintf('A = %.4g\nK_c = %.4g (fixed)', A_fit, Kc);
xl = xlim; yl = ylim;
text(xl(1) + 0.02*(xl(2)-xl(1)), yl(1) + 0.05*(yl(2)-yl(1)), txt, ...
     'FontSize', 11, 'BackgroundColor',[0.95 0.95 0.98], ...
     'EdgeColor',[0.8 0.8 0.9]);

saveas(gcf, 'fit_floquet_fixedKc.png');
fprintf('Plot complete. Saved as fit_floquet_fixedKc.png\n');
