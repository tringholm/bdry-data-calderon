clear 
close all

%%%%%%%%%%%%% Specify parameters %%%%%%%%%%%%%

% Thresholds under which H, N, A are considered zero
H2_thresh = 1E-5; 
N_thresh = 1E-5;
A_thresh = 1E-5;
n_diff_thresh = 1E-1;

% Upper/lower bounds on sigma
sigma_max = 20;
sigma_min = 0.1;

%Constant specifying the fineness of the FEM mesh
meshres = 0.01;

%Constant specifying number of electrodes
M = 100;
h = 2*pi/M;
theta = 0:h:2*pi-h;


%%%%%%%%%%%%% Solve forward problem, get measurements %%%%%%%%%%%%%

% Generate FEM model and solve it. Solution is stored in "results"
[model,results] = fwdProblem(meshres);

% Extract boundary data from model and results
[N,A,H1,H2] = getBoundaryDataCrimeless(model,results,M);


%%%%%%%%%%%%% Reconstruct conductivity data %%%%%%%%%%%%%


% Indices for the cases in Theorem 3
zeroH_ind = (abs(H2) < H2_thresh);                      % Case a)
zeroN_ind = logical((abs(N) < N_thresh).*(~zeroH_ind)); % Case b)
zeroA_ind = logical((abs(A) < A_thresh).*(~zeroH_ind)); % Case c)
double_ind = logical(ones(size(zeroH_ind)) - zeroH_ind - zeroN_ind - zeroA_ind); % case d)
fixed_ind = zeros(size(zeroH_ind)); % Indices of unambiguous points

% In the following, q = 2, p = 2
sig_q2 = zeros(size(N)); % For the final estimate
sig_q2(zeroH_ind) = -1; % Not getting any information in case a)
sig_q2(zeroN_ind) = H2(zeroN_ind)./A(zeroN_ind).^2; % Case b)
fixed_ind(zeroN_ind) = 1; % Mark as unambiguous
sig_q2(zeroA_ind) = N(zeroA_ind).^2./H2(zeroA_ind); % Case c)
fixed_ind(zeroA_ind) = 1; % Mark as unambiguous

plotData(theta,sig_q2,'Before double points','\theta','\sigma','beforeDouble')

%%%%% Tackle the ambiguous case %%%%

% Upper/lower estimates for n and sigma
n_plus = 0.5*(H2./N + sqrt(H2.^2./N.^2 - 4*A.^2));
n_minus = 0.5*(H2./N - sqrt(H2.^2./N.^2 - 4*A.^2));
sig_plus = zeros(size(N));
sig_minus = zeros(size(N));
sig_plus(double_ind) = real(2*N(double_ind).^2./(H2(double_ind) + sqrt(H2(double_ind).^2 - 4*(A(double_ind).^2).*(N(double_ind).^2))));
sig_minus(double_ind) = real(2*N(double_ind).^2./(H2(double_ind) - sqrt(H2(double_ind).^2 - 4*(A(double_ind).^2).*(N(double_ind).^2))));

% Indices where the two estimates agree are alright
agrees = logical((abs(sig_plus - sig_minus) < 1E-12).*double_ind);
sig_q2(agrees) = sig_plus(agrees);
fixed_ind(agrees) = 1;

% Check for over/undershoots in the estimates, fix unambiguous points
oob_plus_ind = (sig_plus > sigma_max)+(sig_plus < sigma_min);
oob_minus_ind = (sig_minus > sigma_max)+(sig_minus < sigma_min);
use_minus_ind = logical(oob_plus_ind.*(1 - oob_minus_ind).*double_ind); % If sig_plus is out of bounds but sig_minus is not
use_plus_ind = logical(oob_minus_ind.*(1 - oob_plus_ind).*double_ind);  % If sig_minus is out of bounds but sig_plus is not
sig_q2(use_minus_ind) = sig_minus(use_minus_ind);
sig_q2(use_plus_ind) = sig_plus(use_plus_ind);
fixed_ind(use_minus_ind) = 1;
fixed_ind(use_plus_ind) = 1;

% Find points signifying n^2 = A^2 nearby. These cannot be labeled with
% either plus or minus
n_diff = abs(n_plus) - abs(n_minus);
n_diff = [n_diff(end); n_diff; n_diff(1)];
n_diff_min_ind = zeros(size(n_plus));
for i = 2:(length(n_diff)-1)
    if n_diff(i) < n_diff(i+1) && n_diff(i) < n_diff(i-1) && abs(n_diff(i)) < n_diff_thresh
       n_diff_min_ind(i-1) = 1;
    end
end

plotData(theta,sig_q2,'Before population','\theta','\sigma','beforePopulation')

% Populate the remaining undecided points
num_pts = length(N);
for i = 1:num_pts
    if ~fixed_ind(i) && ~zeroH_ind(i) && ~n_diff_min_ind(i)  % Point is not fixed yet and should be fixed
        
        %%%% Find right stopping point
        cyclic_right = 1;
        for k = (i+1):num_pts
            if zeroH_ind(k) || ~double_ind(k) || n_diff_min_ind(k) || fixed_ind(k)% Search stops if meeting undecideable point, non-double point or unlabelable
                k_right = k;
                cyclic_right = 0;
                break
            end
        end
        
        % Need to take care of cyclic possibility
        if cyclic_right
            for k = 1:num_pts
                if zeroH_ind(k) || ~double_ind(k) || n_diff_min_ind(k) || fixed_ind(k) % Search stops if meeting undecideable point, non-double point or unlabelable
                    cyc_k_right = 0;
                    break
                end
            end
        else
            cyc_k_right = 0;
        end
        
       %%%% Find left stopping point
        
        cyclic_left = 1;
        for k = (i-1):-1:1
            if zeroH_ind(k) || ~double_ind(k) || n_diff_min_ind(k) || fixed_ind(k)  % Search stops if meeting undecideable point, non-double point or unlabelable
                k_left = k;
                cyclic_left = 0;
                break
            end
        end
        
        % Need to take care of cyclic possibility
        if cyclic_left
            for k = num_pts:-1:1
                if zeroH_ind(k) || ~double_ind(k) || n_diff_min_ind(k) || fixed_ind(k) % Search stops if meeting undecideable point, non-double point or unlabelable
                    cyc_k_left = 0;
                    break
                 end
            end
        else
            cyc_k_left = 0;
        end
        
        % Check if right and left stopping points don't disagree, act
        % accordingly
        
        if cyc_k_left % Cyclic case, leftward search
            k_left_check = cyc_k_left;
        else
            k_left_check = k_left;
        end
        
        if cyc_k_right % Cyclic case, leftward search
            k_right_check = cyc_k_right;
        else
            k_right_check = k_right;
        end
        if (use_plus_ind(k_right_check) && ~use_minus_ind(k_left_check)) || (~use_minus_ind(k_right_check) && use_plus_ind(k_left_check))
            sig_q2(k_left+1:k_right-1) = sig_plus(k_left+1:k_right-1);
            use_plus_ind(k_left+1:k_right-1) = 1;
            fixed_ind(k_left+1:k_right-1) = 1;
            if cyc_k_left % Cyclic case, leftward search
                sig_q2(cyc_k_left+1:num_pts) = sig_plus(cyc_k_left+1:num_pts);
                use_plus_ind(cyc_k_left+1:num_pts) = 1;
                fixed_ind(cyc_k_left+1:num_pts) = 1;
            end
            if cyc_k_right % Cyclic case, rightward search
                sig_q2(1:cyc_k_right-1) = sig_plus(1:cyc_k_right-1);
                use_plus_ind(1:cyc_k_right-1) = 1;
                fixed_ind(1:cyc_k_right-1) = 1;
            end
        elseif (use_minus_ind(k_right_check) && ~use_plus_ind(k_left_check)) || (~use_plus_ind(k_right_check) && use_minus_ind(k_left_check))
            sig_q2(k_left+1:k_right-1) = sig_minus(k_left+1:k_right-1);
            use_minus_ind(k_left+1:k_right-1) = 1;
            fixed_ind(k_left+1:k_right-1) = 1;
            if cyc_k_left % Cyclic addition, leftward search
                sig_q2(cyc_k_left+1:num_pts) = sig_minus(cyc_k_left+1:num_pts);
                use_minus_ind(cyc_k_left+1:num_pts) = 1;
                fixed_ind(cyc_k_left+1:num_pts) = 1;
            end
            if cyc_k_right % Cyclic addition, rightward search
                sig_q2(1:cyc_k_right-1) = sig_minus(1:cyc_k_right-1);
                use_minus_ind(1:cyc_k_right-1) = 1;
                fixed_ind(1:cyc_k_right-1) = 1;
            end
        end
        
    end
end

plotData(theta,sig_q2,'Before interpolation','\theta','\sigma','beforeInterpolation')


% Interpolate the remaining points using the other established points
for i = 1:num_pts
   if ~fixed_ind(i)
       for j = 1:(num_pts - i)
            k = i + j;
            if fixed_ind(k)
               k_right = k; 
               break
            end
        end
        for j = -1:-1:(1 - i)
            k = i + j;
            if fixed_ind(k)
                k_left = k;
                break
            end
        end
        gap_size = k_right - k_left;
        l_wt = gap_size -1;
        r_wt = 1;
        l_val = sig_q2(k_left);
        r_val = sig_q2(k_right);
        for j = 1:gap_size-1
            sig_q2(i + j -1) = l_val*l_wt/gap_size + r_val*r_wt/gap_size;
            l_wt = l_wt - 1;
            r_wt = r_wt + 1;
            fixed_ind(i + j - 1) = 1;
        end
   end
end

plotData(theta,sig_q2,'After interpolation','\theta','\sigma','afterInterpolation')
plotData(theta,smoothGauss(sig_q2),'After smoothing','\theta','\sigma','afterSmoothing')
plotDataComp(theta,smoothGauss(sig_q2),'Comparison','\theta','\sigma','afterSmoothingComparison')

% figure 
% plot(sig_plus);
% hold on;
% plot(sig_minus)
% axis([0, length(sig_plus), 0, 5])
% plot(n_plus - n_minus) 
% 

