%% Domain parameters
clear; clc;

% Defines resolution of the domain (num of grid points)
gridnum=251; % evaluated at 101, 151, 201, and 251
Nx = gridnum; Ny = gridnum;
x = linspace(-0.05, +0.05, Nx);   % creates x- and y- coordinates evenly spanning domain
y = linspace(-0.05, +0.05, Ny);
dx = x(2)-x(1);                   % determines dx and dy from first two points
dy = y(2)-y(1);

% defines mu properties
mu_0 = 4*pi*1e-7;
mu_r_iron = 5000;

%create mu matrix of mu_0's
mu = mu_0*ones(Ny,Nx);

%create current density matrix of 0's
Jz = zeros(Ny,Nx);

% converts x and y vectors into 2D matrix [nb: X(j,i)= x at column i]
[X,Y] = meshgrid(x,y);

%% Place iron triangles
iron1 = (X >= -0.05 & X <= 0 & Y >= 0 & Y <= 0.05 & Y >= X + 0.05);
iron2 = flipud(fliplr(iron1));

% mu = mu_0*mu_r_iron wherever iron1 or iron2 exist
mu(iron1 | iron2) = mu_0*mu_r_iron;

%% Place conductors
cond1 = (X>=0.02  & X<0.03  & Y>=0.02  & Y<0.03);
cond2 = flipud(fliplr(cond1));

% Normalize Jz so the discrete conductors carry exactly +/-10 A
Acond1_discrete = nnz(cond1)*dx*dy;
Acond2_discrete = nnz(cond2)*dx*dy;

Jz(cond1) = +10/Acond1_discrete;
Jz(cond2) = -10/Acond2_discrete;

%% Check discrete current representation
fprintf('Conductor 1 nodes: %d\n', nnz(cond1));
fprintf('Discrete current 1: %.6f A\n', sum(Jz(cond1))*dx*dy);
fprintf('Discrete current 2: %.6f A\n', sum(Jz(cond2))*dx*dy);

%% Solve for Az using Gauss-Seidel
Az = zeros(Ny,Nx); %set all values to zero initially
maxiterations = 50000;
tol = 1e-12; % convergence tolerance on update

for iter = 1:maxiterations
    Az_old = Az; % stores the previous iteration's Az so the convergence check can compare old vs updated values
    for i = 2:Nx-1  % for both i and j: update only interior points (leave one point around entire domain on boundary equal to 0)
        for j = 2:Ny-1
            
            % averaging of mu in x,y directions
            muR = (mu(j,i) + mu(j,i+1)) / 2;
            muL = (mu(j,i) + mu(j,i-1)) / 2;
            muU = (mu(j,i) + mu(j+1,i)) / 2;
            muD = (mu(j,i) + mu(j-1,i)) / 2;
            M = 1/muR + 1/muL + 1/muU + 1/muD;

            Az(j,i) = ((Az(j,i+1)/muR) + (Az(j,i-1)/muL) + (Az(j+1,i)/muU) + (Az(j-1,i)/muD) + (dx^2)*Jz(j,i)) / M;
        end
    end

    change = abs(Az - Az_old);
    max_change = max(change(:));
    if max_change < tol %check the max absolute change between iterations, 
        break;  % if within tolerance, consider acceptable and break
    end
end

%% Check solver convergence
max_update = max(abs(Az(:) - Az_old(:)));

fprintf('Gauss-Seidel iterations: %d\n', iter);
fprintf('Final maximum update: %.3e\n', max_update);

if iter == maxiterations && max_update >= tol
    warning('Solver reached maximum iterations without meeting tolerance.');
end

%% Find Az gradient
dAz_dx = zeros(Ny, Nx);
dAz_dy = zeros(Ny, Nx);

for i = 2:Nx-1 % range is for same reason as above in Gauss-Seidel
    for j = 2:Ny-1       
        dAz_dx(j,i) = (Az(j, i+1) - Az(j, i-1)) / (2*dx);   % ∂Az/∂x
        dAz_dy(j,i) = (Az(j+1, i) - Az(j-1, i)) / (2*dy);   % ∂Az/∂y
    end
end

%% Compute B field vectors and magnitude
Bx = dAz_dy;
By = -dAz_dx;

Bmag = sqrt(Bx.^2 + By.^2);

%% Subsample quiver plots for clarity
step = 10;
Xq = X(1:step:end,1:step:end);
Yq = Y(1:step:end,1:step:end);
Bxq = Bx(1:step:end,1:step:end);
Byq = By(1:step:end,1:step:end);
dAz_dxq = dAz_dx(1:step:end,1:step:end);
dAz_dyq = dAz_dy(1:step:end,1:step:end);

%% Plot numerical solution
figure('Position',[200 100 1100 850]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% B vector field
nexttile
quiver(Xq,Yq,Bxq,Byq,'AutoScale','on','AutoScaleFactor',1)
axis equal
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('x (m)')
ylabel('y (m)')
title('B-Field Vectors')

% B magnitude
nexttile
contourf(X,Y,Bmag,50,'LineColor','none')
axis equal
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('x (m)')
ylabel('y (m)')
colorbar
title('B-Field Magnitude')

% Az gradient
nexttile
quiver(Xq,Yq,dAz_dxq,dAz_dyq,'AutoScale','on','AutoScaleFactor',1)
axis equal
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('x (m)')
ylabel('y (m)')
title('Gradient of A_z')

% Az
nexttile
contourf(X,Y,Az,50,'LineColor','none')
axis equal
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('x (m)')
ylabel('y (m)')
colorbar
title('A_z')

exportgraphics(gcf,'figures/numerical_solution.png','Resolution',300)

%% Define circular conductors 
radius = 0.005;
x1 = 0.025;  y1 = 0.025; % center of +10 A conductor
x2 = -0.025; y2 = -0.025; % center of -10 A conductor
I1 = +10;
I2 = -10;

%% Define distances from conductor centers
radial_x1 = X - x1; % x-component of radial vector from conductor 1
radial_y1 = Y - y1; % y-component of radial vector from conductor 1
rho1 = sqrt(radial_x1.^2 + radial_y1.^2);  % rho = |radial vector|

radial_x2 = X - x2; % repeated for other conductor
radial_y2 = Y - y2; 
rho2 = sqrt(radial_x2.^2 + radial_y2.^2);

% define each conductor by comparing rho to the radius
cond1_in  = rho1 <= radius & rho1 > 0;
cond1_out = rho1 > radius;
cond2_in = rho2 <= radius & rho2 > 0;
cond2_out = rho2 > radius;

%% Bphi magnitudes
% create matrices for the B-fields of each conductor
Bphi1 = zeros(size(X));
Bphi2 = zeros(size(X));

Bphi1(cond1_in)  = (mu_0*I1 .* rho1(cond1_in)) ./ (2*pi*radius^2);
Bphi1(cond1_out) = (mu_0*I1) ./ (2*pi*rho1(cond1_out));

Bphi2(cond2_in)  = (mu_0*I2 .* rho2(cond2_in)) ./ (2*pi*radius^2);
Bphi2(cond2_out) = (mu_0*I2) ./ (2*pi*rho2(cond2_out));

%% Convert to Cartesian components
% Copy rho arrays
rho1_copy = rho1;
rho2_copy = rho2;

% Anywhere the radius is 0, replace with 1 to avoid division by 0
rho1_copy(rho1_copy == 0) = 1;
rho2_copy(rho2_copy == 0) = 1;

% Create unit vectors in phi direction
phi_dir_x1 = -radial_y1 ./ rho1_copy;
phi_dir_y1 =  radial_x1 ./ rho1_copy;
phi_dir_x2 = -radial_y2 ./ rho2_copy;
phi_dir_y2 =  radial_x2 ./ rho2_copy;

% Create B vectors
Bx1 = Bphi1 .* phi_dir_x1;
By1 = Bphi1 .* phi_dir_y1;
Bx2 = Bphi2 .* phi_dir_x2;
By2 = Bphi2 .* phi_dir_y2;

% Force exact-zero field at the conductor centers
Bx1(rho1==0)=0; By1(rho1==0)=0;
Bx2(rho2==0)=0; By2(rho2==0)=0;

% Sum the conductors and find |B|
Bx_total = Bx1 + Bx2;
By_total = By1 + By2;
Bmag = sqrt(Bx_total.^2 + By_total.^2);

%% Plot results for analytical solution
Xq = X(1:step:end,1:step:end);
Yq = Y(1:step:end,1:step:end);
Bxq = Bx_total(1:step:end,1:step:end);
Byq = By_total(1:step:end,1:step:end);

figure('Position',[200 200 1600 700]);

subplot(1,2,1)
quiver(Xq,Yq,Bxq,Byq,'AutoScale','on','AutoScaleFactor',1)
axis equal
xlabel('x (m)'); ylabel('y (m)');
title('B Quiver plot (Analytical Solution)')
ax1 = gca;

subplot(1,2,2)
contourf(X,Y,Bmag,50,'LineColor','none')
axis equal
xlabel('x (m)'); ylabel('y (m)');
colorbar
title('|B| (Analytical Solution)')
ax2 = gca;


%% Compute absolute and relative error between numerical and analytical B

% Absolute error components
dBx = Bx - Bx_total;
dBy = By - By_total;

% Absolute error magnitude
d_Bmag = sqrt(dBx.^2 + dBy.^2);

% Relative error to avoid division by zero
rel_err = zeros(size(Bmag));
nonzero_B = Bmag > 1e-12;
rel_err(nonzero_B) = d_Bmag(nonzero_B) ./ Bmag(nonzero_B);


%% Error summary values

max_abs_error = max(d_Bmag(:));

valid = rel_err(nonzero_B);
max_rel_error = max(valid);

min_abs_error = min(d_Bmag(:));
[min_rel_error, ~] = min(valid);

fprintf('Maximum absolute error between numerical and analytical solutions: %.2g\n', max_abs_error);
fprintf('Minimum absolute error between numerical and analytical solutions: %.2g\n', min_abs_error);
fprintf('\n');
fprintf('Maximum relative error between numerical and analytical solutions: %.2g\n', max_rel_error);
fprintf('Minimum relative error between numerical and analytical solutions: %.2g\n', min_rel_error);

%% Plot absolute error
figure('Position',[300 300 1600 700]);

subplot(1,2,1)
contourf(X, Y, d_Bmag, 50, 'LineColor','none')
axis equal
xlabel('x (m)')
ylabel('y (m)')
colorbar
title('|B_{num} - B_{analytic}|')

subplot(1,2,2)
contourf(X, Y, rel_err, 50, 'LineColor','none')
axis equal
xlabel('x (m)')
ylabel('y (m)')
colorbar
title('Relative Error |ΔB| / |B_{analytic}|')

%% Validate Ampere's law
% find H
Hx = Bx ./ mu;
Hy = By ./ mu;

% Each row: [x-center, y-center, expected current]
conductors = [ 0.025   0.025   10;
              -0.025  -0.025  -10];

% Radii path must be larger than conductor half-diagonal (~7.07 mm)
% so the path encloses the full conductor current
radii = [0.008 0.010 0.012 0.015 0.018 0.020]; % in meters (8-20 mm)

% Number of points used around each circular path
Ntheta_values = [250 500 1000 2000 4000 8000];

results = [];

fprintf('\n===== AMPERE VALIDATION SWEEP =====\n');

for c = 1:size(conductors,1) % select conductors and expected current

    xc = conductors(c,1);
    yc = conductors(c,2);
    I_expected = conductors(c,3);

    fprintf('\nConductor %d (expected %.1f A)\n', c, I_expected);

    for n = 1:length(Ntheta_values) % selects number of points on circular path

        Ntheta = Ntheta_values(n);
        errors = zeros(size(radii));

        for k = 1:length(radii) % selects radius of circular path

            r = radii(k);

            % Constructs circular path around the conductor
            theta = linspace(0, 2*pi, Ntheta+1);
            xq = xc + r*cos(theta);
            yq = yc + r*sin(theta);

            % Interpolate H from the rectangular grid onto points along the circular path
            % H_x​(x_i​,y_i​) → H_x​(x(θ),y(θ))
            Hxq = interp2(x, y, Hx, xq, yq, 'linear');
            Hyq = interp2(x, y, Hy, xq, yq, 'linear');

            % Components of dl/dtheta along the circular path
            dlx = -r*sin(theta);
            dly =  r*cos(theta);

            % Ampere's law
            integrand = Hxq.*dlx + Hyq.*dly; % form H dot (dl/dtheta)
            % Numerically evaluate integral H dot dl around the closed path
            I_enc = trapz(theta, integrand);

            % Percentage error
            error_pct = 100*abs(I_enc - I_expected)/abs(I_expected);
            errors(k) = error_pct;

            results = [results; c Ntheta r I_enc error_pct];
        end

        fprintf('Ntheta = %5d: mean error = %.5f %%, max error = %.5f %%\n', ...
            Ntheta, mean(errors), max(errors));
    end
end

%% Overall error results
fprintf('\n===== OVERALL =====\n');
fprintf('Mean error over all tests: %.5f %%\n', mean(results(:,5)));
fprintf('Maximum error over all tests: %.5f %%\n', max(results(:,5)));
fprintf('Minimum error over all tests: %.5f %%\n', min(results(:,5)));

%% Plot Ampere's law validation
Ntheta_final = max(Ntheta_values);

figure('Position',[200 200 700 500])
hold on

for c = 1:size(conductors,1)
    rows = results(results(:,1) == c & ...
                   results(:,2) == Ntheta_final, :);

    plot(rows(:,3)*1000, rows(:,5), '-o')
end

xlabel('Contour Radius (mm)')
ylabel('Ampere''s Law Error (%)')
title('Ampere''s Law Validation')
legend('Conductor 1', 'Conductor 2', 'Location', 'best')
grid on

exportgraphics(gcf,'figures/ampere_validation.png','Resolution',300)
