clear; close all; clc

Naca = input("NACA Airfoil Number: ",'s'); % Prompting User for NACA number
alpha = input("Angle of Attack (degrees): "); % Prompting User for angle of attack
alpha = deg2rad(-alpha);
U_inf  = input("Flow Velocity (m/s): "); % Setting the flow velocity

[X, Y, ~, ~] = NACABuild(Naca, 300); % Building the NACA
X = X - 0.5; % Moving the airfoil so its centered at x=0
z = X + 1i*Y; % Building the function
z = flip(z);
z(end) = []; % Getting rid of the ovelaping points

p = polygon(z); % Setting up the polygon
f = extermap(p); % Mapping it to a disk


figure() 
plot(f) % Plotting sc transformed plane
ylabel('Imaginary w plane')
xlabel('Real w place')
title('Schwarz–Christoffel Mapped Plane')
print("Figures/Mapped_Plane_of_" + convertCharsToStrings(Naca),'-dpng')

%% Build polar grid INSIDE unit disk (exterior airfoil = interior disk)
Nr     = 200; % Setting the number of radial grid points
Ntheta = 400; % Setting the angular grid points
r      = linspace(0.001, 0.95, Nr); % setting r (avoid r=0 (maps to Inf))
theta  = linspace(0, 2*pi, Ntheta); % setting theta 
[R, Theta] = meshgrid(r, theta); % Creating a meshgrid with all the points
W = R .* exp(1i*Theta); % Converting polar to complex cartesian

%% Forward map w -> z
Z_flat = feval(f, W(:)); % Flattens Z back into a vector
Z = reshape(Z_flat, size(W)); % Reshapes the vector back into the 2d grid
Xp = real(Z); % Seperating the real component
Yp = imag(Z); % Seperating the imaginary component

%% Complex potential
C_scale = abs(-(1e-5)^2 * evaldiff(f, 1e-5)); % Finding the infinity scaling factor
U_comp = U_inf * C_scale; % Calculating the new velocity due to the scaling factor
[~, idx_TE] = max(real(vertex(p))); % Finding the index of the trailing edge

w_pre = prevertex(f); 
beta_TE = angle(w_pre(idx_TE)); % Calclating the angle the trailing edge makes


Gamma = 4 * pi * U_comp * sin(alpha - beta_TE); % Calculating the circulation

Fw  = U_comp*(exp(-1i*alpha)./W + exp(1i*alpha)*W) ...
      - 1i*Gamma/(2*pi)*log(W); % Finds the complex potential at every point
% exp(-1i*alpha)./W is the uniform freestream flow. exp(1i*alpha).*W is the 
% reflection term that enforces the no-penetration condition on the unit 
% circle. -1i*Gamma/2*pi*log(W) adds the vortex that generates lift.
Psi = imag(Fw); % Extracts the stream function from the complex potential

%% Clip extreme values near airfoil surface
Psi(abs(W) > 0.98) = NaN; % Sets values close to airfoil surface to NaN 
% so the data isn't noisy. (Numerical Instability)

%% Plot
figure; hold on;

psi_levels = linspace(-2, 2, 50); % Setting the stream lines
contour(Xp, -Yp, Psi, psi_levels, 'r', 'LineWidth', 0.7); % Plotting the streamlines
fill(real(z), -imag(z), [0.3 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1.5); % Potting the airfoil

axis equal % Setting axis to be equal so airfoil doesn't stretch
xlim([-1 1]); ylim([-1.2 1.2]) % Limiting axis
xlabel('x/c'); ylabel('y/c') % Labeling axis
title(['NACA 0012 Streamlines (\alpha = ', num2str(rad2deg(-alpha)), '°)']) % Labeling plot
grid on; box on
print("Figures/Streamlines_Over_" + convertCharsToStrings(Naca),'-dpng')

%% Calcualting parameters
% Finding the change in our function
dF_dw = U_comp * (-exp(-1i*alpha)./W.^2 + exp(1i*alpha)) - 1i*Gamma./(2*pi*W);

% Flattening the change in the function to be a vector
dz_dw_flat = evaldiff(f, W(:));

% Finding the change in z due to the distortion
dz_dw = reshape(dz_dw_flat, size(W));

% Finding the complex velocity using chain rule
V_complex = dF_dw ./ dz_dw;

% Calculating the magnitude of the complex velocity
Vmag = abs(V_complex);

% Calculating the grid form of the coefficient of pressure
Cp_grid = 1 - (Vmag / U_inf).^2;

% Setting these variables to surfaces
X_surf  = Xp(:, end);
Cp_surf = Cp_grid(:, end);


% Finding the index of the leading edge
[~, LE_idx] = min(X_surf);

% Plotting the results
figure('Name', 'Cp Distribution'); hold on; grid on;

% Plot Upper surface (from TE forward to LE)
plot(X_surf(1:LE_idx), Cp_surf(1:LE_idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Upper Surface');

% Plot Lower surface (from LE backward to TE)
plot(X_surf(LE_idx:end), Cp_surf(LE_idx:end), 'r--', 'LineWidth', 2, 'DisplayName', 'Lower Surface');


% Plotting the results
xlabel('x/c');
ylabel('C_p');
title(['Surface Pressure Coefficient | \alpha = ', num2str(rad2deg(-alpha)), '°']);
legend('Location', 'northeast');
set(gca, 'YDir', 'reverse');
print("Figures/C_p_Over_" + convertCharsToStrings(Naca),'-dpng')

% Calculating chord
chord = max(X_surf) - min(X_surf);

% Finding the Coefficient of lift
Cl_integration = trapz(X_surf, Cp_surf) / chord;

