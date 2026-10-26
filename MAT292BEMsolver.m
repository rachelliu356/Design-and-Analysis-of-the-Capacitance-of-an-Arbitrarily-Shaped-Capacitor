clear; close all;

% Permissivity
eps0 = 8.854187817e-12;

% Grid & domain
Nx = 80; Ny = 80; % number of grid points in each direction
% Note that in most of the results 
Lx = 100.0; Ly = 100.0; Lz = 1; % grid dimensions (meters)
                                 % Lz is approximate as the z component of
                                 % each plate is actually defined elsewhere
dx = Lx/Nx; dy = Ly/Ny; % spacing of mesh points
xl = linspace(0, Lx, Nx); % used for visualization of the geometries
yl = linspace(0, Ly, Ny); % xl and yl effectively form the x and y axes (and provide scaling)
zlb = zeros(Nx, Ny); % zlt holds the z component at mesh point (x, y) for the top surface
zlt = zeros(Nx, Ny); % zlb holds the z component at mesh point (x, y) for the bottom surface

% Here is where we define the plate geometry. I encourage you to play
% around with it to see the different visualizations! Tip: If you are
% trying to make a specific shape but it isn't working, this may be because the x and y values 
% don't have a wide enough "span" to be seen. You can try
% increasing the plate length, dividing the x or y values in your equation
% to make them smaller, or increasing top_valid_num and bot_valid_num.

% Here are a few geometries we have pre-defined for you! If you would like
% to visualize them, you can just write the function name in the z_bot or
% z_top equation as shown below. 
z_bot_flat = @(x, y) 0; % function that describes z of top and bottom plates
z_top_flat = @(x, y) Lz;

z_bot_wave1 = @(x, y) sin(x) + cos(y);
z_top_wave1 = @(x, y) sin(x) + cos(y) + Lz;

z_bot_wave2 = @(x, y) - x^1.5/100 + y^2/300;
z_top_wave2 = @(x, y) sin(x/5)*2 - cos(y/3) + Lz;

z_bot_wave3 = @(x, y) sin(x/5+2);
z_top_wave3 = @(x, y) cos(x/10) + Lz;

xy_circle = @(x, y) (x-Lx/2)^2 + (y-Ly/2)^2;
xy_triange = @(x, y) 15*x +10*y;

n_bumps = 10; % number of bumps for the following surfaces
z_bot_bump = @(x, y) sin(x/Lx*pi*n_bumps)/2;
z_top_bump = @(x, y) sin(x/Lx*pi*n_bumps)/10 + Lz;

% Here is where the geometry that is actually plotted is defined
z_bot = @(x, y) z_bot_bump (x, y); % function that describes z of top and bottom plates
z_top = @(x, y) z_top_bump (x, y); 
top_valid_num = 1000; % we define the "allowed xy region" for the surface through the equation xy_top <= top_valid_num
xy_top = @(x, y) xy_circle(x, y); % for example, xy_top = (x-Lx/2)^2 + (y-Ly/2)^2 <= r^2 = top_valid_num is the equation of a circle
bot_valid_num = 1000;
xy_bot= @(x, y) xy_circle(x, y);
delt_V = 10; % Voltage across the plates

% Initalizing the voltage at each mesh point
V = transpose(zeros(1, Nx*2*Ny)); % there are Nx*Ny mesh points on each surface, for a total of 2*Nx*Ny mesh points
V(Nx*Ny+1:Nx*Ny*2) = delt_V; % The first Nx*Ny mesh points are those on the bottom plate, while the second Nx*Ny mesh points (those in the range Nx*Ny + 1 to 2*Nx*Ny)
                            % are those on the top. These top mesh points
                            % will all be at 10 V, as this is the boundary
                            % condition we have chosen to enforce. Feel
                            % free to play around with the voltage and see
                            % how it infulences the results!

% Initalizing the r (position) matrix to contain the position vector at each mesh point
r = zeros(Nx*Ny*2, 3); % Here, the rows represent mesh point index, and the columns represent the x, y, and z components of the mesh point's position, respectively 
                        % Mesh point indexing is as follows: index = x*Nx + y for the bottom plate and index = x*Nx + y + Nx*Ny for the top
for i = 1:Nx
    for j = 1:Ny
        x_coord = (i-0.5)*dx; % x and y coordinates
        y_coord = (j-0.5)*dy;
        valid_top_coord = xy_top(x_coord, y_coord); % plug xy coordinate into the equation for the top and bottom surfaces to check if this point is "allowed"
        valid_bot_coord = xy_bot(x_coord, y_coord);

        % Top plate:
        if (valid_top_coord <= top_valid_num) % for a circle, x^2 + y^2 <= r^2 - this is what this checks
            zlt(i, j) = z_top(x_coord, y_coord); % set the z coordinate of each mesh point for visualization
            r((i-1)*Nx + j + Nx*Ny, :) = [x_coord, y_coord, z_top(i, j)]; % setting the coordinate of each mesh point inside the position matrix
        else  % if the mesh point falls outside of the allowed xy region, assign it NaN to have it be ignored by the rest of the program
            zlt(i, j) = NaN;
            r((i-1)*Nx + j + Nx*Ny, :) = [NaN, NaN, NaN];
        end

        % Bottom plate:
        if (valid_bot_coord <= bot_valid_num) % Checking if the coordinates are valid for the bottom plate
            zlb(i, j) = z_bot(x_coord, y_coord); % setting the z coordinate of each mesh point on the bottom plate for visualization
            r((i-1)*Ny + j, :) = [x_coord, y_coord, z_bot(i, j)]; % setting position matrix at the bottom plate mesh point
        else  % if the mesh point falls outside of the allowed xy region, set to NaN
            zlb(i, j) = NaN;
            r((i-1)*Ny + j, :) = [NaN, NaN, NaN];
        end
    end
end

% Identify valid (non-NaN) patches
valid = ~any(isnan(r), 2)

A_patch = zeros(size(valid));

for i = 1:length(valid)
    if ~valid(i)
        continue
    end

    neighbors = [];

    % left/right (±Ny)
    if i-Ny >= 1 && valid(i-Ny)
        neighbors(end+1) = norm(r(i,:) - r(i-Ny,:));
    end
    if i+Ny <= length(valid) && valid(i+Ny)
        neighbors(end+1) = norm(r(i,:) - r(i+Ny,:));
    end

    % forward/back (±1)
    if i-1 >= 1 && valid(i-1)
        neighbors(end+1) = norm(r(i,:) - r(i-1,:));
    end
    if i+1 <= length(valid) && valid(i+1)
        neighbors(end+1) = norm(r(i,:) - r(i+1,:));
    end

    % Approximate area from mean neighbor spacing
    if numel(neighbors) >= 2
        A_patch(i) = mean(neighbors)^2;
    else
        A_patch(i) = dx * dy;  % fallback (edge/corner)
    end
end


% G is the influence matrix used to relate V, the voltage, to rho, charge
% density, at each point
% Gij = Aj/(4 pi epsilon |ri - rj|), ri is vector (from position matrix r)
% to center of patch i, rj is vector to center of patch j

% there are Nx*Ny patches on each plate, so 2*Nx*Ny patches in total
% each patch influences every other patch, so G will be a 2*Nx*Ny by
% 2*Nx*Ny matrix
% For geometries where some areas have been set to NaN, this causes issues, so we
% build G only for valid patches - non-NaN

r_valid = r(valid, :);
A_valid = A_patch(valid);
V_valid = V(valid);

N = sum(valid);

G = zeros(N, N);

for i = 1:N
    for j = 1:N
        dR = r_valid(i,:) - r_valid(j,:);
        if i ~= j
            G(i,j) = A_valid(j) / (4*pi*eps0*norm(dR));
        else
            % Self-term
            G(i,j) = sqrt(A_valid(j)) / (4*pi*eps0);
        end
    end
end

% Removing the nullspace by grounding one patch
G(1,:) = 0;
G(1,1) = 1;
V_valid(1) = 0;

% Solve for charge density
rho_valid = G \ V_valid
rho = zeros(size(V))
rho(valid) = rho_valid


rho_top = abs(rho(Ny*Nx+1:2*Ny*Nx)) % Charge densities on top plate
rho_bot = abs(rho(1:Ny*Nx)) % Charge densities on bottom plate
A_patch_top = A_patch(Nx*Ny+1:2*Nx*Ny)
A_patch_bot = A_patch(1:Nx*Ny)

% Total charge on a plate is the sum of charge density on each patch times area of patch
Q_tot_top = dot(rho_top, A_patch_top) % equal patch sizes, so we can "bring out" the A_patch multiplication
C_top = Q_tot_top/delt_V % Capacitance is Q/V
Q_tot_bot = dot(rho_bot, A_patch_bot) % inaccurate We take the absolute value because one plate will have positive charge and the other will have negative,
C_bot = Q_tot_bot/delt_V                % taking the absolute value makes the charge density comparision simpler

C = (C_bot + C_top)/2 % averaging the capacitance calculated from the top and bottom plates
Percent_Error_Estimate = abs(Q_tot_top-Q_tot_bot)/(Q_tot_top+Q_tot_bot)*200 % Theoretically, the total charge on the top and bottom should be the same
                                                                            % We are using their difference when calculated numerically as an apporimate error heuristic
C_ana = Lx*Ly/Lz*eps0 % Capacitance of a parallel plate capacitor assuming plates are long enough for edge effects to be negligible. We calulate this as a reference

% Initalize the charge density colormap
Colour_top = zeros(Nx, Ny);
Colour_bot = zeros(Nx, Ny);

for i = 1:(Nx)
    for j = 1:(Nx)
        % Colour at point x, y will correspond to the charge density at the same point
        Colour_top(i, j) = abs(rho_top((i-1)*Nx + j));
        Colour_bot(i, j) = abs(rho_bot((i-1)*Nx + j)); 
    end
end

% Plot the surfaces with the charge density
figure;
surf(xl, yl, zlb, Colour_bot, 'EdgeColor','none');
hold on;
surf(xl, yl, zlt, Colour_top, 'EdgeColor','none');
colorbar;
xlabel('x');
ylabel('y');
zlabel('z');
title('Charge Distribution');
hold off;