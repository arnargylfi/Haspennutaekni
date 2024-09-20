clear; close all;
% Problem 3 - Floating conductor simulation.
% Power line consists of a 4 wire and represents one phase.
% Floating conductor as a metal sphere with diameter = 65cm (radius=0.325m).

% Create a PDE model for electrostatics
emagmodel = createpde("electromagnetic","electrostatic");

% Manually define the geometry using "outerRect" command. 3 specifies 2D
% shape, and 4 specifies the number of corner points (vertices), so a
% rectangle.
% Outer rectangle (3x6) and inner rectangle (0.3x1.5)
% Outer rectangle (larger frame): first four are x-coordinates in counterclockwise order
% The next four 1 are the corresponding y-coordinates.
%outerRect = [3; 4; -1.5; 1.5; 1.5; -1.5; -3 ; -3; 3; 3];
outerRect = [3; 4; -15; 15; 15; -15; 0 ; 0; 18; 18];

% Power line
% First four are x-coordinates in counterclockwise order
% The next four are the corresponding y-coordinates. 
%innerRect = [3; 4; -0.3; 0.3; 0.3; -0.3; -1.5; -1.5; 1.5; 1.5]; 
%innerRect = [3; 4; -0.0003; 0.0003; 0.0003; -0.0003; -0.0015; -0.0015; 0.0015; 0.0015]; 

% PowerLine conductor 1 (one phase) is defined by a circle by 1, x_center, y_center, radius
%Conductor = [1; 0; 6.8; 0.225; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% PowerLine conductor 1 (one phase) is defined by a circle by 1, x_center, y_center, radius
Conductor1 = [1; -0.225; 7.025; 0.01325; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% PowerLine conductor 2 (one phase) is defined by a circle by 1, x_center, y_center, radius
Conductor2 = [1; 0.225; 7.025; 0.01325; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% PowerLine conductor 3 (one phase) is defined by a circle by 1, x_center, y_center, radius
Conductor3 = [1; 0.225; 6.575; 0.01325; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% PowerLine conductor 4 (one phase) is defined by a circle by 1, x_center, y_center, radius
Conductor4 = [1; -0.225; 6.575; 0.01325; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

sphereRadius=0.325;
sphereX=2;  % X coordinate for the floating conductor sphere
sphereY=6.8;
% Floating conductor (metal sphere) is defined by a circle by 1, x_center, y_center, radius
FloatCond = [1; sphereX; sphereY; sphereRadius; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% Combine all items into the geometry matrix
%gd = [outerRect Conductor FloatCond]; % With a single large conductor
gd = [outerRect Conductor1 Conductor2 Conductor3 Conductor4 FloatCond];
%gd = [outerRect Conductor1 Conductor2 Conductor3 Conductor4];   % Without sphere

% Set up the name-space for the shapes
%ns = char('outer', 'cond', 'floatcond'); % With a single large conductor
ns = char('outer', 'cond1', 'cond2', 'cond3', 'cond4', 'floatcond');
%ns = char('outer', 'cond1', 'cond2', 'cond3', 'cond4');         % Without sphere
ns = ns';

% Specify the set formula (how the shapes are combined)
% 'outer - inner' means the outer rectangle with the inner rectangle removed (gap)
%sf = 'outer - cond - floatcond'; % With a single large conductor
sf = 'outer - cond1 -cond2 - cond3 - cond4 + floatcond';
%sf = 'outer - cond1 -cond2 - cond3 - cond4';                    % Without sphere


% Create the decomposed geometry using decsg
g = decsg(gd, sf, ns);

% Assign the geometry to the PDE model
geometryFromEdges(emagmodel, g);

% Visualize the geometry with edge labels
figure;
pdegplot(emagmodel, "EdgeLabels", "on")
hold on;
pdegplot(emagmodel, "FaceLabels", "on")
axis equal
xlabel('[m]') 
ylabel('[m]') 

% Specify the vacuum permittivity value in the SI system of units
emagmodel.VacuumPermittivity = 8.8541878128E-12;

% Specify the relative permittivity of the material
%electromagneticProperties(emagmodel, "RelativePermittivity", 1.00059);

% Inner region (air gap with relative permittivity of approximately 1)
electromagneticProperties(emagmodel, "RelativePermittivity", 1);

% Inner region (air gap with relative permittivity of approximately 1)
electromagneticProperties(emagmodel, "RelativePermittivity", 10^9, "Face", 2);

% Specify the electrostatic potential at conductors
% In this case, the inner boundary consists of edges [3 4 5 8] (based on edge labels)
electromagneticBC(emagmodel, "Voltage", 200000, "Edge", [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]);

% Specify the electrostatic potential at the outer boundaries
% In this case, the outer boundary consists of edges [1 2 7 6] (based on edge labels)
electromagneticBC(emagmodel, "Voltage", 0, "Edge", [1 2 3 4]);

% Generate the mesh
%generateMesh(emagmodel);
%generateMesh(emagmodel, 'Hmax', 0.1)

meshEdges = {[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], 0.05};  % Edges selected for increased mesh density  

% Generate the mesh and information about it
generateMesh(emagmodel,'Hmax', 0.2,'Hedge',meshEdges);
%generateMesh(emagmodel);
mesh= emagmodel.Mesh;

% Display number of nodes and elements (cells) in the command window
fprintf('Mesh Information:\n');
fprintf('-----------------\n');
fprintf('Number of nodes: %d\n', size(mesh.Nodes, 2));    % Number of mesh nodes
fprintf('Number of elements: %d\n', size(mesh.Elements, 2)); % Number of mesh elements (cells)

% Check if the mesh is quadratic or linear
numNodesPerElement = size(mesh.Elements, 1);
if numNodesPerElement == 6
    fprintf('Mesh type: Quadratic triangular elements (2D)\n');
elseif numNodesPerElement == 3
    fprintf('Mesh type: Linear triangular elements (2D)\n');
else
    fprintf('Unknown mesh type\n');
end

% Visualize the mesh
figure;
pdeplot(emagmodel, 'Mesh', 'on');  % Turn on mesh display
title('Mesh Visualization');
axis tight;
xlabel('[m]') 
ylabel('[m]')
% Display number of mesh nodes and elements on the top right of the figure
text_nodes = ['Num. of nodes: ', num2str(size(mesh.Nodes,2),'%.f')];  % E_field magnitude at point from contour
text(outerRect(3)+0.2, outerRect(10)-0.5, text_nodes, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude   
text_elements = ['Num. of elem. : ', num2str(size(mesh.Elements,2),'%.f')];  % E_field magnitude at point from contour
text(outerRect(3)+0.2, outerRect(10)-2, text_elements, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude  


% Solve the model
R = solve(emagmodel);

% Calculate the electric potential
u = R.ElectricPotential;

% Get the mesh information
% p: Points matrix containing the X and Y coordinates of mesh points
% e: Edges matrix (not needed)
% t: Triangles matrix (not needed)
[p, e, t] = meshToPet(emagmodel.Mesh);  % Get mesh points, edges, and triangles

% Create an interpolant for the electric potential
F = scatteredInterpolant(p(1,:)', p(2,:)', u, 'linear', 'none');  % X, Y, and potential values

% Define the measurement point
FieldmeasurepointX = sphereX - sphereRadius - 0.001; %front left of the sphere
FieldmeasurepointY = sphereY;                        %front left of the sphere
% FieldmeasurepointX = sphereX;             % Measurement above the sphere
% FieldmeasurepointY = sphereY + sphereRadius + 0.001;  % above the sphere

% Calculate the electric potential at this point
%V_at_point = F(FieldmeasurepointX, FieldmeasurepointY);

% Compute the gradient (electric field) using finite differences
%dx = 1e-5;  % Small distance for numerical differentiation

% Approximate the partial derivatives (gradient)
%dV_dx = (F(FieldmeasurepointX + dx, FieldmeasurepointY) - F(FieldmeasurepointX - dx, FieldmeasurepointY)) / (2 * dx);
%dV_dy = (F(FieldmeasurepointX, FieldmeasurepointY + dx) - F(FieldmeasurepointX, FieldmeasurepointY - dx)) / (2 * dx);

% Electric field components (negative gradient of potential) in that point
%Ex_at_point = -dV_dx;
%Ey_at_point = -dV_dy;

% Calculate the total magnitude of the electric potential in that point
%V_magnitude = sqrt(Ex_at_point^2 + Ey_at_point^2);

% Extract and plot the electric potential
figure;
pdeplot(emagmodel, "XYData", u, "Contour", "on")
axis equal
xlabel('[m]') 
ylabel('[m]') 

% Interpolate the potential at the measurement point from the computed solution
V_at_contour_point = F(FieldmeasurepointX, FieldmeasurepointY);

% Define the text position and display the magnitude of the electric field
text_x = FieldmeasurepointX;  % X-coordinate for the text
text_y = 9;   % Y-coordinate for the text

% Convert the potential and electric field magnitude to strings for display
text_V_contour = ['V = ', num2str(V_at_contour_point, '%.2f'), ' V (contour)'];  % Potential at point from contour
%text_Vmagnitude = ['E = ', num2str(round(V_magnitude), '%d'), ' V/m (calculated)']; % Electric potential magnitude at point

% Display the values on the plot
text(text_x, text_y, text_V_contour, 'FontSize', 12, 'Color', 'k');     % Display potential
%text(text_x, text_y - 0.2, text_Vmagnitude, 'FontSize', 12, 'Color', 'k'); % Display electric field magnitude
hold on;
plot(FieldmeasurepointX, FieldmeasurepointY, 'ko', 'MarkerSize', 6, 'LineWidth', 2); % Mark the measurement point

% Display the results in the command window
disp(['Electric potential at (', num2str(FieldmeasurepointX), ', ', num2str(FieldmeasurepointY), '):']);
disp(['V = ', num2str(V_at_contour_point), 'V (contour)']);
disp(['Electric field at (', num2str(FieldmeasurepointX), ', ', num2str(FieldmeasurepointY), '):']);
% disp(['Ex = ', num2str(Ex_at_point), ' V/m']);
% disp(['Ey = ', num2str(Ey_at_point), ' V/m']);
% disp(['|E|= ', num2str(V_magnitude), ' V/m (calculatued)']);

% Define the range of x-values (m) for which to calculate the electric potential field
%x_values = [0, 1, 1.4, 1.5, 1.6, FieldmeasurepointX, FieldmeasurepointX + sphereRadius * 2 + 0.002, 2.4, 2.5, 2.6, 3, 4, 5, 6, 7,];
x_values = [1,1.2,1.4, FieldmeasurepointX, FieldmeasurepointX+2*sphereRadius+0.001, 2.6, 2.8, 3, 4, 5, 6, 7,];

% Initialize array to store electric potential values (V)
V_at_contour_points = zeros(size(x_values));

% Loop over the x-values to compute the electric potential at each point
for i = 1:length(x_values)
    % Use the current x-value for the field measurement point
    FieldmeasurepointX_current = x_values(i);
    
    % Compute the electric potential at this point using interpolation
    V_at_contour_point = F(FieldmeasurepointX_current, FieldmeasurepointY);
    
    % Store the electric potential value
    V_at_contour_points(i) = V_at_contour_point;
end

% Plot the relationship between distance (m) and electric potential (V)
figure;
plot(x_values, V_at_contour_points, '-o', 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Electric Potential (V)');
title('Electric Potential vs Distance');
grid on;

% Display the computed electric potential values
disp('Electric potential (V) at corresponding distances (m):');
disp(table(x_values', V_at_contour_points', 'VariableNames', {'Distance_m', 'Electric_Potential_V'}));

% Define the range of x-values (m) for which to calculate the electric field
%x_values = [0.1, 1, 1.4, 1.5, 1.6, FieldmeasurepointX, FieldmeasurepointX + sphereRadius * 2 + 0.002, 2.4, 2.5, 2.6, 3, 4, 5, 6, 7,];
x_values = [1, 1.5, 1.6, FieldmeasurepointX-0.014, FieldmeasurepointX+2*sphereRadius+0.031, 2.4,2.5, 3, 3.5, 4, 4.5, 5,];

% Initialize array to store electric field magnitudes (V/m)
E_magnitudes = zeros(size(x_values));

%disp(R.ElectricField)
disp(R.Mesh)

% Extract the node coordinates (x in row 1, y in row 2)
x_nodes = R.Mesh.Nodes(1, :);  % x-coordinates of the nodes
y_nodes = R.Mesh.Nodes(2, :);  % y-coordinates of the nodes

% Loop over the x-values to compute the electric field at each point
for i = 1:length(x_values)
    % Use the current x-value for the field measurement point
    FieldmeasurepointX_current = x_values(i);
    
    % % Compute the electric potential at this point using interpolation
    % %E_at_point = F(FieldmeasurepointX_current, FieldmeasurepointY);
    % 
    % % Compute the gradient (electric field) at this point using finite differences
    % %dV_dx = (F(FieldmeasurepointX_current + dx, FieldmeasurepointY) - F(FieldmeasurepointX_current - dx, FieldmeasurepointY)) / (2 * dx);
    % %dV_dy = (F(FieldmeasurepointX_current, FieldmeasurepointY + dx) - F(FieldmeasurepointX_current, FieldmeasurepointY - dx)) / (2 * dx);
    % 
    % % Electric field components (negative gradient of potential) at this point
    % %Ex_at_point = -dV_dx;
    % %Ey_at_point = -dV_dy;
    % 
    % x_index = find_closest_index(R.ElectricField.X_grid, FieldmeasurepointX_current);
    % y_index = find_closest_index(R.ElectricField.Y_grid, FieldmeasurepointY); 

    % Find the closest node index for the current x and y measurement point
    distances = sqrt((x_nodes - FieldmeasurepointX_current).^2 + (y_nodes - FieldmeasurepointY).^2);
    [~, closest_node_index] = min(distances);  % Find the index of the closest node

    % Extract the electric field components at this point
    % Ex_at_point = R.ElectricField.Ex(x_index, y_index);
    % Ey_at_point = R.ElectricField.Ey(x_index, y_index);
    
    % Extract the electric field components at this node
    Ex_at_point = R.ElectricField.Ex(closest_node_index);
    Ey_at_point = R.ElectricField.Ey(closest_node_index);

    % Calculate the total magnitude of the electric field at this point
    E_magnitude_current = sqrt(Ex_at_point^2 + Ey_at_point^2);
    
    % Store the electric field magnitude
    %E_magnitudes(i) = E_magnitude_current;
    E_magnitudes(i) =  E_magnitude_current;
    
end

j=5;
E_at_contour_point=E_magnitudes(j); % I am selecting the point nr. 2 which is the one in front of the sphere.

% Plot the relationship between distance (m) and electric field magnitude (V/m)
figure;
plot(x_values, E_magnitudes, '-o', 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Electric Field (V/m)');
title('Electric Field Magnitude vs Distance');
grid on;

Ex = R.ElectricField.Ex;
Ey = R.ElectricField.Ey; 
E_magnitude = sqrt(Ex.^2 + Ey.^2);

% Find the maximum magnitude and its location
[max_E_magnitude, max_index] = max(E_magnitude(:));
% Display the maximum value and its location
fprintf('Maximum Electric Field Magnitude: %.4f\n', max_E_magnitude);

figure;
pdeplot(emagmodel, 'XYData', E_magnitude, 'Contour', 'on');
%pdeplot(emagmodel, 'XYData', E_magnitude_log, 'Contour', 'on');
title('Electric Field Magnitude');
axis equal;
xlabel('[m]');
ylabel('[m]');
clim([0 50000]);  % Set the color axis limits
hold on;
pdegplot(emagmodel, 'EdgeLabels', 'off', 'FaceAlpha', 0);
% Convert the electric field magnitude to a string for display
text_E_contour = ['E = ', num2str(E_at_contour_point, '%.2f'), ' V/m (contour)'];  % E_field magnitude at point from contour
% Display the values on the plot
%text(text_x, text_y, text_E_contour, 'FontSize', 12, 'Color', 'k');     % Display electric field magnitude
text(x_values(j), text_y, text_E_contour, 'FontSize', 12, 'Color', 'k');     % Display electric field magnitude
%text(text_x, text_y - 0.2, text_Vmagnitude, 'FontSize', 12, 'Color', 'k'); % Display electric field magnitude
hold on;
%plot(FieldmeasurepointX, FieldmeasurepointY, 'ko', 'MarkerSize', 6, 'LineWidth', 2); % Mark the measurement point
plot(x_values(j), FieldmeasurepointY, 'ko', 'MarkerSize', 6, 'LineWidth', 2); % Mark the measurement point

% Display the computed electric field values
disp('Electric field magnitudes (V/m) at corresponding distances (m):');
disp(table(x_values', E_magnitudes', 'VariableNames', {'Distance_m', 'Electric_Field_Vm'}));