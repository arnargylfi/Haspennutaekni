clear; close all;
% Create a PDE model for electrostatics
emagmodel = createpde("electromagnetic", "electrostatic");

% Define the outer rectangle
%outerRect = [3; 4; -0.25; 0.25; 0.25; -0.25; 0; 0; 0.5; 0.5];
outerRect = [3; 4; -0.25; 0.25; 0.25; -0.25; 0; 0; 0.25; 0.25];

% Define the inner rectangle (to be subtracted)
%innerRect = [3; 4; -0.025; 0.025; 0.025; -0.025; 0.125; 0.125; 0.5; 0.5];
innerRect = [3; 4; -0.025; 0.025; 0.025; -0.025; 0.125; 0.125; 0.25; 0.25];

% Define the inner circle (to be subtracted)
% Circle is defined by 1, x_center, y_center, radius
innerCirc = [1; 0; 0.125; 0.025; 0; 0; 0; 0; 0; 0];  % Circle with center (0.125, 0.125) and radius 0.025

% Combine the geometry into the geometry description matrix
gd = [outerRect innerRect innerCirc];  % Padding innerCirc with zeros

% Define names for each shape
ns = char('outer', 'innerRect', 'innerCirc');
ns = ns';

% Define the set formula to subtract the inner rectangle and circle from the outer rectangle
sf = 'outer - innerRect - innerCirc';

% Create the decomposed geometry
g = decsg(gd, sf, ns);

% Assign the geometry to the PDE model
geometryFromEdges(emagmodel, g);

% Visualize the geometry with edge labels
figure;
pdegplot(emagmodel, "EdgeLabels", "on");
axis equal;
xlabel('[m]') 
ylabel('[m]') 
title('Geometry with Inner Rectangle and Circle Subtracted');

% Specify the vacuum permittivity value in the SI system of units
emagmodel.VacuumPermittivity = 8.8541878128E-12;

% Specify the relative permittivity of the material
electromagneticProperties(emagmodel, "RelativePermittivity", 1);

% Specify the electrostatic potential at the inner boundary
electromagneticBC(emagmodel, "Voltage", 0, "Edge", [1]);

% Specify the electrostatic potential at the outer boundary
electromagneticBC(emagmodel, "Voltage", 1000, "Edge", [4 5 8 9]);

% Generate the mesh
%generateMesh(emagmodel);

meshEdges = {[4,5,8,9], 0.001};  % Edges selected for increased mesh density  

% Generate the mesh and information about it
generateMesh(emagmodel,'Hmax', 0.01,'Hedge',meshEdges);
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
text(outerRect(3)+0.01, outerRect(10)-0.02, text_nodes, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude   
text_elements = ['Num. of elem. : ', num2str(size(mesh.Elements,2),'%.f')];  % E_field magnitude at point from contour
text(outerRect(3)+0.01, outerRect(10)-0.05, text_elements, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude  


% Solve the model
R = solve(emagmodel);


% Extract and plot the electric potential
figure('Position', [100, 100, 800, 400]); % Set figure size to be wider (800x400)

%figure;
u = R.ElectricPotential;
pdeplot(emagmodel, "XYData", u, "Contour", "on");
axis equal;
xlabel('[m]') 
ylabel('[m]') 
title('Electric Potential Magnitude');
hold on;
pdegplot(emagmodel, 'EdgeLabels', 'off', 'FaceAlpha', 0);
h = findall(gca, 'Type', 'Line');
set(h, 'Color', 'k', 'LineWidth', 1.0); % Set edge color to black and line width


Ex = R.ElectricField.Ex;
Ey = R.ElectricField.Ey; 
E_magnitude = sqrt(Ex.^2 + Ey.^2);
figure;
pdeplot(emagmodel, 'XYData', E_magnitude, 'Contour', 'on');
title('Electric Field Magnitude');
axis equal;
xlabel('[m]');
ylabel('[m]');
hold on;
pdegplot(emagmodel, 'EdgeLabels', 'off', 'FaceAlpha', 0);
h = findall(gca, 'Type', 'Line');
set(h, 'Color', 'k', 'LineWidth', 1.0); % Set edge color to black and line width

