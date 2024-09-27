clear; close all;
% Create a PDE model for electrostatics
emagmodel = createpde("electromagnetic","electrostatic");

% Manually define the geometry using "outerRect" command. 3 specifies 2D
% shape, and 4 specifies the number of corner points (vertices), so a
% rectangle.
% Outer rectangle (3x6) and inner rectangle (0.3x1.5)
% Outer rectangle (larger frame): first four are x-coordinates in counterclockwise order
% The next four 1 are the corresponding y-coordinates.
%outerRect = [3; 4; -1.5; 1.5; 1.5; -1.5; -3 ; -3; 3; 3];
outerRect = [3; 4; -0.0015; 0.0015; 0.0015; -0.0015; -0.003 ; -0.003; 0.003; 0.003];

% Inner rectangle (the smaller gap) 
% First four are x-coordinates in counterclockwise order
% The next four are the corresponding y-coordinates. 
%innerRect = [3; 4; -0.3; 0.3; 0.3; -0.3; -1.5; -1.5; 1.5; 1.5]; 
innerRect = [3; 4; -0.0003; 0.0003; 0.0003; -0.0003; -0.0015; -0.0015; 0.0015; 0.0015]; 

% Combine both rectangles into the geometry matrix
gd = [outerRect innerRect];

% Set up the name-space for the shapes
ns = char('outer', 'inner');
ns = ns';

% Specify the set formula (how the shapes are combined)
% 'outer - inner' means the outer rectangle with the inner rectangle removed (gap)
sf = 'outer + inner';

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
%electromagneticProperties(emagmodel, "RelativePermittivity", 4);

% Outer region (dielectric material with relative permittivity of 4)
electromagneticProperties(emagmodel, "RelativePermittivity", 4, "Face",1);
% Inner region (air gap with relative permittivity of approximately 1)
electromagneticProperties(emagmodel, "RelativePermittivity", 1, "Face", 2);

% Specify the electrostatic potential at the inner boundary
% In this case, the inner boundary consists of edges [3 4 5 8] (based on edge labels)
electromagneticBC(emagmodel, "Voltage", 0, "Edge", [6]);

% Specify the electrostatic potential at the outer boundary
% In this case, the outer boundary consists of edges [1 2 7 6] (based on edge labels)
electromagneticBC(emagmodel, "Voltage", 14000, "Edge", [1]);

meshEdges = {[3,4,5,8], 0.00002};  % Edges selected for increased mesh density  

% Generate the mesh and information about it
generateMesh(emagmodel,'Hmax', 0.00005,'Hedge',meshEdges);
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
text(outerRect(3)+0.0001, outerRect(10)-0.0002, text_nodes, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude   
text_elements = ['Num. of elem. : ', num2str(size(mesh.Elements,2),'%.f')];  % E_field magnitude at point from contour
text(outerRect(3)+0.0001, outerRect(10)-0.0005, text_elements, 'FontSize', 12, 'Color', 'k','BackgroundColor', 'white','Margin', 0.5, 'EdgeColor', 'none');     % Display electric field magnitude  

% Solve the model
R = solve(emagmodel);

% Extract and plot the electric potential
figure;
u = R.ElectricPotential;
pdeplot(emagmodel, "XYData", u, "Contour", "on")
title('Electric Potential Magnitude');
axis equal
xlabel('[m]') 
ylabel('[m]') 
hold on;
pdegplot(emagmodel, 'EdgeLabels', 'off', 'FaceAlpha', 0);
h = findall(gca, 'Type', 'Line');
set(h, 'Color', 'k', 'LineWidth', 1.0); % Set edge color to black and line width

Ex = R.ElectricField.Ex;
Ey = R.ElectricField.Ey; 
E_magnitude = sqrt(Ex.^2 + Ey.^2);  % Calculating the electric field for all points

% For three visible measurement points, define the range of x-values
x_values = [-0.0009, 0, 0.0009,0];
y_values = [0, 0, 0,-0.0028];

% Initialize array to store electric field magnitudes for measurement points for display 
E_magnitudes = zeros(size(x_values));

% Extract the node coordinates (x in row 1, y in row 2)
x_nodes = R.Mesh.Nodes(1, :);  % x-coordinates of the nodes
y_nodes = R.Mesh.Nodes(2, :);  % y-coordinates of the nodes

% Plot the electric field contours
figure;
pdeplot(emagmodel, 'XYData', E_magnitude, 'Contour', 'on');
title('Electric Field Magnitude');
axis equal;
xlabel('[m]');
ylabel('[m]');
hold on;    % Hold on to add all the extra features to the figure
% Set the edge color (for all lines) to black so that the gap region is visible
pdegplot(emagmodel, 'EdgeLabels', 'off', 'FaceAlpha', 0);
h = findall(gca, 'Type', 'Line');
set(h, 'Color', 'k', 'LineWidth', 1.0); % Set edge color to black and line width
%currentAxis = axis;

% Loop over the x-values to compute the electric field at each point
for i = 1:length(x_values)
  
    % Find the closest mesh node point index for the current x and y measurement point
    distances = sqrt((x_nodes - x_values(i)).^2 + (y_nodes - y_values(i)).^2);
    [~, closest_node_index] = min(distances);  % Find the index of the closest node

    % Extract the electric field component values at this node
    Ex_at_point = R.ElectricField.Ex(closest_node_index);
    Ey_at_point = R.ElectricField.Ey(closest_node_index);

    % Calculate the total magnitude of the electric field at this point
    E_magnitudes(i) = sqrt(Ex_at_point^2 + Ey_at_point^2);
    
    % Set the position coordinates for the E annotation value (at the measurement point)
    text_x_E=x_values(i)+0.00005;
    text_y_E=y_values(i);
    % Set the position coordinates for the E magnitude value (at the top)
    text_x_Emagnitude=outerRect(3)+0.0001;
    text_y_Emagnitude=outerRect(10)-0.0002-(i-1)*0.0002; % Decrement the y position value for each case 

    % Convert the electric field values to a string on the figure
    % starting with the point of measurement
    plot(x_values(i), y_values(i), 'ko', 'MarkerSize', 3, 'LineWidth', 2); % Mark the measurement point
    % E annotation and the associated number next to the measurement point
    text_E = ['E', num2str(i)];  % E_field magnitude at point from contour
    text(text_x_E, text_y_E, text_E, 'FontSize', 12, 'Color', 'k');     
    % E magnitude value displayed on the top of the figure
    text_E_magnitude = ['E',num2str(i),' = ', num2str(E_magnitudes(i), '%.f'), ' V/m'];  % E_field magnitude at point from contour
    text(text_x_Emagnitude, text_y_Emagnitude, text_E_magnitude, 'FontSize', 12, 'Color', 'k');     % Display electric field magnitude   
   
end
%axis(currentAxis);
hold off;
