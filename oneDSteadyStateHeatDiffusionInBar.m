close all; clear; clc

% User Inputs 

barLength = 5; % m
crossSectionalAreaofBar = 0.1; % m^2
numberOfCells = 10; % dimensionless

thermalConductivity = 100; % W/m-K
heatSourcePerVolume = 1000; % W/m^3

temperatureWestBoundary = 300; % deg C
temperatureEastBoundary = 250; % deg C

printSetup = true;
printSolution = true;
plotOutput = true;


% Create the mesh of cells

fprintf(' Creating Mesh ...\n');
fprintf('------------------------------------------------\n');

xCoordinatesOfFaces = linspace(0, barLength, numberOfCells+1);

xCoordinatesOfCentroids = 0.5 * (xCoordinatesOfFaces(1:end-1) + xCoordinatesOfFaces(2:end));

distanceBetweenCentroids = xCoordinatesOfCentroids(2:end) - xCoordinatesOfCentroids(1:end-1);    
distanceWestBoundary = 2 * ((xCoordinatesOfCentroids(1)) - (xCoordinatesOfFaces(1)));
distanceEastBoundary = 2 * ((xCoordinatesOfFaces(end)) - (xCoordinatesOfCentroids(end)));
distanceBetweenCentroids = [distanceWestBoundary, distanceBetweenCentroids, distanceEastBoundary];

cellLength = xCoordinatesOfFaces(2:end) - xCoordinatesOfFaces(1:end-1);
cellVolumes = cellLength * crossSectionalAreaofBar;

% Calculate the Matrix Coefficients

fprintf(' Calculating Matrix Coefficients\n');
fprintf('------------------------------------------------\n');

DA = crossSectionalAreaofBar * thermalConductivity ./ distanceBetweenCentroids;
aW = [0,DA(2:end-1)];
aE = [DA(2:end-1),0];

Sp = [-2*DA(1), zeros(1, numberOfCells-2), -2*DA(end)];
aP = aW + aE - Sp;

Su =[2*DA(1)*temperatureWestBoundary + heatSourcePerVolume*cellVolumes(1), ...
    heatSourcePerVolume*cellVolumes(2:end-1), ...
    2*DA(end)*temperatureEastBoundary + heatSourcePerVolume*cellVolumes(end)];

% Create the matrices

fprintf(' Assembling Matrices\n');
fprintf('------------------------------------------------\n');

Amatrix = diag(-aW(2:end), -1) + diag(aP, 0) + diag(-aE(1:end-1), 1);

BVector = Su';

% Print the setup

if printSetup
    fprintf(' Matrix Coefficients\n');
    fprintf('------------------------------------------------\n');
    disp(table(aW', aP', aE', Su', 'VariableNames', {'aW', 'aP', 'aE', 'Su'}));
    fprintf('------------------------------------------------\n');
    disp(table(Amatrix, BVector, 'VariableNames', {'Amatrix', 'BVector'}));
    fprintf('------------------------------------------------\n');
end

% Solve the matrices

fprintf(' Solving ...\n')
fprintf('------------------------------------------------\n')

temperatureVector = Amatrix \ BVector;

fprintf(' Equations Solved.\n')
fprintf('------------------------------------------------\n')

% Print the solution

if printSolution
    fprintf(' Solution\n');
    fprintf('------------------------------------------------\n');
    disp(table(temperatureVector', 'VariableNames', {'Temperature'}));
end

% Plot the solution

if plotOutput
    figure;
    xPlotting = [xCoordinatesOfFaces(1), xCoordinatesOfCentroids, xCoordinatesOfFaces(end)];
    temperatureVectorPlotting = [temperatureWestBoundary; temperatureVector; temperatureEastBoundary];
    plot(xPlotting, temperatureVectorPlotting, 'o-');
    grid on;
    xlabel('Position (m)');
    ylabel('Temperature (C)');
    title('Temperature Distribution');
end

% Code Ends Here
