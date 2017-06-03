% SCRIPT to check lobe area computed in lober using the MC type computation
% of area

% Load and plot the lobes
lobe1 = load('lobeb.dat.00002');
lobe2 = load('lobeb.dat.00004');

% Generate random points inside the box that bounds the lobe1
xMin = min(lobe1(:,1));
xMax = max(lobe1(:,1));

yMin = min(lobe1(:,2));
yMax = max(lobe1(:,2));

testPts = 1e6;
xPts = xMin + (xMax - xMin)*rand(testPts,1);
yPts = yMin + (yMax - yMin)*rand(testPts,1);

in = inpolygon(xPts,yPts,lobe1(:,1),lobe1(:,2));

% Number of pts that landed inside the polygon
numPtsIn = length(xPts(in));
fracIn = numPtsIn/testPts;

% Area of the polygon = fraction that landed inside * area of the box([xMin, xMax] x [yMin, yMax]) 
areaLobe = fracIn * (yMax - yMin)*(xMax - xMin)

% Plotting
plot(lobe1(:,1),lobe1(:,2),'-r')
hold on
plot(lobe2(:,1),lobe2(:,2),'-b')
% plot(xPts, yPts, '.k')


plot(lobe1(:,1),lobe1(:,2),xPts(in),yPts(in),'.r',xPts(~in),yPts(~in),'.k')









