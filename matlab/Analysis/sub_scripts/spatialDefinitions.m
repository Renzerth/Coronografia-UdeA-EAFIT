%% Spatial definitions
sSize = 2^k-1; % Number of samples; odd number so that vortex gets
               % centered (spatial size); Spatial size. ref: 2^k-1
SpatialSupport = SpatialSupport/2; % Half support of the SLM window in cm
spaceVector = -SpatialSupport:2*SpatialSupport/(sSize-1):SpatialSupport; 
% Symmetric space
[X,Y] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian coordinates
[phi,r] = cart2pol(X,Y); % Polar coordinates
x = spaceVector; % Cartesian x-vector
y = x; % Cartesian y-vector