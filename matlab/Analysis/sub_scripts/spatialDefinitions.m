%% Spatial definitions
sSize = 2^k-1; % Number of samples; odd number so that vortex gets
               % centered (spatial size); Spatial size. ref: 2^k-1
SpatialSupport = SpatialSupport/2; % Half support of the SLM window in cm
spaceVector = -SpatialSupport:2*SpatialSupport/(sSize-1):SpatialSupport; 
% Symmetric space
[X,Y] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian coordinates

%% Polar coordinates with a shift of the mask
if shiftBool == 1
 shiftCart = SpatialSupport*shiftCart/100; % Percentage w.r.t the half size 
 shiftX = shiftCart(2); % Cartesian shift in x
 shiftY = shiftCart(1); % Cartesian shift in y
else
 shiftX = 0; shiftY = 0; % Shift deactivated   
end
[phi,r] = cart2pol(X-shiftX,Y+shiftY); % Polar coordinates with an added
                                       % shift. The signs compensate the 
                                       % normal cartesian convention for 
                                       % displacing the phase mask
x = spaceVector; % Cartesian x-vector
y = x; % Cartesian y-vector: square grid