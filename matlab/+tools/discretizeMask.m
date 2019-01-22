function discreteMask = discretizeMask(mask,levels)
% This function performs a discretization of a given mask where the 
% number of levels and values are given as an array of real values.
% Right now we are assuming that the values must be between -pi..pi.

discreteMask = mask;
[n, m] = size(discreteMask);

for i = 1:n,
    for j = 1:m,
        found = 1;
        z = 1;
        while found ~= 0,
            if mask(i, j) <= levels(z+1);
                discreteMask(i,j) = levels(z);
                found = 0;
            else
                z = z + 1;
            end
        end
    end
end

end