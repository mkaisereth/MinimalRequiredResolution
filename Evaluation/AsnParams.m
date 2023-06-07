% class for add sine noise parameters
classdef AsnParams
    properties
        % the number of peaks to add within height and width
        NumOfPeaks
        % the amplitude of the peaks [in % to zMax-zMin in [0,1]]
        Amplitude
        % additional random noise [in % to zMax-zMin in [0,1]]
        RandomNoise
    end
end