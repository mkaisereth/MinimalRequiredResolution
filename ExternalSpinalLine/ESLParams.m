classdef ESLParams
    properties
        DownSample
        AbsoluteSampling
        SmoothRadius
        DiscreteKernel
        DiscreteSteps
        DidShift0_1
        DidConvertFromMToMM
        ImgForegroundMin
        ImgForegroundMax
        ManualPcRoi
        LimitZResolution
        FrequencyFilter
        % not used anymore? 09.09.2022, use real xAvg
        %FrequencyFilterSamplingPeriod
        FilePath
    end
end