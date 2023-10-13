function plotKDE(data, support, bandwidth)
%%  Summary
%
%   Inputs:
%       data: Data to which a KDE will be fit.
%       support: Support bounds for the KDE.
%       bandwidth: Smoothing parameter.
%
%   Outputs:
%       None (figures)
%
%   Parent functions:
%       displaySSLC
%
%   Child functions:
%       None
%
%%  Create and plot KDE

    [pdfValues, xValues] = ksdensity(data, 'Support', support, 'BoundaryCorrection', 'Reflection', 'Bandwidth', bandwidth); 
    area(xValues, pdfValues, 'FaceAlpha', 0.5)
    box on
    set(gca, 'linewidth', 2)
end

