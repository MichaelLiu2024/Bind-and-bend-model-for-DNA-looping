function plotDNAProperties(planarCoordsNoProt, planarCoordsWithProt, type)
%%  Summary
%
%   Inputs:
%       planarCoordsNoProt: Matrix of planar DNA coordinates without
%       protamine.
%       planarCoordsWithProt: Matrix of planar DNA coordinates with
%       protamine.
%       type: Whether the data is experimental or simulated.
%
%   Outputs:
%       None (figures)
%
%   Parent functions:
%       Main
%
%   Child functions:
%       combo_calc
%       centerPlanarCoords
%
%%  Create Figures

%   No protamine figure
    noProt = figure;
    title([type, ' Contours Without Protamine'])
    xlabel('x (nm)')
    ylabel('y (nm)')
    axis([-100, 150, -100, 100])
    axis equal
    box on
    hold on

%   With protamine figure
    withProt = figure;
    title([type, ' Contours With Protamine'])
    xlabel('x (nm)')
    ylabel('y (nm)')
    axis([-100, 150, -100, 100])
    axis equal
    box on
    hold on

%   Cosine figure
    cosine = figure;
    title([type, ' Tangent-Tangent Cosine vs. Arclength'])
    xlabel('Arclength (nm)')
    ylabel('Cosine')
    axis([0, 50, -0.4, 1])
    xticks(0:10:50)
    yticks(-0.4:0.2:1)
    box on
    hold on
%     s = linspace(0, 50);
%     plot(s, exp(-s./120), s, exp(-s./20), s, cos(s./20), 'LineWidth', 5)
    
%   MSD figure
    msd = figure;
    title([type, ' Mean Squared Displacement vs. Arclength'])
    xlabel('Arclength (nm)')
    ylabel('MSD (squared nm)')
    axis([0, 60, 0, 3000])
    xticks(0:10:60)
    yticks(0:1000:3000)
    box on
    hold on
%     s = linspace(0, 60);
%     f = @(x, Lp) 2*Lp*x .* (1 - Lp./x .* (1 - exp(-x./Lp)));
%     plot(s, f(s, 120), s, f(s, 20), s, 2*18^2*(1 - cos(s./18)), 'LineWidth', 5)
    
%   Radius of curvature figure
    radius = figure;
    title([type, ' Distribution of Radii of Curvature'])
    xlabel('Radius of Curvature (nm)')
    ylabel('Probability Density')
    axis([0, 80, 0, 0.04])
    xticks(0:10:80)
    box on
    hold on
    
%   Store curvatures
    noProtCurvatures = NaN(size(planarCoordsNoProt));
    withProtCurvatures = NaN(size(planarCoordsWithProt));

%   Plot colors
    unloopedColor = [0.5, 0.5, 0.5];
    loopedColor = [0.9290, 0.6940, 0.1250];
    
%%  Add Plots
    
%   Plots for no protamine molecules
    for i = 1:size(planarCoordsNoProt, 1)
%       Coordinates
        pc = planarCoordsNoProt(i, :);
        
%       Centered Molecule
        figure(noProt)
        plot(centerPlanarCoords(pc), 'Color', unloopedColor)
        
%       Properties
        dnaProperties = combo_calc(pc);
        
%       Cosine
        figure(cosine)
        plot(dnaProperties.contours, dnaProperties.cosines, 'Color', unloopedColor)
        
%       MSD
        figure(msd)
        plot(dnaProperties.contour, dnaProperties.end2end, 'Color', unloopedColor)
        
%       Store radii
        noProtCurvatures(i, 1:length(pc)) = LineCurvature2D(pc);
    end
    
%   Plots for with protamine molecules
    for i = 1:size(planarCoordsWithProt, 1)
%       Coordinates
        pc = planarCoordsWithProt(i, :);
        
%       Centered Molecule
        figure(withProt)
        plot(centerPlanarCoords(pc), 'Color', loopedColor)
        
%       Properties
        dnaProperties = combo_calc(pc);
        
%       Cosine
        figure(cosine)
        plot(dnaProperties.contours, dnaProperties.cosines, 'Color', loopedColor)
        
%       MSD
        figure(msd)
        plot(dnaProperties.contour, dnaProperties.end2end, 'Color', loopedColor)
        
%       Store radii
        withProtCurvatures(i, 1:length(pc)) = LineCurvature2D(pc);
    end
    
    noProtCurvatures(noProtCurvatures < 1/80) = NaN;
    withProtCurvatures(withProtCurvatures < 1/80) = NaN;
    
    noProtRadii = rmmissing(reshape(1./abs(noProtCurvatures), 1, []));
    withProtRadii = rmmissing(reshape(1./abs(withProtCurvatures), 1, []));
    
%   Radius figure
    figure(radius)
    bandwidth = 4;
    plotKDE(noProtRadii, [0, 80], bandwidth)
    plotKDE(withProtRadii, [0, 80], bandwidth)
    yticks(0:0.01:0.04)
    legend({'No Protamine', 'With Protamine'})
end