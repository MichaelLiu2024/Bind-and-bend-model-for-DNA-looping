function displaySSLC(simulationData, experimentalData, DNA)
%%  Summary
%
%   Inputs:
%       simulationData: Described in performSimulation.
%       experimentalData: Imported experimental looping data.
%       DNA: Described in performSimulation.
%
%   Outputs:
%       None (figures)
%
%   Parent functions:
%       Main
%
%   Child functions:
%       plotKDE
%
%%  Function

%   We fit only the data from single-looped molecules
    simulatedStartSites = simulationData.startSites(isfinite(simulationData.startSites)) / DNA.contourLength;
    simulatedLoopCircs = simulationData.loopCircs(isfinite(simulationData.loopCircs)) / DNA.contourLength;
    
%   Proportion of unlooped and single-looped molecules
    disp("Proportion of unlooped molecules: " + sum(simulationData.startSites == Inf) / length(simulationData.startSites));
    disp("Proportion of single-looped molecules: " + length(simulatedStartSites) / length(simulationData.startSites));

%   Small allowance for KDE support bounds
    epsilon = 1e-10;
    
%   KDE bandwidth
    bandwidth = 0.03;

%   New figure
    figure
    hold on

%   Distribution of simulated start sites
    plotKDE(simulatedStartSites, [-epsilon, 0.5+epsilon], bandwidth)

%   Distribution of experimental start sites
    plotKDE(experimentalData.(['startSites_', num2str(DNA.contourLength)]), [-epsilon, 0.5+epsilon], bandwidth)

%   Title, etc.
    title('Distribution of Fractional Start Sites')
    ylabel('Probability Density')
    xlabel('Fractional Start Site')
	xlim([0, 0.5])
    
    y = round(ylim*2)/2;
    xticks(0:0.1:0.5)
    yticks(y(1):1:y(2))
    
%   New figure
    figure
    hold on
    
%   Distribution of simulated loop circumferences
    plotKDE(simulatedLoopCircs, [-epsilon, 1+epsilon], bandwidth)

%   Distribution of experimental loop circumferences
    plotKDE(experimentalData.(['circumferences_', num2str(DNA.contourLength)]), [-epsilon, 1+epsilon], bandwidth)
    
%   Title, etc.
    title('Distribution of Fractional Loop Circumferences')
    ylabel('Probability Density')
    xlabel('Fractional Loop Circumference')
	xlim([0, 1])
    
    y = round(ylim*2)/2;
    xticks(0:0.2:1)
    yticks(y(1):1:y(2))
end