function [startSite, loopCirc, planarCoords] = performTrial(DNA, AGENT, SIM, C, trialNumber)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%       AGENT: Described in parent function.
%       SIM: Described in parent function.
%       C: Physical constants. Created in parent function.
%       trialNumber: Keeps track of the trial number for naming images.
%
%   Outputs:
%       startSite: The simulated molecule's start site.
%       loopCirc: The simulated molecule's loop circumference.
%       planarCoords: The planar coordinates of each endpoint and binding
%       site on the DNA.
%
%   Parent functions:
%       performSimulation
%
%   Child functions:
%       computePlanarCoords
%       determineOrientation
%       computeBindingProbabilities
%       computeCrossoverData
%       computeClosenessData
%       computeStartSiteLoopCirc
%       plotPlanarCoords
%
%%  Initialize the DNA state

%   Generate a smooth worm-like chain representation of the DNA contour
    DNA.linearCoords = [0, linspace(AGENT.length/2, DNA.contourLength - AGENT.length/2, DNA.numSites), DNA.contourLength];
    DNA.vertexAngles = [0, normrnd(0, sqrt(2 * AGENT.length/DNA.persistenceLength), 1, DNA.numSites), NaN];    
    DNA.planarCoords = computePlanarCoords(DNA);
    
%   Apply a cubic smoothing spline to the DNA contour
    DNA.planarCoords = csaps(DNA.linearCoords, DNA.planarCoords, 0.1, DNA.linearCoords);
    DNA.linearCoords = cumsum([0, abs(diff(DNA.planarCoords))]);
    DNA.vertexAngles = diff([0, unwrap(angle(diff(DNA.planarCoords))), NaN]);
    
%   For every 10th molecule of DNA on average, add half as many condensing agents
    if(rand < 0.1)
        DNA.numAgents = round(DNA.numAgents/2);
    end

%   Assign a random DNA orientation to start
    DNA.orientation = randsample([-1, 1], 1);
    
%   Make the condensing agent bending in phase with the existing DNA bend orientation
    DNA.orientation = determineOrientation(DNA);

%%  Simulate the agent binding events

%   There are numAgents binding events total
    for i = 1:DNA.numAgents
%       Compute the relative probability of agent binding at each binding site
        bindingProbabilities = computeBindingProbabilities(DNA, AGENT, SIM, C);
        
%       Plot the DNA contour and the probability distribution for agent binding
        if(SIM.displayEvolutions)
            plotPlanarCoords(DNA, []);

            figure
            plot(DNA.linearCoords, bindingProbabilities);

            if(SIM.saveImages)
                saveas(gcf, ['C:\Users\carterlaboratory\Desktop\Contours\', 'Contour_', num2str(trialNumber), '_', num2str(i), '.svg'])
            end
        end
        
%       Choose the binding site at which the next agent binds
        nextSite = randsample(DNA.numSites+2, 1, true, bindingProbabilities);
        
%       Place a agent at the next binding site
        DNA.occupiedVertices(nextSite) = 1;
        
%       Update the vertex charges
        DNA.vertexCharges(nextSite) = DNA.vertexCharges(nextSite) + AGENT.charge*AGENT.power;

%       Update the DNA contour to include agent bending
        additionalBend = DNA.orientation*AGENT.betaMax - min(AGENT.betaMax, max(-AGENT.betaMax, DNA.vertexAngles(nextSite)));
        DNA.vertexAngles(nextSite) = DNA.vertexAngles(nextSite) + additionalBend;
        DNA.planarCoords = computePlanarCoords(DNA);
    end

%%  Calculate the DNA crossover points or closest DNA contacts

%   Find crossover points
    crossoverData = computeCrossoverData(DNA.planarCoords);
    
%   Store the closest DNA contacts
    closenessData = [];

%   Determine if the DNA molecule is multi-looped, single-looped, or unlooped
    if(size(crossoverData, 1) > 1)
%       The DNA is multi-looped
        startSite = -Inf;
        loopCirc = -Inf;
    elseif(size(crossoverData, 1) == 1)
%       The DNA is single-looped via a crossover event
        [startSite, loopCirc] = computeStartSiteLoopCirc(DNA, [crossoverData, crossoverData(2)]);
    else
%       Calculate the closest DNA contacts
        closenessData = computeClosenessData(DNA);
        
        if(~isempty(closenessData))
%           The DNA is single-looped via thermal loop closing
            [startSite, loopCirc] = computeStartSiteLoopCirc(DNA, repelem(closenessData,2));
        else
%           The DNA is unlooped
            startSite = Inf;
            loopCirc = Inf;
        end
    end

%   Store the planar coordinates
    planarCoords = DNA.planarCoords;
    
%   Plot the DNA contour, agent locations, and crossover/closeness
%   points, if applicable
    if(SIM.displayMolecules)
        plotPlanarCoords(DNA, [crossoverData(:, 2), closenessData]);

        if(SIM.saveImages)
            saveas(gcf, ['C:\Users\carterlaboratory\Desktop\Contours\', 'Contour_', num2str(trialNumber), '.svg'])
        end

%       Save images if requested
        if(SIM.saveImages && size(crossoverData, 1) > 1)
            saveas(gcf, ['C:\Users\micha\OneDrive\Documents\MATLAB\Thermal Fluctuation Model 7_X_2021\Simulated Data\Multiloop Images 8_12_2021\Image ', num2str(trialNumber), '.png'])
        end
    end
end