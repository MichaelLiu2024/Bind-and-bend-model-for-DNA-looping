function simulationData = performSimulation(DNA, AGENT, SIM)
%%  Summary
%
%   Inputs:
%       DNA: A structure that represents the DNA molecule being simulated.
%       Initially, it must contain the following fields:
%           contourLength: The contour length of the DNA.
%           persistenceLength: The persistence length of the DNA.
%           numAgents: The total number of condensing agents that are to
%           bind to the DNA.
%           power: The fraction of the total DNA charge that is felt by
%           other charges in the simulation.
%           kLinear: Parameter for thermal loop closing, defined in the accompanying
%           Biophysical Journal article.
%           kPlanar: Parameter for thermal loop closing, defined in the accompanying
%           Biophysical Journal article.
%       AGENT: A structure that represents the condensing agent being
%       simulated. Initially, it must contain the following fields:
%           length: The total length along the DNA contour that the
%           condensing agent occupies when bound.
%           radius: The effective radius of the condensing agent.
%           power: The fraction of the total condensing agent charge that
%           is felt by other charges in the simulation.
%       SIM: A structure that holds various simulation options. Initially,
%       it must contain the following fields:
%           power: The screening factor p_site, defined in the accompanying
%           Biophysical Journal article.
%           numTrials: The total number of DNA molecules to simulate.
%           displayEvolutions: Whether or not to display the evolution of
%           the DNA as protamine is added.
%           displayMolecules: Whether or not to display the final DNA contours.
%           saveImages: Whether or not to save the final DNA contours as
%           png images.
%
%   Outputs:
%       simulationData: A structure that holds information about the
%       simulated DNA molecules. It contains three fields:
%           startSites: A vector of start sites.
%           loopCircs: A vector of loop circumferences.
%           planarCoords: A matrix of the planar coordinates of the DNA
%           molecules. Note that the ith simulated DNA molecule,
%           represented by planarCoords(i, :), has start site startSites(i)
%           and loop circumference loopCircs(i).
%
%   Parent functions:
%       Main
%
%   Child functions:
%       performTrial
%
%%  Constants

%   Elementary charge
    C.E = 1.602176634e-19;
%   Vacuum permittivity
    C.E0 = 8.8541878128e-12;

%   Boltzmann constant
    C.KB = 1.38064852e-23;
%   Room temperature
    C.T = 300;

%   Maximum distance between the bound condensing agent and the DNA (in nm)
    C.R0 = 1;
%   Radius of the DNA double helix (in nm)
    C.R = 1;
%   Phosphate group radius (in nm)
    C.rPhos = 0.29;

%   Bulk dielectric constant of water
    C.D = 80;
%   Scaling constant
    C.C = 2.674;
%   Half-saturation length of water (in nm)
    C.H = 0.75;

%   Minimum agent-DNA separation = phosphate group radius + condensing agent radius
    C.rMin = C.rPhos + AGENT.radius;

%%  Condensing Agent
    
%   Calculate the maximum sterically allowed bend angle
    AGENT.betaMax = 2/C.R * (C.R0 - C.rMin);
    
%   Calculate the total agent charge
    AGENT.charge = 21 * C.E;

%%  DNA
    
%   Calculate the number of agent binding sites on the DNA
    DNA.numSites = floor(DNA.contourLength / AGENT.length);

%   Adjust the binding site length according to the total number of agent binding sites calculated above
    AGENT.length = DNA.contourLength / DNA.numSites;
    
%   Calculate the total DNA charge per binding site
    DNA.charge = -5.88 * C.E * AGENT.length;
    
%   Store the charges at each DNA vertex
    DNA.vertexCharges = [NaN, repelem(DNA.charge*DNA.power, DNA.numSites), NaN];
    
%   Keep track of the binding sites that are occupied by condensing agent
    DNA.occupiedVertices = false(1, DNA.numSites+2);
    
%%  Trials

%   Create arrays to store the start sites, loop circumferences, and planar
%   DNA contours from the trials
    startSites = NaN(SIM.numTrials, 1);
    loopCircs = NaN(SIM.numTrials, 1);
    planarCoords = NaN(SIM.numTrials, DNA.numSites+2);

%   Perform the trials
    for i = 1:SIM.numTrials
        [startSites(i), loopCircs(i), planarCoords(i,:)] = performTrial(DNA, AGENT, SIM, C, i);
    end
    
%%  Return Values

%   Simulation results are stored in a structure
    simulationData.startSites = startSites;
    simulationData.loopCircs = loopCircs;
    simulationData.planarCoords = planarCoords;
end