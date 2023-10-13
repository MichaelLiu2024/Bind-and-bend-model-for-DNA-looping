function bindingProbabilities = computeBindingProbabilities(DNA, AGENT, SIM, C)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%       AGENT: Described in parent function.
%       SIM: Described in parent function.
%       C: Described in parent function.
%
%   Outputs:
%       bindingProbabilities: A vector containing the relative probability
%       of condensing agent binding at each binding site.
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       uElec
%
%%  Compute the total energy change due to agent binding (uBind) and DNA bending (uBend) at each unoccupied binding site
    
%   Store the binding and bending energies
    uBind = NaN(1, DNA.numSites+2);
    uBend = NaN(1, DNA.numSites+2);
    
%   Numerator of Coulomb's Law for the j != i binding sites
    qAgent_qJ = AGENT.charge * DNA.vertexCharges;
    
%   Calculate the energy for each unoccupied binding site
    for i = find(~DNA.occupiedVertices)
%%      Condensing agent cannot bind at the DNA endpoints
        if(i == 1 || i == length(DNA.occupiedVertices))
            continue
        end
        
%%      uBind
        
%       Angle at the ith binding site
        beta = abs(DNA.vertexAngles(i));

%       Distance between sites for the j != i binding sites (min value is 3.4nm)
        r = max(AGENT.length/2, abs(DNA.planarCoords - DNA.planarCoords(i)));
        
%       Electrostatic energy contributions from the j != i binding sites (uJ)
        uJ = uElec(r, qAgent_qJ, C);
        
%       Effective agent-DNA separation for the j = i binding site
        rI = max(C.rMin, C.R0 - C.R/2 * beta);
        
%       Numerator of Coulomb's Law for the j = i binding site
        qAgent_qDNA = AGENT.charge * DNA.charge * SIM.power;
        
%       Electrostatic energy contribution from the j = i binding site (uI)
        uJ(i) = uElec(rI, qAgent_qDNA, C);
        
%       Include the entropic energy cost from the j = i binding site
        uJ(i) = uJ(i) + C.KB*C.T/2 * log(-uJ(i));
        
%       Remove j = i binding site energy (debug)
%       uJ(i) = 0;
        
%       Total stabilization energy
        uBind(i) = sum(uJ(2:end-1));
        
%       Remove j != i binding site energy (debug)
%       uBind(i) = uJ(i);

%%      uBend

%       Bend the DNA only if the existing thermal bend is smaller than
%       betaMax
        if(beta > AGENT.betaMax)
            uBend(i) = 0;
        else
            uElas = C.KB*C.T * DNA.persistenceLength/(4*AGENT.length) * (AGENT.betaMax^2 - beta^2);
            finalUElec = uElec(C.rMin, qAgent_qDNA, C);
            initialUElec = uElec(rI, qAgent_qDNA, C);
            
            uBend(i) = uElas + finalUElec - initialUElec + C.KB*C.T/2 * log(finalUElec/initialUElec);
            
%           Remove bending energy (debug)
%           uBend(i) = 0;
        end
    end

%%  Shift the total energies by a constant to avoid exponentiating large numbers in the Boltzmann distribution

%   Total energy in KB*T units
    uTotal = (uBind + uBend) / (C.KB*C.T);
    
%   Shift amount
    shift = quantile(uTotal, 0.9) - diff(quantile(uTotal, [0.1, 0.9]))/2;
    
%   Calculate relative binding probabilities
	bindingProbabilities = exp(-(uTotal - shift));

%   Invalid binding sites have zero binding probability
    bindingProbabilities(isnan(bindingProbabilities)) = 0;
    
%%  Debug plots

%   figure
%   plot(uBind/C.KB/C.T)
%   figure
%   plot(uBend/C.KB/C.T)
%   figure
%   plot(uTotal)
end