function sweepData = performParameterSweep(DNA, AGENT, SIM, parameterType, parameterName, parameterValues, experimentalData_217_398_1023)
%%  Summary

%%  Parameter Sweep

    for i = 1:length(parameterValues)
        
        switch parameterType
            case 'DNA'
                DNA.(parameterName) = parameterValues(i);
            case 'AGENT'
                AGENT.(parameterName) = parameterValues(i);
            case 'SIM'
                SIM.(parameterName) = parameterValues(i);
        end
        
        sweepData.(['simulationData', num2str(i)]) = performSimulation(DNA, AGENT, SIM);
        
        displaySSLC(sweepData.(['simulationData', num2str(i)]), experimentalData_217_398_1023, DNA)
    end
end

