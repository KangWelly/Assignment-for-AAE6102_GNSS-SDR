function showAcqResult(Acquired, signal)
%Prints the Acq. results
%
%showAcqResult(Acquired, signal)
%
%   Inputs:
%       Acquired - Acq. results
%       signal - signal settings

%--------------------------------------------------------------------------

fprintf('Acquisition results\n');
fprintf(  '|   PRN   |   Doppler  | Code Offset | \n');
for channelNr = 1 : length(Acquired.sv)
        fprintf('|    %2d   |   %5.0f    |    %6d   |\n', ...
                Acquired.sv(channelNr), ...
                Acquired.fineFreq(channelNr) - signal.IF, ...
                Acquired.codedelay(channelNr) ...
            );
end
fprintf('\n\n');
