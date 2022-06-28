function [Data_good] = qc_filter(Data,variables,qc_flags)

% DESCRIPTION:
%   This function generates a new data structure composed of chosen variables
%   based on provided QC flag values
%
% INPUTS:
%   Data : struct that must contain the PRES field and the given
%               variables (_ADJUSTED fields are used if available) 
%   variables: cell array with names of the measured fields (e.g., BBP700)
%   qc_flags: numerical array of QC flag values
%
% AUTHORS: 
%   N. Buzby (UW)
%
% DATE: October 6, 2021

Data_good = struct('CYCLE_NUMBER',Data.CYCLE_NUMBER,'TIME',Data.TIME,...
    'LATITUDE',Data.LATITUDE,'LONGITUDE',Data.LONGITUDE);

for v=1:length(variables)
    if isfield(Data,[variables{v}, '_ADJUSTED'])
        field = Data.([variables{v}, '_ADJUSTED']);
        idx=ismember(Data.([variables{v}, '_ADJUSTED_QC']),qc_flags);
    else
        warning(['adjusted values for %s for are not available,',...
                ' using raw values instead'], variables{v})
        field = Data.(variables{v});
        idx=ismember(Data.([variables{v}, '_QC']),qc_flags);
    end
    Data_good.(variables{v}) = field.*idx;
    Data_good.(variables{v})(Data_good.(variables{v}) == 0) = NaN;
end