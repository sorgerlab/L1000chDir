function [ dist ] = projectedMetric( ctrl,expm,unitV )
%calculate projected distance from ctrl to expm on chdir
%  input:
%     ctrl: 978*m matrix, replicate matrix of control.
%     expm: 978*k matrix, replicate matrix of an expm.
%     unitV: 978*1 vector, characteristic direction for the expm.

ctrlCount = size(ctrl,2);
expmCount = size(expm,2);

ctrlProjection = dot(ctrl,repmat(unitV,1,ctrlCount));
expmProjection = dot(expm,repmat(unitV,1,expmCount));

dist = abs(mean(expmProjection) - mean(ctrlProjection));

end

