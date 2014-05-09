function [x,xD] = simState(f,x0,u,t,forward)
% The "simState" function integrates the state in time.
%
% SYNTAX:
%   [x,xD] = optimal.simState(f,x0,u,t)
%   [x,xD] = optimal.simState(f,x0,u,t,forward)
% 
% INPUTS:
%   f - (function_handle)
%       State dynamics.
%       SYNTAX:
%           xDot = f(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           xDot - (n x tn number) State derivative.
%
%   x0 - (n x 1 number) 
%       Initial state.
%
%   u - (m x tn-1 number)
%       Input trajectory.
%
%   t - (1 x tn number)
%       Time trajectory.
%
%   forward - (1 x l logical) [true]
%       If true the simulates forward in time, if false simulates backward
%       in time.
% 
% OUTPUTS:
%   x - (n x tn number) 
%       State trajectory.
%
%   xD - (n x tn number) 
%       State time derivative trajectory.
%
% EXAMPLES: TODO: Add examples
%
% NOTES:
%
% NECESSARY FILES:
%
% SEE ALSO:
%    optimal.simCostate
%
% AUTHOR:
%    Rowland O'Flaherty (www.rowlandoflaherty.com)
%
% VERSION: 
%   Created 02-MAY-2014
%-------------------------------------------------------------------------------

%% Apply default values
if nargin < 5, forward = true; end

%% Simulate
tn = numel(t);
n = size(x0,1);
x = nan(n,tn);
xD = nan(n,tn-1);

if forward
    x(:,1) = x0;
    for k = 1:tn-1
        xD(:,k) = f(x(:,k),u(:,k),t(k));
        ts = t(k+1) - t(k);
        x(:,k+1) = x(:,k) + xD(:,k)*ts;
    end
else
    x(:,tn) = x0;
    for k = 1:tn-1
        xD(:,k) = f(x(:,k),u(:,k),t(k));
        ts = t(k) - t(k-1);
        x(:,k+1) = x(:,k) - xD(:,k)*ts;
    end
end

end
