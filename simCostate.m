function [lambda,lambdaD] = simCostate(g,lambda0,x,u,t,forward)
% The "simCostate" function integrates the costate in time.
%
% SYNTAX:
%   [lambda,lambdaD] = optimal.simCostate(g,lambda0,x,u,t)
%   [lambda,lambdaD] = optimal.simCostate(g,lambda0,x,u,t,forward)
% 
% INPUTS:
%   g - (function_handle)
%       State dynamics.
%       SYNTAX:
%           lambdaDot = g(x,u,lambda,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           lambda - (n x tn number) Costate.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           lambdaDot - (n x tn number) Costate derivative.
%
%   lambda0 - (n x 1 number)
%       Initial costate.
%
%   x - (n x tn number) 
%       State trajectory.
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
%   lambda - (n x tn number) 
%       Costate trajectory.
%
%   lambdaD - (n x tn number) 
%       Costate time derivative trajectory.
%
% EXAMPLES: TODO: Add examples
%
% NOTES:
%   To simulate backwards in time, all trajectories must be given as
%   trajectories moving backwards in time.
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
n = size(lambda0,1);
lambda = nan(n,tn);
lambdaD = nan(n,tn);

if forward
    lambda(:,1) = lambda0;
    for k = 1:tn-1
        lambdaD(:,k) = g(x(:,k),u(:,k),lambda(:,k),t(k));
        ts = t(k+1) - t(k);
        lambda(:,k+1) = lambda(:,k) + lambdaD(:,k)*ts;
    end
else
    lambda(:,tn) = lambda0;
    for k = tn:-1:2
        lambdaD(:,k) = g(x(:,k),u(:,k-1),lambda(:,k),t(k));
        ts = t(k) - t(k-1);
        lambda(:,k-1) = lambda(:,k) - lambdaD(:,k)*ts;
    end
end
