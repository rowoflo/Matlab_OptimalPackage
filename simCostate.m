function [lambda,lambdaD] = simCostate(g,lambda0,x,u,t,varargin)
% The "simCostate" function integrates the costate
%
% SYNTAX:
%   [lambda,lambdaD] = optimal.simCostate(g,lambda0,x,u,t)
%   [lambda,lambdaD] = optimal.simCostate(g,lambda0,x,u,t,'PropertyName',PropertyValue,...)
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
% PROPERTIES:
%   'direction' - ('forward' or 'backward') ['forward']
%       The direction of simulation from the intial state.
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

%% Check Inputs

% Check number of inputs
narginchk(5,inf)

% Check input arguments for errors
assert(isa(g,'function_handle'),...
    'optimal:simCostate:g',...
    'Input argument "g" must be a function handle.')

assert(isnumeric(t) && isreal(t) && isvector(t),...
    'optimal:simCostate:t',...
    'Input argument "t" must be a vector of real numbers.')
t = t(:)';
tn = numel(t);

assert(isnumeric(lambda0) && isvector(lambda0),...
    'optimal:simCostate:lambda0',...
    'Input argument "lambda0" must be a vector.')
lambda0 = lambda0(:);
n = size(lambda0,1);

assert(isnumeric(x) && isequal(size(x),[n,tn]),...
    'optimal:simCostate:x',...
    'Input argument "x" must be a %d x %d matrix.',n,tn)

assert(isnumeric(u) && numel(size(u)) == 2 && size(u,2) == tn-1,...
    'optimal:simCostate:u',...
    'Input argument "u" must be a matrix with a length of %d.',tn-1)

% Get and check properties
propargin = size(varargin,2);

assert(mod(propargin,2) == 0,...
    'optimal:simCostate:properties',...
    'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')

propStrs = varargin(1:2:propargin);
propValues = varargin(2:2:propargin);

for iParam = 1:propargin/2
    switch lower(propStrs{iParam})
        case lower('direction')
            direction = propValues{iParam};
        otherwise
            error('optimal:simCostate:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('direction','var'), direction = 'forward'; end

% Check property values for errors TODO: Add property error checks
assert(ismember(direction,{'forward','backward'}),...
    'optimal:simCostate:direction',...
    'Property "direction" must be either ''forward'' or ''backward''.')

%% Simulate
lambda = nan(n,tn);
lambdaD = nan(n,tn);
lambda(:,tn) = lambda0;

%% Simulate backward
for k = tn:-1:2
    lambdaD(:,k) = g(x(:,k),u(:,k-1),lambda(:,k),t(k));
    ts = t(k) - t(k-1);
    lambda(:,k-1) = lambda(:,k) - lambdaD(:,k)*ts;
end


% if strcmp(direction,'forward')
%     [lambda,lambdaD] = optimal.simulate(@(lambda_,t_) g(xFFunc(t_,x,t),uFFunc(t_,u,t),lambda_,t_),lambda0,t);
% %     [lambda,lambdaD] = optimal.simulate(@(lambda_,k_) g(x(k_),u(k_),lambda_,t_),lambda0,t);
% else
%     [lambda,lambdaD] = optimal.simulate(@(lambda_,t_) g(xBFunc(t_,x,t),uBFunc(t_,u,t),lambda_,t_),lambda0,t(end:-1:1));
% %     [lambda,lambdaD] = optimal.simulate(@(lambda_,k_) g(x(k_),u(k_),lambda_,t_),lambda0,t);
%     lambda = lambda(:,end:-1:1);
%     lambdaD = lambdaD(:,end:-1:1);
% end
% 
% end
% 
% function x = xFFunc(t,xTraj,tTraj)
% k = find(t >= tTraj,1,'last');
% x = xTraj(k);
% end
% 
% function x = xBFunc(t,xTraj,tTraj)
% k = find(t > tTraj,1,'last');
% x = xTraj(k);
% end
% 
% function u = uFFunc(t,uTraj,tTraj)
% k = find(t >= tTraj,1,'last');
% u = uTraj(k);
% end
% 
% function u = uBFunc(t,uTraj,tTraj)
% k = find(t > tTraj,1,'last');
% u = uTraj(k);
% end
