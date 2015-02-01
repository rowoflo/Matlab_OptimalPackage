function [x,u,xDot] = simulate(f,g,t,x0,xm,xM,varargin)
% The "simulate" function integrates the state in time.
%
% SYNTAX:
%   [x,u,xDot] = optimal.simulate(f,g,x0,t)
%   [x,u,xDot] = optimal.simulate(f,g,x0,t,xm,xM)
%   [x,u,xDot] = optimal.simulate(f,g,x0,t,xm,xM,'PropertyName',PropertyValue,...)
% 
% INPUTS:
%   f - (function_handle)
%       State dynamics.
%       SYNTAX:
%           xDot = f(x,u,t,k);
%       INPUTS:
%           x - (n x 1 number) State.
%           u - (m x 1 number) Input.
%           t - (1 x 1 number) Time.
%           k - (1 x 1 integer) Step number (Matlab indexing).
%       OUTPUTS:
%           xDot - (n x tn number) State derivative.
%
%   g - (function_handle)
%       Policy.
%       SYNTAX:
%           u = g(x,t,k);
%       INPUTS:
%           x - (n x 1 number) State.
%           t - (1 x 1 number) Time.
%           k - (1 x 1 integer) Step number (Matlab indexing).
%       OUTPUTS:
%           u - (m x 1 number) Input.
%
%   t - (1 x tn number)
%       Time trajectory.
%
%   x0 - (n x 1 number) 
%       Initial state.
%
%   xm - (n x 1 number)
%       Minimum state constraint.
%
%   xM - (n x 1 number)
%       Maximum state constraint.
%   
%
% PROPERTIES: TODO: Add properties
%   'direction' - ('forward' or 'backward') ['forward']
%       Direction of simulation. If 'forward' then "x0" is the state at
%       time "t[0]", 'backward' then "x0" is the state at time "t[tn]".
% 
% OUTPUTS:
%   x - (n x tn number) 
%       Solution state trajectory.
%
%   u - (m x tn-1 number)
%       Input trajectory
%
%   xDot - (n x tn number) 
%       Solution state time derivative trajectory.
%
% EXAMPLES: TODO: Add examples
%
% NOTES:
%   To simulate backwards in time, the time trajectory must be given as
%   trajectory moving backwards in time.
%
% NECESSARY FILES:
%
% SEE ALSO:
%    optimal.simState | optimal.simCostate
%
% AUTHOR:
%    Rowland O'Flaherty (www.rowlandoflaherty.com)
%
% VERSION: 
%   Created 09-MAY-2014
%-------------------------------------------------------------------------------

%% Check Inputs

% Check number of inputs
narginchk(4,inf)

% Check input arguments for errors
assert(isa(f,'function_handle'),...
    'optimal:simulate:f',...
    'Input argument "f" must be a function handle.')

assert(isa(g,'function_handle'),...
    'optimal:simulate:g',...
    'Input argument "g" must be a function handle.')
m = size(g(0,0,1),1);

assert(isnumeric(t) && isreal(t) && isvector(t),...
    'optimal:simulate:t',...
    'Input argument "t" must be a vector of real numbers.')
t = t(:)';
tn = numel(t);

assert(isnumeric(x0) && isvector(x0),...
    'optimal:simulate:x0',...
    'Input argument "x0" must be a vector.')
x0 = x0(:);
n = size(x0,1);

if nargin < 5, xm = -inf*ones(n,1); end
assert(isnumeric(xm) && isreal(xm) && isvector(xm) && length(xm) == n,...
    'optimal:simulate:xm',...
    'Input argument "xm" must be a %d element vector of real numbers.',n)

if nargin < 6, xM = inf*ones(n,1); end
assert(isnumeric(xM) && isreal(xM) && isvector(xM) && length(xM) == n && all(xM >= xm),...
    'optimal:simulate:xM',...
    'Input argument "xM" must be a %d element vector of real numbers all greater or equal to xm.',n)

%% Get and check properties
propargin = size(varargin,2);

assert(mod(propargin,2) == 0,...
    'optimal:simulate:properties',...
    'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')

propStrs = varargin(1:2:propargin);
propValues = varargin(2:2:propargin);

for iParam = 1:propargin/2
    switch lower(propStrs{iParam})
        case lower('direction')
            direction = propValues{iParam};
        otherwise
            error('optimal:simulate:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('direction','var'), direction = 'forward'; end

% Check property values for errors TODO: Add property error checks
assert(ischar(direction) && ismember(direction,{'forward','backward'}),...
    'optimal:simulate:direction',...
    'Property "direction" must be either ''forward'' or ''backward''.')

%% Run
if strcmp(direction,'forward') % Simulate forward
    x = nan(n,tn);
    u = nan(m,tn);
    xDot = nan(n,tn-1);
    x(:,1) = x0;
    
    for k = 1:tn-1
        u(:,k) = g(x(:,k),t(k),k);
        xDot(:,k) = f(x(:,k),u(:,k),t(k),k);
        ts = t(k+1) - t(k);
        x(:,k+1) = constrain(x(:,k) + xDot(:,k)*ts,xm,xM);
    end
    u(:,tn) = g(x(:,tn),t(tn),tn);
    
else % Simulate backward
% If you simulate forward then backwards the answers will not be the same
% because the u is based on left-side state in forward simulation and
% right-side state in backward simulation. There might be way to make them
% equal by formately picking the u at each step the backward simulation as
% the one such that in one step forward simuation will result in the
% right-hand state. Probably the same with xDot.

    x = nan(n,tn);
    u = nan(m,tn);
    xDot = nan(n,tn-1);
    x(:,end) = x0;
    
    for k = tn:-1:2
        u(:,k) = g(x(:,k),t(k),k);
        xDot(:,k-1) = f(x(:,k),u(:,k-1),t(k),k);
        ts = t(k) - t(k-1);
        x(:,k-1) = constrain(x(:,k) - xDot(:,k-1)*ts,xm,xM);
    end
    u(:,1) = g(x(:,1),t(1),1);

end

end

%% Helper functions
function x = constrain(x,xm,xM)
x = min(max(x,xm),xM);
end
