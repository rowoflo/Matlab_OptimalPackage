function [x,u,xDot] = simulate(t,x0,f,g,varargin)
% The "simulate" function integrates the state in time.
%
% SYNTAX:
%   [x,u,xDot] = optimal.simulate(t,x0,f,g)
%   [x,u,xDot] = optimal.simulate(t,x0,f,g,'PropertyName',PropertyValue,...)
% 
% INPUTS:
%   t - (1 x tn increasing)
%       Time trajectory.
%
%   x0 - (n x 1) 
%       Initial state.
%
%   f - (function_handle)
%       State dynamics.
%       SYNTAX:
%           xDot = f(x,u,t,k);
%       INPUTS:
%           x - (n x 1) State.
%           u - (m x 1) Input.
%           t - (1 x 1) Time.
%           k - (1 x 1 integer) Step number (Matlab indexing).
%       OUTPUTS:
%           xDot - (n x tn number) State derivative.
%
%   g - (function_handle)
%       Policy.
%       SYNTAX:
%           u = g(x,t,k);
%       INPUTS:
%           x - (n x 1) State.
%           t - (1 x 1) Time.
%           k - (1 x 1 integer) Step number (Matlab indexing).
%       OUTPUTS:
%           u - (m x 1) Input.
%
%   xm - (n x 1)
%       Minimum state constraint.
%
%   xM - (n x 1)
%       Maximum state constraint.
%   
%
% PROPERTIES:
%   'statemin' - (n x 1) [-inf*ones(n,1)]
%       Minimum state constraint.
%
%   'statemax' - (n x 1) [inf*ones(n,1)]
%       Maximum state constraint.
%
%   'inputmin' - (m x 1) [-inf*ones(m,1)]
%       Minimum input constraint.
%
%   'inputmax' - (m x 1) [inf*ones(m,1)]
%       Maximum input constraint.
%
%   'direction' - ('forward' or 'backward') ['forward']
%       Direction of simulation. If 'forward' then "x0" is the state at
%       time "t[0]", 'backward' then "x0" is the state at time "t[tn]".
% 
% OUTPUTS:
%   x - (n x tn) 
%       Solution state trajectory.
%
%   u - (m x tn-1)
%       Input trajectory
%
%   xDot - (n x tn) 
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

assert(isa(f,'function_handle'),...
    'optimal:simulate:f',...
    'Input argument "f" must be a function handle.')

assert(isa(g,'function_handle'),...
    'optimal:simulate:g',...
    'Input argument "g" must be a function handle.')
m = size(g(0,0,1),1);

%% Get and check properties
propargin = size(varargin,2);

assert(mod(propargin,2) == 0,...
    'optimal:simulate:properties',...
    'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')

propStrs = varargin(1:2:propargin);
propValues = varargin(2:2:propargin);

for iParam = 1:propargin/2
    switch lower(propStrs{iParam})
        case lower('statemin')
            xm = propValues{iParam};
        case lower('statemax')
            xM = propValues{iParam};
        case lower('inputmin')
            um = propValues{iParam};
        case lower('inputmax')
            uM = propValues{iParam};
        case lower('direction')
            direction = propValues{iParam};
        otherwise
            error('optimal:simulate:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('xm','var'), xm = -inf*ones(n,1); end
if ~exist('xM','var'), xM = inf*ones(n,1); end
if ~exist('um','var'), um = -inf*ones(m,1); end
if ~exist('uM','var'), uM = inf*ones(m,1); end
if ~exist('direction','var'), direction = 'forward'; end

% Check property values for errors
assert(isnumeric(xm) && isreal(xm) && isvector(xm) && length(xm) == n,...
    'optimal:simulate:xm',...
    'Input argument "xm" must be a %d element vector of real numbers.',n)

assert(isnumeric(xM) && isreal(xM) && isvector(xM) && length(xM) == n && all(xM >= xm),...
    'optimal:simulate:xM',...
    'Input argument "xM" must be a %d element vector of real numbers all greater or equal to minimum state constraint.',n)

assert(isnumeric(um) && isreal(um) && isvector(um) && length(um) == m,...
    'optimal:simulate:um',...
    'Input argument "um" must be a %d element vector of real numbers.',m)

assert(isnumeric(uM) && isreal(uM) && isvector(uM) && length(uM) == m && all(uM >= um),...
    'optimal:simulate:uM',...
    'Input argument "uM" must be a %d element vector of real numbers all greater or equal to minimum input constraint.',m)

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
        u(:,k) = constrain(g(x(:,k),t(k),k),um,uM);
        xDot(:,k) = f(x(:,k),u(:,k),t(k),k);
        ts = t(k+1) - t(k);
        x(:,k+1) = constrain(x(:,k) + xDot(:,k)*ts,xm,xM);
    end
    u(:,tn) = constrain(g(x(:,tn),t(tn),tn),um,uM);
    
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
        u(:,k) = constrain(g(x(:,k),t(k),k),um,uM);
        xDot(:,k-1) = f(x(:,k),u(:,k),t(k),k);
        ts = t(k) - t(k-1);
        x(:,k-1) = constrain(x(:,k) - xDot(:,k-1)*ts,xm,xM);
    end
    u(:,1) = constrain(g(x(:,1),t(1),1),um,uM);

end

end

%% Helper functions
function x = constrain(x,xm,xM)
x = min(max(x,xm),xM);
end
