function gamma = armijo(x,u,lambda,t,f,J,dHdu,alpha,beta,varargin)
% The "armijo" function calculates an appropriate step size for gradiant
% descent.
%
% SYNTAX:
%   gamma = optimal.armijo(x,u,lambda,t,f,J,dHdu)
%   gamma = optimal.armijo(x,u,lambda,t,f,J,dHdu,alpha,beta)
%   gamma = optimal.armijo(x,u,lambda,t,f,J,dHdu,alpha,beta,'PropertyName',PropertyValue,...)
% 
% INPUTS:
%   x - (n x tn number) 
%       State trajectory.
%
%   u - (m x tn-1 number)
%       Input trajectory.
%
%   lambda - (n x tn number)
%       Costate trajectory.
%
%   t - (1 x tn number)
%       Time trajectory.
%
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
%   J - (function handle)
%       Trajectory cost.
%       SYNTAX:
%           C  = J(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           C - (1 x 1 number) Trajectory cost.
%
%   dHdu - (function handle)
%       Hamiltonian partial to input.
%       SYNTAX:
%           hu  = dHdu(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           hu - (1 x n x tn number) Hamiltonian partial to input.
%
%   alpha - (1 x 1 0<=number<=1) [0.5]
%       Armijo parameter.
%
%   beta - (1 x 1 0<=number<=1) [0.5]
%       Armijo parameter.
%
% PROPERTIES:
%   'kappaMax' - (1 x 1 positive integer) [30]
%       Max number of iterations for the Armijo algorithm.
% 
% OUTPUTS:
%   gamma - (1 x 1 number) 
%       Calculated step size.
%
% EXAMPLES: TODO: Add examples
%
% NOTES:
%
% NECESSARY FILES:
%
% SEE ALSO: TODO: Add see alsos
%    relatedFunction1 | relatedFunction2
%
% AUTHOR:
%    Rowland O'Flaherty (www.rowlandoflaherty.com)
%
% VERSION: 
%   Created 02-MAY-2014
%-------------------------------------------------------------------------------

%% Check Inputs

% Check number of inputs
narginchk(7,inf)

% Apply default values
if nargin < 7, alpha = 0.5; end
if nargin < 8, beta = 0.5; end

% Check input arguments for errors
assert(isnumeric(t) && isreal(t) && isvector(t),...
    'optimal:armijo:t',...
    'Input argument "t" must be a vector of real numbers.')
t = t(:)';
tn = numel(t);

assert(isnumeric(x) && numel(size(u)) == 2 && size(x,2) == tn,...
    'optimal:armijo:x',...
    'Input argument "x" must be a matrix witha a length of %d.',tn)
n = size(x,1);

assert(isnumeric(u) && numel(size(u)) == 2 && size(u,2) == tn-1,...
    'optimal:armijo:u',...
    'Input argument "u" must be a matrix with a length of %d.',tn-1)

assert(isnumeric(lambda) && isequal(size(lambda),[n,tn]),...
    'optimal:armijo:lambda',...
    'Input argument "lambda" must be a %d x %d matrix of real numbers.',n,tn)

assert(isa(f,'function_handle'),...
    'optimal:armijo:f',...
    'Input argument "f" must be a function handle.')

assert(isa(J,'function_handle'),...
    'optimal:armijo:J',...
    'Input argument "J" must be a function handle.')

assert(isa(dHdu,'function_handle'),...
    'optimal:armijo:dHdu',...
    'Input argument "dHdu" must be a function handle.')

assert(isnumeric(alpha) && isreal(alpha) && numel(alpha) == 1 && alpha >= 0 && alpha <= 1,...
    'optimal:armijo:alpha',...
    'Input argument "alpha" must be a number between 0 and 1.')

assert(isnumeric(beta) && isreal(beta) && numel(beta) == 1 && beta >= 0 && beta <= 1,...
    'optimal:armijo:beta',...
    'Input argument "beta" must be a number between 0 and 1.')

% Get and check properties
propargin = size(varargin,2);

assert(mod(propargin,2) == 0,...
    'optimal:armijo:properties',...
    'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')

propStrs = varargin(1:2:propargin);
propValues = varargin(2:2:propargin);

for iParam = 1:propargin/2
    switch lower(propStrs{iParam})
        case lower('kappaMax')
            kappaMax = propValues{iParam};
        otherwise
            error('optimal:armijo:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('kappaMax','var'), kappaMax = 30; end

% Check property values for errors
assert(isnumeric(kappaMax) && isreal(kappaMax) && numel(kappaMax) == 1 && kappaMax > 0 && mod(kappaMax,1) == 0,...
    'optimal:armijo:kappaMax',...
    'Property "kappaMax" must be a positive integer.')

%% Initialize
dHduT = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
norm2dHduT = sum(sum(dHduT.*dHduT,1));
kappa = 0;
J0 = J(x,u,t);
J1 = inf;

%% Calculate step size
while J1 - J0 > -alpha*beta^kappa*norm2dHduT && kappa < kappaMax
    gamma = beta^kappa;
    u1 = u - gamma*dHduT;
    x1 = optimal.simState(f,x(:,1),u1,t);
    J1 = J(x1,u1,t);
    kappa = kappa + 1;
end


end
