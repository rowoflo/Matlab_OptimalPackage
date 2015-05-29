function gamma = armijo(t,x,u,lambda,H,Phi,dHdu,alpha,beta,gamma0,kappaMax)
% The "armijo" function calculates an appropriate step size for gradiant
% descent.
%
% SYNTAX:
%   gamma = optimal.armijo(t,x,u,lambda,H,Phi,dHdu,alpha,beta,gamma0,kappaMax)
% 
% INPUTS:
%   t - (1 x tn increasing)
%       Time trajectory.
%
%   x - (n x tn) 
%       State trajectory.
%
%   u - (m x tn)
%       Input trajectory.
%
%   lambda - (n x tn)
%       Costate trajectory.
%
%   H - (function_handle)
%       Trajectory cost.
%       SYNTAX:
%           C_H  = H(x,u,lambda,t);
%       INPUTS:
%           x - (n x tn number) State trajecory.
%           u - (m x tn number) Input trajecory.
%           lambda - (n x tn number) Costate trajectory.
%           t - (1 x tn number) Time trajecory.
%       OUTPUTS:
%           C - (1 x 1 number) Trajectory cost.
%
%   Phi - (function_handle)
%       Next state and input trajectories given a step size and current
%       trajectories and step size.
%       SYNTAX:
%           [x1,u1] = Phi(x,u,lambda,gamma);
%       INPUTS:
%           x - (n x tn number) State trajecory.
%           u - (m x tn number) Input trajecory.
%           lambda - (n x tn number) Costate trajectory.
%           gamma - (1 x 1 number) Step size.
%       OUTPUTS:
%           x1 - (n x tn number) Next state trajectory.
%           u1 - (m x tn number) Next input trajectory.
%           lambda1 - (n x tn number) Next costate trajectory.
%
%   dHdu - (function_handle)
%       Cost partial to input.
%       SYNTAX:
%           C_dHdu  = dHdu(x,u,lambda,t);
%       INPUTS:
%           x - (n x tn number) State trajecory.
%           u - (m x tn number) Input trajecory.
%           lambda - (n x tn number) Costate trajectory.
%           t - (1 x tn number) Time trajecory.
%       OUTPUTS:
%           C_dHdu - (m x tn number) Partial cost to input.
%
%   alpha - (1 x 1 0<=number<=1) [0.5]
%       Armijo parameter.
%
%   beta - (1 x 1 0<=number<=1) [0.5]
%       Armijo parameter.
%
%   gamma0 - (1 x 1 0<=number) [1]
%       Initial step size.
%
%   kappaMax - (1 x 1 positive integer) [10]
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

% Check Inputs

% Check number of inputs
narginchk(7,inf)

% Apply default values
if nargin < 8, alpha = 0.5; end
if nargin < 9, beta = 0.5; end
if nargin < 10, gamma0 = 1; end
if nargin < 11, kappaMax = 10; end

% Check input arguments for errors
assert(isnumeric(t) && isreal(t) && isvector(t),...
    'optimal:armijo:t',...
    'Input argument "t" must be a vector of real numbers.')
t = t(:)';
tn = numel(t);

assert(isnumeric(x) && numel(size(x)) == 2 && size(x,2) == tn,...
    'optimal:armijo:x',...
    'Input argument "x" must be a matrix with a length of %d.',tn)

assert(isnumeric(u) && numel(size(u)) == 2 && size(u,2) == tn,...
    'optimal:armijo:u',...
    'Input argument "u" must be a matrix with a length of %d.',tn)

assert(isnumeric(lambda) && numel(size(lambda)) == 2 && size(lambda,2) == tn,...
    'optimal:armijo:lambda',...
    'Input argument "lambda" must be a matrix with a length of %d.',tn)

assert(isa(H,'function_handle'),...
    'optimal:armijo:H',...
    'Input argument "H" must be a function handle.')

assert(isa(Phi,'function_handle'),...
    'optimal:armijo:Phi',...
    'Input argument "Phi" must be a function handle.')

assert(isa(dHdu,'function_handle'),...
    'optimal:armijo:dHdu',...
    'Input argument "dHdu" must be a function handle.')

assert(isnumeric(alpha) && isreal(alpha) && numel(alpha) == 1 && alpha >= 0 && alpha <= 1,...
    'optimal:armijo:alpha',...
    'Input argument "alpha" must be a number between 0 and 1.')

assert(isnumeric(beta) && isreal(beta) && numel(beta) == 1 && beta >= 0 && beta < 1,...
    'optimal:armijo:beta',...
    'Input argument "beta" must be a number between 0 and 1.')

assert(isnumeric(gamma0) && isreal(gamma0) && numel(gamma0) == 1 && gamma0 > 0,...
    'optimal:armijo:gamma0',...
    'Input argument "gamma0" must be a number positive number.')

assert(isnumeric(kappaMax) && isreal(kappaMax) && numel(kappaMax) == 1 && kappaMax > 0 && mod(kappaMax,1) == 0,...
    'optimal:armijo:kappaMax',...
    'Property "kappaMax" must be a positive integer.')

%% Initialize
tDel = mean(diff(t));
dHdu_norm = norm(dHdu(x,u,lambda,t))*tDel;
kappa = 0;
H0 = H(x,u,lambda,t);
H1 = inf;
gamma = gamma0;

%% Run
while H1 - H0 > -alpha*dHdu_norm*beta^kappa && kappa < kappaMax
    gamma = gamma0*beta^kappa;
    [x1, u1, lambda1] = Phi(x,u,lambda,gamma);
    H1 = H(x1,u1,lambda1,t);
    kappa = kappa + 1;
end

end
