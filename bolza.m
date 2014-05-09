function [x,u,lambda,J,JTape,gammaTape] = bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,varargin)
% The "bolza" function solves the Bolza optimal control problem. The Bolza
% problem is defined as
%
% $$\min_u \int_{t_0}^{t_f} L(x,u,t) dt + \Psi(x_f,t_f)$$
% $$s.t.~ \dot{x} = f(x,u,t),~x(0) = x_0$$
%
% SYNTAX:
%   [x,u,lambda,JTape,gammaTape] = optimal.bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx)
%   [x,u,lambda,JTape,gammaTape] = optimal.bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,'PropertyName',PropertyValue,...)
%
% NOTATION:
%   n - State dimension.
%   m - Input dimension.
%   tn - Number of time samples.
%   kn - Number of iterations.
% 
% INPUTS:
%   t - (1 x tn number) 
%       Row vector of all time points.
%
%   x0 - (n x 1 number)
%       Initial state.
%
%   u - (m x tn-1 number)
%       Initial input trajectory.
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
%   dfdx - (function handle)
%       State dynamics partial derivative to state.
%       SYNTAX:
%           A = dfdx(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           A - (n x n x tn number) Partial to state.
%
%   dfdu - (function handle)
%       State dynamics partial derivative to input.
%       SYNTAX:
%           B = dfdx(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           B - (n x m x tn number) Partial to input.
%
%   L - (function handle)
%       Instantaneous cost.
%       SYNTAX:
%           c  = L(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           c - (1 x tn number) Cost.
%
%   dLdx - (function handle)
%       Instantaneous cost partial to state.
%       SYNTAX:
%           cx  = dLdx(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           cx - (1 x n x tn number) Cost partial to state.
%
%   dLdu - (function handle)
%       Instantaneous cost partial to input.
%       SYNTAX:
%           cu  = dLdu(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           cu - (1 x m x tn number) Cost partial to input.
%
%   Psi - (function handle)
%       Final cost.
%       SYNTAX:
%           cf  = Psi(xf,tf);
%       INPUTS:
%           xf - (n x tn number) Final state.
%           tf - (1 x tn number) Final time.
%       OUTPUTS:
%           cf - (1 x 1 number) Final cost.
%
%   dPsidx - (function handle)
%       Final cost partial to final state
%       SYNTAX:
%           cfx  = dPsidx(xf,tf);
%       INPUTS:
%           xf - (n x tn number) Final state.
%           tf - (1 x tn number) Final time.
%       OUTPUTS:
%           cfx - (1 x n number) Final cost partial to final state.
%           
%
% PROPERTIES:
%   'armijoParams' - (1 x 2 0<=number<=1) [ [0.5 0.5] ]
%       Armijo parameters [alpha beta].
%
%   'stoppingCondition' - (function handle)
%       Iteration stopping condition
%       SYNTAX:
%           stopFlag  = stop(x,u,t,k,dHduT);
%       INPUTS:
%           x - (n x tn number) Current state trajectory.
%           u - (m x tn-1 number) Current input trajectory.
%           t - (1 x tn number) Current time trajectory.
%           k - (1 x 1 number) Current iteration number.
%           dHduT - (m x tn-1 number) Current dHdu tranpose trajectory.
%       OUTPUTS:
%           stopFlag - (1 x 1 logical) If true the optimization iteration
%               will stop.
% 
% OUTPUTS:
%   x - (n x tn number) 
%       Locally optimal state trajectory.
%   
%   u - (m x tn-1 number)
%       Locally optimal input trajectory.
%
%   lambda - (n x tn number)
%       Locally optimal costate trajectory.
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
%   JTape - (1 x kn number)
%       History of cost verse iteration number.
%
%   gammaTape - (1 x kn number)
%       History of step size verse iteration number.
%
% EXAMPLES:
%   See "bolzaExamples.m" script for examples.
%
% NOTES:
%
% NECESSARY FILES:
%   +optimal, armjo.m, simState.m, simCostate.m
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
narginchk(10,inf)

% Apply default values TODO: Add apply defaults
% if nargin < 2, input2 = defaultInputValue; end

% Check input arguments for errors
assert(isnumeric(t) && isvector(t),...
    'optimal:bolza:t',...
    'Input argument "t" must be a vector of numbers.')
t = t(:)';
tn = length(t); % Number of time samples

assert(isnumeric(x0) && isvector(x0),...
    'optimal:bolza:x0',...
    'Input argument "x0" must be a vector of numbers.')
x0 = x0(:);
n = size(x0,1); % Stete dimension

assert(isnumeric(u) && size(u,2) == tn-1,...
    'optimal:bolza:u',...
    'Input argument "u" must be a matrix of size m x %d.',tn-1)
m = size(u,1); % Input dimension

assert(isa(f,'function_handle'),...
    'optimal:bolza:f',...
    'Input argument "f" must be a function handle.')

assert(isa(dfdx,'function_handle'),...
    'optimal:bolza:dfdx',...
    'Input argument "dfdx" must be a function handle.')

assert(isa(dfdu,'function_handle'),...
    'optimal:bolza:dfdu',...
    'Input argument "dfdu" must be a function handle.')

assert(isa(L,'function_handle'),...
    'optimal:bolza:L',...
    'Input argument "L" must be a function handle.')

assert(isa(dLdx,'function_handle'),...
    'optimal:bolza:dLdx',...
    'Input argument "dLdx" must be a function handle.')

assert(isa(dLdu,'function_handle'),...
    'optimal:bolza:dLdu',...
    'Input argument "dLdu" must be a function handle.')

assert(isa(Psi,'function_handle'),...
    'optimal:bolza:Psi',...
    'Input argument "Psi" must be a function handle.')

assert(isa(dPsidx,'function_handle'),...
    'optimal:bolza:dPsidx',...
    'Input argument "dPsidx" must be a function handle.')

% Get and check properties
propargin = size(varargin,2);

assert(mod(propargin,2) == 0,...
    'optimal:bolza:properties',...
    'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')

propStrs = varargin(1:2:propargin);
propValues = varargin(2:2:propargin);

for iParam = 1:propargin/2
    switch lower(propStrs{iParam})
        case lower('armijoParams')
            armijoParams = propValues{iParam};
        case lower('stoppingCondition')
            stoppingCondition = propValues{iParam};
        otherwise
            error('optimal:bolza:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('armijoParams','var'), armijoParams = [0.5 0.5]; end
if ~exist('stoppingCondition','var'), stoppingCondition = @stopDefault; end

% Check property values for errors
assert(isnumeric(armijoParams) && isreal(armijoParams) && length(armijoParams) == 2,...
    'optimal:bolza:armijoParams',...
    'Property "armijoParams" must be a 2 element vector.')
alpha = armijoParams(1);
beta = armijoParams(2);

assert(isa(stoppingCondition,'function_handle'),...
    'optimal:bolza:stoppingCondition',...
    'Property "stoppingCondition" must be a function handle.')

%% Initialize
% Hamiltonian
H = @(x_,u_,lambda_,t_) L(x_,u_,t_) + sum(lambda_.*f(x_,u_,t_),1); % (1 x 1) Hamiltonian
dHdx = @(x_,u_,lambda_,t_) dLdx(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 n 1]).*dfdx(x_,u_,t_),1); % (1 x n) Hamiltonian partial to state
dHdu = @(x_,u_,lambda_,t_) dLdu(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 m 1]).*dfdu(x_,u_,t_),1); % (1 x m) Hamiltonian partial to input

% Costate
lambda = nan(n,tn); % (n x tn) Costate vector record
lambdaf = @(xf_) dPsidx(xf_)'; % (n x 1) Costate at final time
g = @(x_,u_,lambda_,t_) -dHdx(x_,u_,lambda_,t_)'; % (n x 1) Costate dynamics (i.e lambdaDot)

% Cost
J = @(x_,u_,t_) sum(L(x_(:,1:end-1),u_,t_(1:end-1)).*diff(t)) + Psi(x_(:,end),t_(end)); % (1 x 1) Cost

% Records
JTape = nan;
gammaTape = nan;

%% Solve
x = optimal.simState(f,x0,u,t);
xf = x(:,end);
lambda = optimal.simCostate(g,lambdaf(xf),x,u,t,false);
dHduT = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
k = 0;
while ~stoppingCondition(x,u,lambda,t,k,dHduT)
    % Increment counter
    k = k + 1;
    
    % Simulate state forward
    x = optimal.simState(f,x0,u,t);
    xf = x(:,end);
    
    % Simulate costate backward
    lambda = optimal.simCostate(g,lambdaf(xf),x,u,t,false);
    
    % Calculate step size
    gamma = optimal.armijo(x,u,lambda,t,f,g,lambdaf,H,dHdu,alpha,beta);
%     gamma = .05;

    % Update records
    JTape(k) = J(x,u,t);
    gammaTape(k) = gamma;
    
    % Update input
    dHduT = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
    u = u - gamma.*dHduT;
    
end
x = optimal.simState(f,x0,u,t);


end

function stopFlag = stopDefault(~,~,~,~,k,dHduT)
norm2dHduT = sum(sum(dHduT.*dHduT,1))
stopFlag = norm2dHduT < 100;
% stopFlag = k > 50;
end


