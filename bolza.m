function [x,u,lambda,J,dHdu,JTape,gammaTape] = bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,varargin)
% The "bolza" function solves the Bolza optimal control problem. The Bolza
% problem is defined as
%
% $$\min_u \int_{t_0}^{t_f} L(x,u,t) dt + \Psi(x_f,t_f)$$
% $$s.t.~ \dot{x} = f(x,u,t),~x(0) = x_0$$
%
% SYNTAX:
%   [x,u,lambda,J,dHdu,JTape,gammaTape] = optimal.bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx)
%   [x,u,lambda,J,dHdu,JTape,gammaTape] = optimal.bolza(t,x0,u,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,'PropertyName',PropertyValue,...)
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
%           xf - (n x 1 number) Final state.
%           tf - (1 x 1 number) Final time.
%       OUTPUTS:
%           cf - (1 x 1 number) Final cost.
%
%   dPsidx - (function handle)
%       Final cost partial to final state
%       SYNTAX:
%           cfx  = dPsidx(xf,tf);
%       INPUTS:
%           xf - (n x 1 number) Final state.
%           tf - (1 x 1 number) Final time.
%       OUTPUTS:
%           cfx - (1 x n number) Final cost partial to final state.
%           
%
% PROPERTIES:
%   'armijoAlpha' - (1 x 1 0<=number<=1) [0.5]
%       Armijo alpha parameter. Aggresion on change in cost.
%
%   'armijoBeta' - (1 x 1 0<number<1) [0.5]
%       Armijo beta parameter. Step size controller. If output blows up try
%       lowering this value.
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
%   'method' - ('sweep', 'shooting', 'fmin') ['sweep']
%       Algorithm solution method.
%
%   'display' - ('none', 'iter') ['none']
%       Level of display.
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
%   dHdu - (function handle)
%       Hamiltonian partial to input.
%       SYNTAX:
%           hu  = dHdu(x,u,t);
%       INPUTS:
%           x - (n x tn number) State.
%           u - (m x tn number) Input.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           hu - (1 x m x tn number) Hamiltonian partial to input.
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
        case lower('armijoAlpha')
            alpha = propValues{iParam};
        case lower('armijoBeta')
            beta = propValues{iParam};
        case lower('stoppingCondition')
            stoppingCondition = propValues{iParam};
        case lower('method')
            method = propValues{iParam};
        case lower('display')
            displayLevel = propValues{iParam};
        otherwise
            error('optimal:bolza:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('alpha','var'), alpha = 0.5; end
if ~exist('beta','var'), beta = 0.5; end
if ~exist('stoppingCondition','var'), stoppingCondition = @stopDefault; end
if ~exist('method','var'), method = 'sweep'; end
if ~exist('displayLevel','var') displayLevel = 'none'; end

% Check property values for errors
assert(isnumeric(alpha) && isreal(alpha) && numel(alpha) && alpha >= 0 && alpha <= 1,...
    'optimal:bolza:armijoAlpha',...
    'Property "armijoAlpha" must be a number between 0 and 1.')

assert(isnumeric(beta) && isreal(beta) && numel(beta) && beta > 0 && beta < 1,...
    'optimal:bolza:armijoBeta',...
    'Property "armijoBeta" must be a number between 0 and 1.')

assert(isa(stoppingCondition,'function_handle'),...
    'optimal:bolza:stoppingCondition',...
    'Property "stoppingCondition" must be a function handle.')

assert(ismember(method,{'sweep','shooting','fmin'}),...
    'optimal:bolza;method',...
    'Property "method" must be either ''sweep'', ''shooting'', ''fmin''.')

assert(ismember(displayLevel,{'none','iter'}),...
    'optimal:bolza;displayLevel',...
    'Property "displayLevel" must be either ''none'', ''iter''.')

%% Initialize
% Hamiltonian
H = @(x_,u_,lambda_,t_) L(x_,u_,t_) + sum(lambda_.*f(x_,u_,t_),1); % (1 x 1) Hamiltonian
dHdx = @(x_,u_,lambda_,t_) dLdx(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 n 1]).*dfdx(x_,u_,t_),1); % (1 x n) Hamiltonian partial to state

% State
x = optimal.simState(f,x0,u,t); % (n x tn) State trajectory
xf = x(:,end); % (n x 1) Final state

% Costate
lambdaf = @(xf_) dPsidx(xf_)'; % (n x 1) Costate at final time
g = @(x_,u_,lambda_,t_) -dHdx(x_,u_,lambda_,t_)'; % (n x 1) Costate dynamics (i.e lambdaDot)
lambda = optimal.simCostate(g,lambdaf(xf),x,u,t,false); % (n x tn) Costate trjectory)

% Hamiltonian
dHdu = @(x_,u_,lambda_,t_) dLdu(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 m 1]).*dfdu(x_,u_,t_),1); % (1 x m) Hamiltonian partial to input
dHduT = permute(dHdu(x(:,1:end-1),u(:,1:end-1),lambda(:,1:end-1),t(1:end-1)),[2 3 1]);

% Cost
J = @(x_,u_,t_) sum(L(x_(:,1:end-1),u_(:,1:end-1),t_(1:end-1)).*diff(t)) + Psi(x_(:,end),t_(end)); % (1 x 1) Cost

% Trajectory update
Phi = @(x_,u_,lambda_,gamma_) trajUpdate(x_,u_,lambda_,gamma_,dHdu);

% Iteration count
k = 0;

% Records
JTape = nan;
gammaTape = nan;

%% Solve
dHduT = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
dHduTnorm2Init = sum(dHduT(:).^2);
JInit = J(x,u,t);
switch method
    case 'sweep'
        while ~stoppingCondition(x,u,lambda,t,k,dHduT)
            % Increment counter
            k = k + 1;
            
            % Simulate state forward
            x = optimal.simState(f,x0,u,t);
            xf = x(:,end);
            
            % Simulate costate backward
            lambda = optimal.simCostate(g,lambdaf(xf),x,u,t,false);
            
            % Calculate step size
            gamma = optimal.armijo(t,x,u,lambda,H,Phi,dHdu,alpha,beta);
            
            % Update records
            JTape(k) = J(x,u,t);
            gammaTape(k) = mean(gamma(:));
            
            % Update input
            dHduT = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
            dHduTnomr2 = sum(dHduT(:).^2);
            
            if JTape(k) == inf || isnan(dHduTnomr2) || (JTape(k) > JInit && k >= 10)
                fprintf('Unable to converge. %.3f\n\n',dHduTnomr2)
                break
            end
          
            u = u - repmat(gamma,[m,1]).*dHduT;
            
            
            switch displayLevel
                case 'iter'
                    fprintf('%d\tCost: %.3f\tdHdu Norm: %.3f\n',k,JTape(k),sum(dHduT(:).^2))
            end
        end
        
    case 'shooting'
        gamma = .1;
        xBar = fminsearch(@(x_) Psi(x_,t(end)),x(end));
        lambdaEps = 0.1;
        lambda0i = zeros(n,n+1);
        lambda0i(:,2:end) = repmat(lambda0i(:,1),1,n) + diag(lambdaEps*ones(n,1));
        lambdaD = nan(n,tn);
        xD = nan(n,tn-1);
        ui = nan(n,tn-1);
        G = @(lambdaf_) (lambdaf_ - lambdaf(xBar))'*(lambdaf_ - lambdaf(xBar));
        Gi = nan(1,n+1);
        while ~stoppingCondition(x,u,lambda,t,k,dHduT)
            % Increment counter
            k = k + 1;
            
            % Simulate forward
            for i = (n+1):-1:1
                lambda(:,1) = lambda0i(:,i);
                for j = 1:tn-1
                    ts = t(j+1) - t(j);
                    ui(:,j) = fminsearch(@(u_) dHdu(x(:,j),u_,lambda(:,j),t(j))'*dHdu(x(:,j),u_,lambda(:,j),t(j)), u(:,j));
                    xD(:,j) = f(x(:,j),ui(:,j),t(j));
                    x(:,j+1) = x(:,j) + xD(:,j)*ts;
                    lambdaD(:,j) = g(x(:,j),ui(:,j),lambda(:,j),t(j));
                    lambda(:,j+1) = lambda(:,j) + lambdaD(:,j)*ts;
                end
                Gi(i) = G(lambda(:,end));
            end
            u = ui;
            dGdlambda0 = (Gi(2:end) - Gi(1)) / lambdaEps;
            lambda0i(:,1) = lambda0i(:,1) - gamma*dGdlambda0';
            lambda0i(:,2:end) = repmat(lambda0i(:,1),1,n) + diag(lambdaEps*ones(n,1));
            
        end
        
    case 'fmin'
        options = optimset('Display','iter');
        fun = @(u_) norm2dHdu(u_,m,tn,x0,t,f,g,lambdaf,dHdu);
        u = lsqnonlin(fun,reshape(u,[1 m*(tn-1)]),-1*ones(1,m*(tn-1)),1*ones(1,m*(tn-1)),options);        
        u = reshape(u,[m,tn-1]);
end
x = optimal.simState(f,x0,u,t);

end

function stopFlag = stopDefault(~,~,~,~,k,~)
% norm2dHduT = sum(sum(dHduT.*dHduT,1));
% stopFlag = norm2dHduT < size(dHduT,2)/10 | k >= 10;
stopFlag = k >= 10;
end

function value = norm2dHdu(u,m,tn,x0,t,f,g,lambdaf,dHdu)
u = reshape(u,[m,tn-1]);
x = optimal.simState(f,x0,u,t);
xf = x(:,end);
lambda = optimal.simCostate(g,lambdaf(xf),x,u,t,false);
value = repmat(sum(permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]).^2,1),1,m);
end

function [x1,u1,lambda1] = trajUpdate(x,u,lambda,gamma,t,dHdu)
dHduT = permute(dHdu(x(:,1:end-1),u(:,1:end-1),lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
end


