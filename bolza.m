function [x,u,lambda,JTape,gammaTape] = bolza(t,x0,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,varargin)
% The "bolza" function ...  TODO: Add description
%
% SYNTAX: TODO: Add syntax
%   output = bolza(input1)
%   output = bolza(input1,input2)
%   output = bolza(input1,input2,'PropertyName',PropertyValue,...)
% 
% INPUTS: TODO: Add inputs
%   input1 - (size type) 
%       Description.
%
%   input2 - (size type) [defaultInputValue] 
%       Description for optional input.
%
% PROPERTIES: TODO: Add properties
%   'propertiesName' - (size type) [defaultPropertyValue]
%       Description.
% 
% OUTPUTS: TODO: Add outputs
%   output - (size type) 
%       Description.
%
% EXAMPLES: TODO: Add examples
%
% NOTES:
%
% NECESSARY FILES: TODO: Add necessary files
%   +somePackage, someFile.m
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

% Check input arguments for errors TODO: Add error checks
% assert(isnumeric(input1) && isreal(input1) && isequal(size(input1),[1,1]),...
%     'optimal:bolza:input1',...
%     'Input argument "input1" must be a ? x ? matrix of real numbers.')

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
        otherwise
            error('optimal:bolza:options',...
              'Option string ''%s'' is not recognized.',propStrs{iParam})
    end
end

% Set to default value if necessary
if ~exist('armijoParams','var'), armijoParams = [0.5 0.5]; end

% Check property values for errors
assert(isnumeric(armijoParams) && isreal(armijoParams) && length(armijoParams) == 2,...
    'optimal:bolza:armijoParams',...
    'Property "armijoParams" must be a 2 element vector.')
alpha = armijoParams(1);
beta = armijoParams(2);

%% Initialize
% Time - variables
tn = length(t); % (1 x 1) Number of time samples

% State
n = size(x0,1); % (1 x 1) Dimension of the state
x = nan(n,tn); % (n x tn) State vector record

% Input
m = size(dfdu(0,0,0),2);
% u = zeros(m,tn-1); % (m x tn) Input vector record
% u = sin(2*pi/t(end)*t(1:end-1));
u = rand(m,tn-1);

% Hamiltonian
H = @(x_,u_,lambda_,t_) L(x_,u_,t_) + sum(lambda_.*f(x_,u_,t_),1); % (1 x 1) Hamiltonian
dHdx = @(x_,u_,lambda_,t_) dLdx(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 n 1]).*dfdx(x_,u_,t_),1); % (1 x n) Hamiltonian partial to state
dHdu = @(x_,u_,lambda_,t_) dLdu(x_,u_,t_) + sum(repmat(permute(lambda_,[1,3,2]),[1 m 1]).*dfdu(x_,u_,t_),1); % (1 x m) Hamiltonian partial to input

% Costate
lambda = nan(n,tn); % (n x tn) Costate vector record
lambdaf = @(xf_) dPsidx(xf_)'; % (n x 1) Costate at final time
g = @(x_,u_,lambda_,t_) -dHdx(x_,u_,lambda_,t_)'; % (n x 1) Costate dynamics (i.e lambdaDot)

% Cost
J = @(x_,u_,t_) sum(L(x_(:,1:end-1),u_(1:end-1),t_(1:end-1))) + Psi(x_(end),t_(end)); % (1 x 1) Cost

% Records
JTape = nan;
gammaTape = nan;

%% Solve
x = optimal.simStateForward(f,x0,u,t);
k = 0;
C = 10;
while ~stop(x,u,t,k,C)
    % Increment counter
    k = k + 1;
    
    % Simulate state forward
    x = optimal.simStateForward(f,x0,u,t);
    xf = x(end);
    
    % Simulate costate backward
    lambda = optimal.simCostateBackward(g,lambdaf(xf),x,u,t);
    
    % Calculate step size
    gamma = optimal.armijo(x,u,lambda,t,f,J,dHdu,alpha,beta);
    
    % Update records
    JTape(k) = J(x,u,t);
    gammaTape(k) = gamma;
    
    % Update input
    dHduVec = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[1 3 2]);
    u = u - gamma*dHduVec;
    
end
x = optimal.simStateForward(f,x0,u,t);


end

function stopFlag = stop(x,u,t,k,C)
stopFlag = k > C;
end


