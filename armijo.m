function gamma = armijo(x,u,lambda,t,f,J,dHdu,alpha,beta,varargin)
% The "armijo" function ...  TODO: Add description
%
% SYNTAX: TODO: Add syntax
%   output = armijo(input1)
%   output = armijo(input1,input2)
%   output = armijo(input1,input2,'PropertyName',PropertyValue,...)
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
% 
% % Check number of inputs TODO: Add number argument check
% narginchk(1,inf)
% 
% % Apply default values TODO: Add apply defaults
% if nargin < 2, input2 = defaultInputValue; end
% 
% % Check input arguments for errors TODO: Add error checks
% assert(isnumeric(input1) && isreal(input1) && isequal(size(input1),[1,1]),...
%     'optimal:armijo:input1',...
%     'Input argument "input1" must be a ? x ? matrix of real numbers.')
% 
% % Get and check properties
% propargin = size(varargin,2);
% 
% assert(mod(propargin,2) == 0,...
%     'optimal:armijo:properties',...
%     'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')
% 
% propStrs = varargin(1:2:propargin);
% propValues = varargin(2:2:propargin);
% 
% for iParam = 1:propargin/2
%     switch lower(propStrs{iParam})
%         case lower('propertyName')
%             propertyName = propValues{iParam};
%         otherwise
%             error('optimal:armijo:options',...
%               'Option string ''%s'' is not recognized.',propStrs{iParam})
%     end
% end
% 
% % Set to default value if necessary TODO: Add property defaults
% if ~exist('propertyName','var'), propertyName = defaultPropertyValue; end
% 
% % Check property values for errors TODO: Add property error checks
% assert(isnumeric(propertyName) && isreal(propertyName) && isequal(size(propertyName),[1,1]),...
%     'optimal:armijo:propertyName',...
%     'Property "propertyName" must be a ? x ? matrix of real numbers.')

%% Initialize
dHduTraj = permute(dHdu(x(:,1:end-1),u,lambda(:,1:end-1),t(1:end-1)),[2 3 1]);
norm2dHdu = sum(sum(dHduTraj.*dHduTraj,1));
kappa = 0;
J0 = J(x,u,t);
J1 = inf;
kappaMax = 30;

%% Calculate step size
while J1 - J0 > -alpha*beta^kappa*norm2dHdu && kappa < kappaMax
    gamma = beta^kappa;
    u1 = u - gamma*dHduTraj;
    x1 = optimal.simStateForward(f,x(:,1),u1,t);
    J1 = J(x1,u1,t);
    kappa = kappa + 1;
end


end
