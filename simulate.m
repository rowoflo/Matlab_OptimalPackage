function [x,xD] = simulate(f,x0,t,varargin)
% The "simulate" function integrates the state in time.
%
% SYNTAX:
%   [x,xD] = optimal.simulate(f,x0,t)
%   [x,xD] = optimal.simulate(f,x0,t,'PropertyName',PropertyValue,...)
% 
% INPUTS:
%   f - (function_handle)
%       State dynamics.
%       SYNTAX:
%           xDot = f(x,t);
%       INPUTS:
%           x - (n x tn number) State.
%           t - (1 x tn number) Time.
%       OUTPUTS:
%           xDot - (n x tn number) State derivative.
%
%   x0 - (n x 1 number) 
%       Initial conidtion.
%
%   t - (1 x tn number)
%       Time trajectory.
%
% PROPERTIES: TODO: Add properties
%   'propertiesName' - (size type) [defaultPropertyValue]
%       Description.
% 
% OUTPUTS:
%   x - (n x tn number) 
%       Solution trajectory.
%
%   xD - (n x tn number) 
%       Solution time derivative trajectory.
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
narginchk(3,inf)

% Check input arguments for errors
assert(isa(f,'function_handle'),...
    'optimal:simulate:f',...
    'Input argument "f" must be a function handle.')

assert(isnumeric(x0) && isvector(x0),...
    'optimal:simulate:x0',...
    'Input argument "x0" must be a vector.')
x0 = x0(:);
n = size(x0,1);

assert(isnumeric(t) && isreal(t) && isvector(t),...
    'optimal:simulate:t',...
    'Input argument "t" must be a vector of real numbers.')
t = t(:)';
tn = numel(t);

% % Get and check properties
% propargin = size(varargin,2);
% 
% assert(mod(propargin,2) == 0,...
%     'optimal:simulate:properties',...
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
%             error('optimal:simulate:options',...
%               'Option string ''%s'' is not recognized.',propStrs{iParam})
%     end
% end
% 
% % Set to default value if necessary TODO: Add property defaults
% if ~exist('propertyName','var'), propertyName = defaultPropertyValue; end
% 
% % Check property values for errors TODO: Add property error checks
% assert(isnumeric(propertyName) && isreal(propertyName) && isequal(size(propertyName),[1,1]),...
%     'optimal:simulate:propertyName',...
%     'Property "propertyName" must be a ? x ? matrix of real numbers.')

%% Initialize
x = nan(n,tn);
xD = nan(n,tn-1);
x(:,1) = x0;

%% Simulate forward
for k = 1:tn-1
    xD(:,k) = f(x(:,k),t(k));
    ts = t(k+1) - t(k);
    x(:,k+1) = x(:,k) + xD(:,k)*ts;
end

end
