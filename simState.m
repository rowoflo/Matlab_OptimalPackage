function [x,xD] = simState(f,x0,u,t,varargin)
% The "simState" function integrates the state in time.
%
% SYNTAX:
%   [x,xD] = optimal.simState(f,x0,u,t)
%   [x,xD] = optimal.simState(f,x0,u,t,'PropertyName',PropertyValue,...)
% 
% INPUTS:
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
%   x0 - (n x 1 number) 
%       Initial state.
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
%   x - (n x tn number) 
%       State trajectory.
%
%   xD - (n x tn number) 
%       State time derivative trajectory.
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

% % Check number of inputs
% narginchk(4,inf)
% 
% % Check input arguments for errors
% assert(isa(f,'function_handle'),...
%     'optimal:simState:f',...
%     'Input argument "f" must be a function handle.')
% 
% assert(isnumeric(t) && isreal(t) && isvector(t),...
%     'optimal:simState:t',...
%     'Input argument "t" must be a vector of real numbers.')
% t = t(:)';
tn = numel(t);
% 
% assert(isnumeric(u) && numel(size(u)) == 2 && size(u,2) == tn-1,...
%     'optimal:simState:u',...
%     'Input argument "u" must be a matrix with a length of %d.',tn-1)
% 
% % Get and check properties
% propargin = size(varargin,2);
% 
% assert(mod(propargin,2) == 0,...
%     'optimal:simState:properties',...
%     'Properties must come in pairs of a "PropertyName" and a "PropertyValue".')
% 
% propStrs = varargin(1:2:propargin);
% propValues = varargin(2:2:propargin);
% 
% for iParam = 1:propargin/2
%     switch lower(propStrs{iParam})
%         case lower('direction')
%             direction = propValues{iParam};
%         otherwise
%             error('optimal:simState:options',...
%               'Option string ''%s'' is not recognized.',propStrs{iParam})
%     end
% end
% 
% % Set to default value if necessary
% if ~exist('direction','var'), direction = 'forward'; end
% 
% % Check property values for errors TODO: Add property error checks
% assert(ismember(direction,{'forward','backward'}),...
%     'optimal:simState:direction',...
%     'Property "direction" must be either ''forward'' or ''backward''.')

%% Simulate
% if strcmp(direction,'forward')
%     [x,xD] = optimal.simulate(@(x_,t_,k_) f(x_,uFFunc(t_,u,t),t_),x0,t);
%     [x,xD] = optimal.simulate(@(x_,k_) f(x_,u(k_),t(k_)),x0,t);

n = size(x0,1);
x = nan(n,tn);
xD = nan(n,tn-1);
x(:,1) = x0;
for k = 1:tn-1
    xD(:,k) = f(x(:,k),u(:,k),t(k));
    ts = t(k+1) - t(k);
    x(:,k+1) = x(:,k) + xD(:,k)*ts;
end
% else
%     [x,xD] = optimal.simulate(@(x_,t_) f(x_,uBFunc(t_,u,t),t_),x0,t(end:-1:1));
% %     [x,xD] = optimal.simulate(@(x_,k_) f(x_,u(k_),t(k_)),x0,t);
%     x = x(:,end:-1:1);
%     xD = xD(:,end:-1:1);
% end

end

% function u = uFFunc(t,uTraj,tTraj)
% % k = find(t >= tTraj,1,'last');
% % u = uTraj(k);
% u = uTraj(diff(t >= tTraj) ~= 0);
% end
% 
% function u = uBFunc(t,uTraj,tTraj)
% k = find(t > tTraj,1,'last');
% u = uTraj(k);
% end
