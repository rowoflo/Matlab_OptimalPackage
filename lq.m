function [x,u,xDot,P,W] = lq(t,x0,f,A,B,Q,R,S,q,r,xr,ur,varargin)
% The "lq" function solves the general Linear Quadratic problem: 
%   min_u J = \int_0^T 1/2 x^\T Q x + 1/2 u^\T R u + q^\T x + r^\T u dt + 1/2 x^\T(tf) S x(tf) 
%   s.t
%   \dot{x} = A x + B u where Q = Q^\T >= 0, R = R^\T > 0, S = S^\T >= 0
%   with the solution as
%   u = -R^{-1} B^\T P (x - xr) - R^{-1} (B^\T w + r) + ur
%
% SYNTAX:
%   [x,u,P,W] = optimal.lq(t,x0,A,B,Q,R,S,q,r,xr,ur)
%   [x,u,P,W] = optimal.simulate(t,x0,A,B,Q,R,S,q,r,'PropertyName',PropertyValue,...)
% 
% INPUTS:
%   t - (1 x tn)
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
%   A - (n x n or n x n x tn or function_handle)
%       State linearization at the reference state.
%       * 2 dimensional matrix means time-invariant.
%       * 3 dimensional specifies matrix at each time point.
%       * Function handle: A(k) returns n x n matrix for given time index.
%
%   B - (n x m or n x m x tn or function_handle)
%       Input linearization at the referenece input.
%       * 2 dimensional matrix means time-invariant.
%       * 3 dimensional specifies matrix at each time point.
%       * Function handle: B(k) returns n x m matrix for given time index.
%
%   Q - (n x n or n x n x tn or function_handle)
%       Semi-positive definite symmetric weighting matrix on the quadratic
%       term of the state.
%       * 2 dimensional matrix means time-invariant.
%       * 3 dimensional specifies matrix at each time point.
%       * Function handle: Q(k) returns n x n matrix for given time index.
%
%   R - (m x m or m x m x tn or function_handle)
%       Positive definite symmetric weighting matrix on the quadratic term
%       of the input.
%       * 2 dimensional matrix means time-invariant.
%       * 3 dimensional specifies matrix at each time point.
%       * Function handle: R(k) returns m x m matrix for given time index.
%
%   S - (n x n)
%       Semi-positive definite symmetric weighting matrix on the quadratic
%       term of the final state.
%
%   q - (n x 1 or n x tn or function_handle) [zeros(n,1)]
%       Weighting vector on the linear term of the state.
%       * 1 dimensional vector means time-invariant.
%       * 2 dimensional matrix specifies vector each time point.
%       * Function handle: q(k) returns n x 1 vector for given time index.
%
%   r - (m x 1 or n x tn or function_handle) [zeros(m,1)]
%       Weighting vector on the linear term of the input.
%       * 1 dimensional vector means time-invariant.
%       * 2 dimensional matrix specifies vector each time point.
%       * Function handle: r(k) returns m x 1 vector for given time index.
%
%   xr - (n x 1 or n x tn or function_handle) [zeros(n,1)]
%       Reference state trajectory, which the system is linearized around.
%       * 1 dimensional vector means time-invariant.
%       * 2 dimensional matrix specifies vector each time point.
%       * Function handle: q(k) returns n x 1 vector for given time index.
%
%   ur - (m x 1 or m x tn or function_handle) [zeros(m,1)]
%       Reference inpout trajectory, which the system is linearized around.
%       * 1 dimensional vector means time-invariant.
%       * 2 dimensional matrix specifies vector at each time point.
%       * Function handle: r(k) returns m x 1 vector for given time index.
%
% PROPERTIES:
%   'statemin' - (n x 1 number) [-inf*ones(n,1)]
%       Minimum state constraint.
%
%   'statemax' - (n x 1 number) [inf*ones(n,1)]
%       Maximum state constraint.
%
%   'inputmin' - (m x 1 number) [-inf*ones(m,1)]
%       Minimum input constraint.
%
%   'inputmax' - (m x 1 number) [inf*ones(m,1)]
%       Maximum input constraint.
%
% 
% OUTPUTS:
%   x - (n x tn)
%       Optimal state trajectory.
%
%   u - (m x tn)
%       Optimal input trajectory.
%
%   xDot - (n x tn)
%       Optimal state time derivative trajectory.
%
%   P - (n x n x tn) 
%       Quadratic term in the solution to Riccati equation.
%
%   W - (n x tn) 
%       Linear term in the solution to Riccati equation.
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
%   Created 28-JAN-2015
%-------------------------------------------------------------------------------

%% Initialize

% Check number of inputs
narginchk(12,inf)

% Check input arguments for errors
assert(isnumeric(t) && isreal(t) && isvector(t) && all(diff(t) > 0),...
    'optimal:lq:t',...
    'Input argument "t" must be a vector of increaseing real numbers.')
t = t(:)';
tn = size(t,2);

assert(isnumeric(x0) && isreal(x0) && isvector(x0),...
    'optimal:lq:x0',...
    'Input argument "x0" must be a vector of real numbers.')
x0 = x0(:);
n = size(x0,1);

assert(isa(f,'function_handle'),...
    'optimal:simulate:f',...
    'Input argument "f" must be a function handle.')

if isa(A,'function_handle')
    fA = A;
else
    assert(isnumeric(A) && isreal(A) && size(A,1) == size(A,2),...
        'optimal:lq:A',...
        'Input argument "A" must be a %d x %d matrix of real numbers.',n,n)
    if length(size(A)) == 3
        assert(size(A,3) == tn,...
            'optimal:lq:A',...
            'Input argument "A" must be a %d x %d x %d matrix of real numbers.',n,n,tn)
        fA = @(k_) A(:,:,k_);
    else
        fA = @(k_) A(:,:,1);
    end
end

if isa(B,'function_handle')
    fB = B;
    m = size(B(1),2);
else
    assert(isnumeric(B) && isreal(B) && size(B,1) == n,...
        'optimal:lq:B',...
        'Input argument "B" must be a %d x m matrix of real numbers.',n)
    m = size(B,2);
    if length(size(B)) == 3
        assert(size(B,3) == tn,...
            'optimal:lq:B',...
            'Input argument "B" must be a %d x m x %d matrix of real numbers.',n,tn)
        fB = @(k_) B(:,:,k_);
    else
        fB = @(k_) B(:,:,1);
    end
end

if isa(Q,'function_handle')
    fQ = Q;
else
    assert(isnumeric(Q) && isreal(Q) && isequal(size(Q),[n n]) && all(all(Q == Q')) && all(real(eig(Q)) >= 0),...
        'optimal:lq:Q',...
        'Input argument "Q" must be a %d x %d semi-positive definite symmetric matrix of real numbers.',n,n)
    if length(size(Q)) == 3
        assert(size(Q,3) == tn,...
            'optimal:lq:Q',...
            'Input argument "Q" must be a %d x %d x %d matrix of real numbers.',n,n,tn)
        fQ = @(k_) Q(:,:,k_);
    else
        fQ = @(k_) Q(:,:,1);
    end
end

if isa(R,'function_handle')
    fR = R;
else
    assert(isnumeric(R) && isreal(R) && isequal(size(R),[m m]) && all(all(R == R')) && all(real(eig(R)) > 0),...
        'optimal:lq:R',...
        'Input argument "R" must be a %d x %d positive definite symmetric matrix of real numbers.',m,m)
    if length(size(R)) == 3
        assert(size(R,3) == tn-1,...
            'optimal:lq:R',...
            'Input argument "R" must be a %d x %d x %d matrix of real numbers.',m,m,tn)
        fR = @(k_) R(:,:,k_);
    else
        fR = @(k_) R(:,:,1);
    end
end

assert(isnumeric(S) && isreal(S) && isequal(size(S),[n n]) && all(all(S == S')) && all(real(eig(S)) >= 0),...
    'optimal:lq:S',...
    'Input argument "S" must be a %d x %d semi-positive definite symmetric matrix of real numbers.',n,n)

if isa(q,'function_handle')
    fq = q;
else
    assert(isnumeric(q) && isreal(q) && size(q,1) == n && ...
        (size(q,2) == tn || size(q,2) == 1),...
        'optimal:lq:q',...
        'Input argument "q" must be a %d x 1 vector or %d x %d matrix of real numbers.',n,n,tn)
    if size(q,2) == tn
        fq = @(k_) q(:,k_);
    else
        fq = @(k_) q(:,1);
    end
end

if isa(r,'function_handle')
    fr = r;
else
assert(isnumeric(r) && isreal(r) && size(r,1) == m && ...
    (size(r,2) == tn || size(r,2) == 1),...
    'optimal:lq:r',...
    'Input argument "r" must be a %d x 1 vector or %d x %d matrix of real numbers.',m,m,tn)
    if size(r,2) == tn
        fr = @(k_) r(:,k_);
    else
        fr = @(k_) r(:,1);
    end
end

if isa(xr,'function_handle')
    fxr = xr;
else
    assert(isnumeric(xr) && isreal(xr) && size(xr,1) == n && ...
        (size(xr,2) == tn || size(xr,2) == 1),...
        'optimal:lq:xr',...
        'Input argument "xr" must be a %d x 1 vector or %d x %d matrix of real numbers.',n,n,tn)
    if size(xr,2) == tn
        fxr = @(k_) xr(:,k_);
    else
        fxr = @(k_) xr(:,1);
    end
end

if isa(ur,'function_handle')
    fur = ur;
else
assert(isnumeric(ur) && isreal(ur) && size(ur,1) == m && ...
    (size(r,2) == tn || size(r,2) == 1),...
    'optimal:lq:ur',...
    'Input argument "ur" must be a %d x 1 vector or %d x %d matrix of real numbers.',m,m,tn)
    if size(ur,2) == tn
        fur = @(k_) ur(:,k_);
    else
        fur = @(k_) ur(:,1);
    end
end

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

%% Solve
g = @(~,~,~) 0;

fP = @(P_,~,~,k_) reshape(-reshape(P_,[n n])*fA(k_) - fA(k_)'*reshape(P_,[n n]) - fQ(k_) + reshape(P_,[n n])*fB(k_)*fR(k_)^-1*fB(k_)'*reshape(P_,[n n]),[n*n 1]);
P = reshape(optimal.simulate(t,reshape(S,[n*n 1]),fP,g,'direction','backward'),[n n tn]);

if any(any(fq(1:tn) ~= 0)) || any(any(fr(1:tn) ~= 0))
    fW = @(W_,~,~,k_) (P(:,:,k_)*fB(k_)*fR(k_)^-1*fB(k_)' - fA(k_)')*W_ - fq(k_) + P(:,:,k_)*fB(k_)*fR(k_)^-1*fr(k_);
    W = optimal.simulate(t,zeros(n,1),fW,g,'direction','backward');
    g = @(x_,t_,k_) -fR(k_)^-1*fB(k_)'*P(:,:,k_)*(x_ - fxr(k_)) - fR(k_)^-1*(fB(k_)'*W(:,k_) + fr(k_)) + fur(k_);
else
    g = @(x_,t_,k_) -fR(k_)^-1*fB(k_)'*P(:,:,k_)*(x_ - fxr(k_)) + fur(k_);
end
[x,u,xDot] = optimal.simulate(t,x0,f,g,'statemin',xm,'statemax',xM,'inputmin',um,'inputmax',uM);

end
