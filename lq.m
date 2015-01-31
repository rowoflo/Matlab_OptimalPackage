function [P,W] = lq(t,A,B,Q,R,S,q,r)
% The "lq" function the general Linear Quadratic problem: 
%   min_u J = \int_0^T 1/2 x^\T Q x + 1/2 u^\T R u + q^\T x + r^\T u dt + 1/2 x^\T(tf) S x(tf) 
%   s.t
%   \dot{x} = A x + B u where Q = Q^\T >= 0, R = R^\T > 0, S = S^\T >= 0
%   with the solution as
%   u = -R^{-1} B^\T P x - R^{-1} (B^\T w + r)
%
% SYNTAX:
%   [P,W] = lq(t,A,B,Q,R,S)
%   [P,W] = lq(t,A,B,Q,R,S,q,r)
% 
% INPUTS:
%   t - (1 x tn)
%       Time trajectory.
%
%   A - (n x n or n x n x tn or function_handle)
%       State to state dynamics matrix.
%       * 2 dimensional matrix means time-invariant.
%       * 3 dimensional specifies matrix at each time point.
%       * Function handle: A(k) returns n x n matrix for given time index.
%
%   B - (n x m or n x m x tn or function_handle)
%       Input to state dynamics matrix.
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
%   S - (n x n number)
%       Semi-positive definite symmetric weighting matrix on the quadratic
%       term of the final state.
%
%   q - (n x 1 or n x tn or function_handle) [zeros(n,1)]
%       Weighting vector on the linear term of the state.
%       * 1 dimensional matrix means time-invariant.
%       * 2 dimensional specifies matrix at each time point.
%       * Function handle: q(k) returns n x 1 vector for given time index.
%
%   r - (m x 1 or n x tn or function_handle) [zeros(m,1)]
%       Weighting vector on the linear term of the input.
%       * 1 dimensional matrix means time-invariant.
%       * 2 dimensional specifies matrix at each time point.
%       * Function handle: r(k) returns m x 1 vector for given time index.
% 
% OUTPUTS:
%   P - (n x n x tn number) 
%       Quadratic term in the solution to Riccati equation.
%
%   W - (n x tn number) 
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
narginchk(6,8)

% Check input arguments for errors
assert(isnumeric(t) && isreal(t) && isvector(t) && all(diff(t) > 0),...
    'optimal:lq:t',...
    'Input argument "t" must be a vector of increaseing real numbers.')
t = t(:)';
tn = size(t,2);

if isa(A,'function_handle')
    fA = A;
    n = size(A(1),1);
else
    assert(isnumeric(A) && isreal(A) && size(A,1) == size(A,2),...
        'optimal:lq:A',...
        'Input argument "A" must be a n x n matrix of real numbers.')
    n = size(A,1);
    if length(size(A)) == 3
        assert(size(A,3) == tn,...
            'optimal:lq:A',...
            'Input argument "A" must be a n x n x %d matrix of real numbers.',tn)
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
    assert(isnumeric(Q) && isreal(Q) && isequal(size(Q),[n n]) && all(Q == Q') && all(real(eig(Q)) >= 0),...
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
    assert(isnumeric(R) && isreal(R) && isequal(size(R),[m m]) && all(R == R') && all(real(eig(Q)) > 0),...
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

assert(isnumeric(S) && isreal(S) && isequal(size(S),[n n]) && all(S == S') && all(real(eig(S)) >= 0),...
    'optimal:lq:S',...
    'Input argument "S" must be a %d x %d semi-positive definite symmetric matrix of real numbers.',n,n)

if nargin < 7, q = zeros(n,1); end
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

if nargin < 8, r = zeros(m,1); end
if isa(r,'function_handle')
    fr = r;
else
assert(isnumeric(r) && isreal(r) && size(r,1) == n && ...
    (size(r,2) == tn || size(r,2) == 1),...
    'optimal:lq:r',...
    'Input argument "r" must be a %d x 1 vector or %d x %d matrix of real numbers.',m,m,tn)
    if size(r,2) == tn
        fr = @(k_) r(:,k_);
    else
        fr = @(k_) r(:,1);
    end
end

%% Solve
g = @(~,~,~) 0;

fP = @(P_,~,~,k_) reshape(-reshape(P_,[n n])*fA(k_) - fA(k_)'*reshape(P_,[n n]) - fQ(k_) + reshape(P_,[n n])*fB(k_)*fR(k_)^-1*fB(k_)'*reshape(P_,[n n]),[n*n 1]);
P = reshape(optimal.simulate(fP,g,t,S,-inf*ones(n,1),inf*ones(n,1),'direction','backward'),[n n tn]);

fW = @(W_,~,~,k_) (P(:,:,k_)*fB(k_)*fR(k_)^-1*fB(k_)' - fA(k_)')*W_ - fq(k_) + P(:,:,k_)*fB(k_)*fR(k_)^-1*fr(k_);
W = optimal.simulate(fW,g,t,zeros(n,1),-inf*ones(n,1),inf*ones(n,1),'direction','backward');

end
