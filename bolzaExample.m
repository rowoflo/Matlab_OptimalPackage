%% bolzaExample.m
% The "bolzaExample" script is example code for how to use the
% optimal.bolza function.
%
% NOTES:
%
% NECESSARY FILES AND/OR PACKAGES:
%   +optimal, bolza.m
%
% SEE ALSO: TODO: Add see alsos
%    relatedFunction1 | relatedFunction2
%
% AUTHOR:
%    Rowland O'Flaherty (www.rowlandoflaherty.com)
%
% VERSION: 
%    Created 02-MAY-2014
%-------------------------------------------------------------------------------

%% Clear
clear; clf; clc;

%% Import
import optimal.*

%% Initialize

% Time - parameters
ts = .1; % (1 x 1) Time step size
t0 = 0; % (1 x 1) Initial time
tf = 10; % (1 x 1) Final time
% Time - variables
t = t0:ts:tf; % (1 x tn) Time vector record for all time

% State - parameters
x0 = 0; % (n x 1) Initial state
xBar = 10; % (n x 1) Desired state

% Dynamics
f = @(x_,u_,t_) u_; % (n x 1) State dynamics (i.e. xDot)
dfdx = @(x_,u_,t_) zeros(size(x_,1),size(x_,1),size(t_,2)); % (n x n) Dynamics partial to state
dfdu = @(x_,u_,t_) ones(size(x_,1),size(u_,1),size(t_,2)); % (n x m) Dynamics partial to input

% Cost
rho = 100; % (1 x 1) Final cost weight
L = @(x_,u_,t_) u_.^2; % (1 x 1) Instantaneous cost
dLdx = @(x_,u_,t_) zeros(1,size(x_,1),size(t_,2)); % (1 x n) Instantaneous cost partial to state
dLdu = @(x_,u_,t_) 2*permute(u_,[1,3,2]); % (1 x m) Instantaneous cost partial to input
Psi = @(xf_,tf_) rho*(xf_ - xBar).^2; % (1 x 1) Final cost
dPsidx = @(xf_,tf_) 2*rho*(xf_ - xBar); % (1 x n) Final cost partial to final state

% Armijo parameters
alpha = 0.25;
beta = 0.75;

% Solve Bolza
[x,u,lambda,JTape,gammaTape] = bolza(t,x0,f,dfdx,dfdu,L,dLdx,dLdu,Psi,dPsidx,'armijoParams',[alpha beta]);

%% Output
% Data
fprintf('Final cost: %.3f\n',JTape(end));
fprintf('Final state: %.3f\n',x(end));

% Plot
figure(1)
subplot(3,1,1)
plot(t(1:end-1),u)
title('Input Trajectory')
xlabel('Time')
ylabel('Input')
grid on
subplot(3,1,2)
plot(t,x)
hold on
plot(t,repmat(xBar,size(t)),'r--')
hold off
title('State Trajectory')
xlabel('Time')
ylabel('State')
grid on
subplot(3,1,3)
plot(t,lambda)
title('Costate Trajectory')
xlabel('Time')
ylabel('Costate')
grid on

figure(2)
subplot(2,1,1)
plot(JTape)
title('Cost Profile')
xlabel('Iteration')
ylabel('Cost')
subplot(2,1,2)
plot(gammaTape)
title('Step Size Profile')
xlabel('Iteration')
ylabel('Step Size')


