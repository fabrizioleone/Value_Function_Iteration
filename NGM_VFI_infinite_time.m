% NGM model with value function iteration
% Code for infinite time DP 
% we set a initial guess for our value and we iterate forward

% Fabrizio Leone - 26-01-2018

%% 1. Housekeeping
clear all
clc

%% 2. Set parameters
sigma=1;
beta=0.96;
alpha=0.33;
delta=0.04;
kstar=(((1/(alpha*beta))-((1-delta)/alpha)))^(1/(alpha-1));

%create a grid for capital
kmin=kstar*.9;
kmax=kstar*1.1;
n=200;
step=(kmax-kmin)/(n-1);
k=(kmin:step:kmax);

%% 3. Set variables

y=k.^alpha; %production function
ytot=repmat(y+(1-delta)*k,n,1);
kprime=repmat(k,n,1);
c=ytot-kprime'; %check how consumption is constructed

%define utility function

if sigma==1
u=log(c);
else
u=(c.^(1-sigma)-1)/(1-sigma);
end
u(c<0) = -Inf;

% Initialize VFI
V0=zeros(n,1); %initial value function guess
diff=10;
toler=10^(-5);
it=1;
maxit=10^10;


while diff>=toler & it<=maxit %it<=maxit is not really necessary
    
    vp=repmat(V0,1,n);
    v=u+beta*vp;
    [V1 p]=max(v,[],1); %compute value and policy function (p): search for the max along each column
    value(it,:)=V1; %store the value of each iteration
    policy(it,:)=p;%store the policy of each iteration
    
    %metric distance
    diff=max(abs(V1'-V0)); 
    
    %update value function
    V0=V1';
    it=it+1;
    disp(it)
    
    % when the loop breaks, V1 will be the value function which solves the
    % problem: V_t=V_{t+1}. it is the fixed point.
    
end

%% 4. Plotting results

figure(1)
plot(k,V1)
title('Value Function against capital stock')
xlabel('capital stock')
ylabel('value function')
%notice: the value function is increasing in the current stock of capital

figure(2)
plot(k,p);
title('Policy function against capital stock')
xlabel('capital stock')
ylabel('value function')
% 
% figure(3)
% surf(value)
% title('Value Function surface')
% xlabel('capital today')
% ylabel('capital tomorrow')
% zlabel('value function')

%% 5. Simulation of transition dynamics

   P=100; % arbitrary length for transition path
   capital_index=ones(1,P);
   capital_transition=ones(1,P);
   capital_index(1)=(3); % arbitrary starting point in terms of the index
   %capital_index(1)=(150); % set an initial value above the ss.
   capital_transition(1)= k(capital_index(1));
   
   for t=2:P
   capital_index(t)=(p(capital_index(t-1)));%  evolution in index space   
   capital_transition(t)=k(capital_index(t)); % evoluation in capital space
   end
    
figure(4)    
plot(capital_transition)
title('Transitional dynamics for capital')
xlabel('time period')
ylabel('capital stock')
