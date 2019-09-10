% NGM model with value function iteration
% Code for finite time DP 
% we set V_{T+1}=0 and solve backwards

% Fabrizio Leone - 03-02-2018

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
T=150; %number of periods. so T-1 is our last choice period 
step=(kmax-kmin)/(T-1);
k=(kmin:step:kmax);

%% 3. Set variables

y=k.^alpha; %production function
ytot=repmat(y+(1-delta)*k,T,1);
kprime=repmat(k,T,1);
c=ytot-kprime'; %check how consumption is constructed

%define utility function

if sigma==1
u=log(c);
else
u=(c.^(1-sigma)-1)/(1-sigma);
end
u(c<0) = -Inf;

% Initialize VFI
V0=zeros(T,T); %initial value function guess
value=zeros(T,T);
policy=zeros(T,T);

for t=1:T-1 
    
    
    v=u+beta*V0;
    [V1 p]=max(v,[],1); %compute value and policy function (p): search for the max along each column
    value(:,T-t)=V1'; %store the value of each iteration
    policy(:,T-t)=p';%store the policy of each iteration
    V0=V1';

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
