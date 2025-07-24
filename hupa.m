% HupA gene regulatory network in V. vulnificus from 
% "Switch-like behavior in the heme receptor for Vibrio vulnificus"

% reproduces the dotted blue curve from Figure 5

%% set up

% nondimensionalized parameters values
extFe=0.1;      % external iron
eps=0.01;       % leakage
ft=1;           % total amount of Fur
a1=1;           % time scale constant for HupA
a2=1;           % time scale constant for HupR
a3=10;          % creation rate of external heme
a4=1;           % importation rate of heme via HupA
a4m=10;         % iron degradation / usage rate
k1=1;           % binding affinity of fur and fe
k2=10;          % importation rate of external fe
k3=1;           % strength of fur inhibition of iron importation
b1=4;           % binding affinity for HupR/heme and DNA binding site
b2=2;           % binding affinity for f* and DNA binding site

params = [a1,a2,a3,a4,a4m,k1,k2,k3,b1,b2,ft,eps,extFe]; % vector to pass to ode 45
 
% simulation time
t0=0;                   % initial / starting time
tfinal=10000;           % final time
tvec=t0:0.01:tfinal;    % vector of times to simulate over

% set up conditions for varying external heme over time
h_max=1.5;  % max h_ext value (minimum set at 0)
stp=0.0015; % steepness of transition 
t_mid1=(t0+tfinal)/4;   % quarter point along total simulation fime
t_mid2=3*(t0+tfinal)/4; % 3/4 point along total simulation time

% function to describe smooth transition between 0 and h_max
var_heme = @(t) h_max/2*tanh(stp*(t-t_mid1))-h_max/2*tanh(stp*(t-t_mid2));
  
   
x0=[0.01,0.2,0.1,1,1];  % initial conditions
 
%% simulation and plots

% numerically simulate system of odes  
[t,x] = ode45(@(t,x)system(t,x,params,var_heme),tv,x0);


% Plots 

figure(1)   % not included in the paper
clf
% plots bifurcation variable (extHeme) and variable of interest (HupA)
% against time

hold on
plot1=plot(tvec,var_heme(tvec),'-','linewidth',2,'Color',"#0072BD");  % plots variable external heme
plot2=plot(tvec,x(:,1),'linewidth',2,'Color','#D95319');               % plots HupA
hold off

% plot formatting
legend([plot1(1),plot2(1)],{'$h_{ext}$','HupA'},'Interpreter','latex') 
xlabel('t','Interpreter','latex')
ylabel('conc.','Interpreter','latex')
fontsize(16,"points")
set(gca,'TickLabelInterpreter','latex')
box on
hold off



figure(2)
clf 
% plots trajectory in HupA for variable heme
% blue dashed curve in Figure 5

hold on
plot1=plot(var_heme(tvec),x(:,1),'--','linewidth',2,'Color',"#0072BD",'DisplayName','trajectory');

% plot formatting
xlim([0,1])
ylim([0,0.4])
xlabel('$h_{ext}$','Interpreter','latex')
ylabel('$a$','Interpreter','latex')
fontsize(16,"points")
set(gca,'TickLabelInterpreter','latex')
box on
hold off

function dxdt = system(t,x,p,var_heme)
% system of ODEs (Eq. 6-10)

% sets variable value of heme each time function is called
extHe=var_heme(t);

% ODEs
dadt=p(1)*((p(9)*x(2)^4*x(3)^2)/(p(9)*x(2)^4*x(3)^2+p(10)*x(4)+1)+p(12)-x(1));  %HupA
drdt=p(2)*((1)/(p(9)*x(2)^4*x(3)^2+p(10)*x(4)+1)- x(2));    % HupR
dhidt=p(3)*extHe/(1+extHe)*x(1)-x(3);       % heme
dfedt= p(4)*x(3) - p(5)*x(5) + p(7)*1/(1+p(8)*x(4))*p(13)/(1+p(13)) - p(6)*((p(11)-x(4))*x(5)-x(4));    % iron
dfsdt=p(6)*((p(11)-x(4))*x(5)-x(4));        % activated fur

% pass all odes to single vector
dxdt = [dadt drdt dhidt dfsdt dfedt]';
end