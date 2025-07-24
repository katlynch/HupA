% HupA gene regulatory network in V. vulnificus from 
% "Switch-like behavior in the heme receptor for Vibrio vulnificus"


% simplified model with constant f*
% plots bifurcation graphs from Figure 3

clear all

%% Set Constants 

% model parameters
b2=1; eps=0.01; b1=1; a3=1;

% values of a overwhich we will calculate the solution curve in V(h_ext)
a = linspace(0,1,10000);

%% Reduced Model

% set up / clean figure
figure(1)
clf
hold on

% values of f* to iterate through
fs_vec = [0, 0.1, 0.5, 1, 2];
for i = 1:length(fs_vec)
    fs = fs_vec(i);     % set current value of f*
    vh = -(b2*fs + 1)^3*(a - eps)./(b1*a3*a.*(a - eps - 1).^3);    % calculate corresponding value of V(h_ext)
    plot(vh,a,'linewidth',2,'DisplayName',num2str(fs))       % plot curve
end

% plot formatting
lgd=legend('Location','northwest','Interpreter','latex');   
xlim([1e-2,1e6])                                            
ylim([0,1])                                                 
xlabel('$V(h_{ext})$','Interpreter','latex')
ylabel('$a$','Interpreter','latex')                         
title(lgd,'$f^*$','Interpreter','latex')
fontsize(16,"points")
set(gca, 'XScale', 'log')
set(gca,'TickLabelInterpreter','latex')
box on

%% n=2

% open and clear new figure
figure(2)
clf
hold on

% defines colors in the same order as the previous plot
clrs = ['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E';'#77AC30'; '#EDB120';];


for i=1:length(fs_vec)  % iterate through the same values of f* as above
    fs = fs_vec(i);     % set current value of f*
    vh2 = -(b2*fs + 1)^5*(a - eps)./(b1*(a - eps - 1).^5*a3^2.*a.^2); % calculate corresponding value of V(h_ext)
    
    
    % we now find limit points so that the unstable and stable regions can
    % be differentiated on the plot
    dx=diff(vh2); % finds difference between consecutive points
    idx=find(dx<0); % location of places where dx is negative, i.e.  curve changes directions at a limit point
    lp1=[vh2(idx(1)),a(idx(1))];    % first limit point
    lp2=[vh2(idx(end)),a(idx(end))];    % second limit point
    
     
    X = [vh2;a]';        % pairs a values (set above) with calculated v(h_ext)
    Xstable_l=X(X(:,1)<lp1(1) & X(:,2)<lp1(2),:);   % pulls out left stable branch
    Xstable_r=X(X(:,1)>lp2(1) & X(:,2)>lp2(2),:);   % right stable branch
    Xunstable=X(X(:,1)<lp1(1) & X(:,1)>lp2(1) & X(:,2)>lp1(2) & X(:,2)<lp2(2),:);   % unstable branch
    
    % plot all three branches with appropriate color and line style
    plot(Xstable_l(:,1),Xstable_l(:,2),'linewidth',2,'color',clrs(i,:),'DisplayName',num2str(fs))
    plot(Xstable_r(:,1),Xstable_r(:,2),'linewidth',2,'color',clrs(i,:),'HandleVisibility','off')
    plot(Xunstable(:,1),Xunstable(:,2),'--','linewidth',2,'color',clrs(i,:),'HandleVisibility','off')
  
end


% overlays the bottom line with one color for visual clarity with log scale
line([1e-10,10^2],[eps,eps],'linewidth',2,'Color','#77AC30','HandleVisibility','off')


% plot formatting
lgd=legend('Location','northwest','Interpreter','latex');
xlim([1e-2,1e6])
ylim([0,1])
xlabel('$V(h_{ext})^2$','Interpreter','latex')
ylabel('$a$','Interpreter','latex')
title(lgd,'$f^*$','Interpreter','latex')
fontsize(16,"points")
set(gca, 'XScale', 'log')
set(gca,'TickLabelInterpreter','latex')
box on


%% two parameter diagram
% for the n=2 case (four hupr two heme)
% plots cusp in fs and hext (implicit)

% opern and clear new figure
figure(3)
clf

hold on 

% sets values of f* for which we calculate the location of saddle nodes
fs_vec = linspace(0,10,10000);

% calculate the right limit point
lp1=-27648*(6180504276*sqrt(426) - 383599420199)*(fs_vec + 1).^5/536067676053505;
% calculate location of the left limit point
lp2=27648*(6180504276*sqrt(426) + 383599420199)*(fs_vec + 1).^5/536067676053505;

% plot both saddle node curves
plot(lp1,fs_vec,'linewidth',2,'color','r');
plot(lp2,fs_vec,'linewidth',2,'color','r');
hold off



% plot formatting
ylim([0,5])
xlim([10e-2,1e6])
xlabel('$V(h_{ext})^2$','Interpreter','latex')
ylabel('$f^*$','Interpreter','latex')
fontsize(16,"points")
set(gca, 'XScale', 'log')
set(gca,'TickLabelInterpreter','latex')
text(1e5,1.5,'I','Interpreter','latex','FontSize',20)   % labels right monostable region
text(1e0,1.5,'I','Interpreter','latex','FontSize',20)   % labels left monostable region
text(4.5e2,1.5,'II $\rightarrow$','Interpreter','latex','FontSize',20)  % labels bistable region
box on




