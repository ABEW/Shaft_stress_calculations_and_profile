function Figuretastic( Func,Limits,Fig_Title, x_label, y_label )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',16);

box(axes1,'on');
hold(axes1,'all');

% Create plot
fplot(Func,Limits)
set(findobj(gca,'Type','Line'),'LineWidth',2)

% Create xlabel
xlabel(x_label,'FontWeight','bold','FontSize',16);

% Create ylabel
ylabel(y_label,'FontWeight','bold','FontSize',16);

% Create title
T=title(Fig_Title);
set(T,'FontSize',16,'FontWeight','bold')


end

