function plotformat
set(gca,'FontSize',12,'FontWeight','bold','linewidth',2,'FontName','Times New Roman') 
set(gca,'fontsize',22,'FontName','Times New Roman');




%% Y-Label
A=get(gca,'YLabel');
a=get(A,'String');

if ~isempty(a)
    ylabel(a,'fontsize',28,'FontName','Times New Roman')
end

%% X-Label
B=get(gca,'XLabel');
b=get(B,'String');

if ~isempty(b)
    xlabel(b,'fontsize',28,'FontName','Times New Roman')
end

%% Title
%get the title
C=get(gca,'Title');
c=get(C,'String');

%new title
if ~isempty(c)
    title(c,'fontsize',30,'FontName','Times New Roman')
end

% %% Legend
% D=get(gca,'Legend');
% d=get(D,'String');
% 
% if ~isempty(d)
%     legend(d,'fontsize',22,'fontname','Calibri')
% end
end