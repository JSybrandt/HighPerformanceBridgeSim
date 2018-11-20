clc; clear all;
lim=361*2; % Number of days monitoring subject bridge
E_damaged=zeros(2,lim);

% Load variable arrays
BridgeVariables=load('BridgeVariables.dat');
E=BridgeVariables(1,3); % Modulus of elasticity of bridge N/m^2

for i=1:2
Damage_Case=i;
if Damage_Case==1
DayDamage1=round(lim*.25)+round(rand(1)*(lim*.33-lim*.25)); % The day damage is iniciated on bridge
ED=[.05*E,.05*ones(1,(lim-DayDamage1+1))*E; .1*E,.0025*ones(1,(lim-DayDamage1+1))*E]; % Damaged Modulus 1

for Day=1:lim
E_damaged(:,Day)=E;
       if Day>=DayDamage1
           if sum(ED(2,1:(Day-DayDamage1+1)))/E >= .5
               E_damaged(1,Day)=E-.05*E;
               E_damaged(2,Day)=E*.5;
           else
E_damaged(1,Day)=E-.05*E;
E_damaged(2,Day)=(E-sum(ED(2,1:(Day-DayDamage1+1))));
           end
       end
end

figure(1)
set(gcf,'color','white')
subplot(2,1,2);
plot(1:lim,E_damaged(1,:),'b','linewidth',2);hold on
plot(1:lim,E_damaged(2,:),'--r','linewidth',2); 
[fillhandle,msg]=jbfill(1:lim,E_damaged(1,:),E_damaged(2,:),'g','g',1,.3);
xlabel('Time (days) ');
ylabel('Elastic Modulus (N/m^2)');
legend('Min Damage','Max Damage','location','southwest','FontName','Timesnewroman')
yyaxis right
plot(1:lim,E_damaged(1,:)/E,'b','linewidth',2,'HandleVisibility','off');
plot(1:lim,E_damaged(2,:)/E,'--r','linewidth',2,'HandleVisibility','off');
ylabel('% Health Remaining');
xlim([1 365*2])
plotformat

elseif Damage_Case==2
DayDamage=[round(lim*.25)+round(rand(1)*(lim*.3-lim*.25)),round(lim*.35)+round(rand(1)*(lim*.45-lim*.35)), round(lim*.55)+round(rand(1)*(lim*.65-lim*.55)), round(lim*.7)+round(rand(1)*(lim*.75-lim*.7)), round(lim*.8)+round(rand(1)*(lim*.9-lim*.8))];
ED=[.05*ones(1,5)*E; .1*ones(1,5)*E]; % Damaged Modulus 1 

for Day=1:lim
   E_damaged(:,Day)=E;
   
    if Day>=DayDamage(5)
E_damaged(:,Day)=E-ED(:,1)-ED(:,2)-ED(:,3)-ED(:,4)-ED(:,5); % Overall Damaged Modulus 1 
    elseif Day>=DayDamage(4)
E_damaged(:,Day)=E-ED(:,1)-ED(:,2)-ED(:,3)-ED(:,4); % Overall Damaged Modulus 1
    elseif Day>=DayDamage(3)
E_damaged(:,Day)=E-ED(:,1)-ED(:,2)-ED(:,3); % Overall Damaged Modulus 1 
    elseif Day>=DayDamage(2)
E_damaged(:,Day)=E-ED(:,1)-ED(:,2); % Overall Damaged Modulus 1  
    elseif Day>=DayDamage(1)
E_damaged(:,Day)=E-ED(:,1); % Overall Damaged Modulus 1
    end  
end
subplot(2,1,1);
set(gcf,'color','white')
plot(1:lim,E_damaged(1,:),'b','linewidth',2);hold on
plot(1:lim,E_damaged(2,:),'--r','linewidth',2);
[fillhandle,msg]=jbfill(1:lim,E_damaged(1,:),E_damaged(2,:),'g','g',1,.3);
ylabel('Elastic Modulus (N/m^2)');
legend('Min Damage','Max Damage','location','southwest','FontName','Timesnewroman')
yyaxis right
plot(1:lim,E_damaged(1,:)/E,'b','linewidth',2,'HandleVisibility','off');
plot(1:lim,E_damaged(2,:)/E,'--r','linewidth',2,'HandleVisibility','off');
ylabel('% Health Remaining');
xlim([1 365*2])
plotformat

end

end




    
    