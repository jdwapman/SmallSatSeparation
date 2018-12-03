close all

%% Color map with red and blue on the edges

C = distinguishable_colors(N);

figure
subplot(2,2,1)


for nplot=1:1:N
polarplot(thetaActT(nplot), rActT(nplot),'o','color',C(nplot,:));
hold on 
end
title("Actual Final State")
% Find the linearized final state
subplot(2,2,2)
for nplot=1:1:N
polarplot(thetaLinT(nplot), rLinT(nplot),'o','color',C(nplot,:));
hold on
end
title("Linearized Final State")
legend; 

%% plotting the constelation by day actual states
subplot(2,2,3)
for  tplot = 1:1:T
for nplot= 1:1:N 

if tplot==T % mark differently for end state
  polarplot(thetaAct(nplot,tplot), 100+T-tplot, 'x','color',C(nplot,:))
else 
polarplot(thetaAct(nplot,tplot), 100+T-tplot, '.-','color',C(nplot,:))
end
hold on
title("Actual State")
end
pause(0.02);
end

%{
plotting the constelation by day linear states
subplot(2,2,4)
for  tplot = 1:1:T
for nplot= 1:1:N  
    if tplot==T % mark differently for end state
         polarplot(thetalin(nplot,tplot), 100+T-tplot, 'x','color',C(nplot,:))
  title("Linear State")
hold on
    else
       polarplot(thetalin(nplot,tplot),100+T-tplot, '.-','color',C(nplot,:))  
       hold on
    end
end
   pause(0.02);
end

%}



%% Plot inputs

figure
subplot(2,2,1)
t = 1:1:T;
plot(t, uOptReshape);

xlabel("Time (Days)")
ylabel('Area (m^2)')
title("Satellite Area vs Time")

%% Plot radius trajectory
subplot(2,2,2)
t = 1:1:T+1;
for nplot=1:1:N
plot(t, rAct(nplot,:)/1000,'-','color',C(nplot,:));
hold on
end
xlabel("Time (Days)")
ylabel("Altitude (km)")
title("Altitude")
subplot(2,2,3)
for nplot=1:1:N
plot(t, wAct(nplot,:)/1000,'-','color',C(nplot,:));
hold on
end
xlabel("Time (Days)")
ylabel("Angular Velocity (rad/sec)")
title("Angular Velocity ")






%% plot angular seperation
 AngSepRef=zeros(N,T);
 AngSepLin=zeros(N,T);
 AngSepAct=zeros(N,T);

 %clalc linear stuff
 rlin=r0;
 wlin=w0;
 thetalin=theta0;
 for n=1:1:N
    for t=1:1:T
    rlin(n,t+1)=rlin(n,t)+dt*Sr(rRef(n,t),wlin(n,t))*uOptReshape(n,t);
    wlin(n,t+1)=wlin(n,t)+dt*Somega(rRef(n,t),wlin(n,t))*uOptReshape(n,t);
    thetalin(n,t+1)=thetalin(n,t)+dt*wRef(n,t)+0.5*dt*dt*Somega(rRef(n,t),wlin(n,t))*uOptReshape(n,t);
    end
 end
 
 
for  tplot = 1:1:T
    for nplot= 1:1:N 
        if nplot==1 % if the first then do something different
            AngSepRef(nplot,tplot)=thetaRef(nplot,tplot)-thetaRef(N,tplot);
            AngSepLin(nplot,tplot)=thetalin(nplot,tplot)-thetalin(N,tplot);
            AngSepAct(nplot,tplot)=thetaAct(nplot,tplot)-thetaAct(N,tplot);
        else 
            AngSepRef(nplot,tplot)=thetaRef(nplot,tplot)-thetaRef(nplot-1,tplot);
            AngSepLin(nplot,tplot)=thetalin(nplot,tplot)-thetalin(nplot-1,tplot);
            AngSepAct(nplot,tplot)=thetaAct(nplot,tplot)-thetaAct(nplot-1,tplot);
        end  
    end
end


t = 1:1:T;
figure
 hold on
for nplot=2:1:N 
    subplot(1,2,1)
    hold on
 a1=plot(t, AngSepRef(nplot,t),':','color',C(nplot,:),'DisplayName','Reference');
 a2=plot(t, AngSepLin(nplot,t),'o:','color',C(nplot,:),'DisplayName','Linear');
 a3=plot(t, AngSepAct(nplot,t),'--','color',C(nplot,:),'DisplayName','Actual');
subplot(1,2,2) 
hold on
plot(t, wRef(nplot,t),':','color',C(nplot,:),'DisplayName','Reference');
 plot(t, wlin(nplot,t),'o:','color',C(nplot,:),'DisplayName','Linear');
 plot(t, wAct(nplot,t),'--','color',C(nplot,:),'DisplayName','Actual');
     legend('Ref','Linear','Act');
end
xlabel("Time (Days)")
ylabel("Angular seperation (rad)")
title("Angular seperation ")


