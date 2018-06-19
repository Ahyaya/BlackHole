function tidal(M,r0,a,v)
%输入格式tidal(M,r0,α,v)
%M为史瓦西黑洞质量（几何制）
%r0为粒子入射位置（几何制）
%α为粒子入射方向与半径方向的夹角（角度制）
%v为粒子入射时当地静态观察者测量的速度（几何制）
risk=0;
h=r0*sqrt(pi*r0/M)/1000;clr=[1 0 0;0 0 0;0 0 1];
tp(1)=0;
t(1)=0;
r(1)=r0;
fi(1)=0;
arf=a./180.*pi;
dt(1)=sqrt(1/(1-2*M/r0)/(1-v*v));
dr(1)=-v*cos(arf)*(1-2*M/r0)*dt(1);
dfi(1)=v*sin(arf)/r0*sqrt(1-2*M/r0)*dt(1);

ddt(1)=-2*M/r(1)/r(1)/(1-2*M/r(1))*dt(1)*dr(1);
ddr(1)=-M/r(1)/r(1)*(1-2*M/r(1))*dt(1)*dt(1)+M/r(1)/r(1)/(1-2*M/r(1))*dr(1)*dr(1)+r(1)*(1-2*M/r(1))*dfi(1)*dfi(1);
ddfi(1)=-2/r(1)*dr(1)*dfi(1);
for s=1:25000
    if r(s)<6*M
        risk=2;
    if r(s)<2*M
        risk=1;
        break;
    end
    end
    x(s)=r(s).*cos(fi(s));
    y(s)=r(s).*sin(fi(s));
    A=-(1-2*M/r(s))*dt(s);
    B=dr(s)*dr(s)/(1-2*M/r(s))+r(s)*r(s)*dfi(s)*dfi(s);
    C=-(1-2*M/r(s));
    kn=-A/B*sqrt(B/(A*A+B*C));
    atidal(s)=-(-2*M/r(s)/r(s)/r(s)*dt(s)*dt(s)*dr(s)*dr(s)-M/r(s)*C*dt(s)*dt(s)*dfi(s)*dfi(s)+M/r(s)/C*dr(s)*dr(s)*dfi(s)*dfi(s))*kn*kn;
    t21=t(s)+dt(s)*0.5*h;
    r21=r(s)+dr(s)*0.5*h;
    fi21=fi(s)+dfi(s)*0.5*h;
    dt21=dt(s)+ddt(s)*0.5*h;
    dr21=dr(s)+ddr(s)*0.5*h;
    dfi21=dfi(s)+ddfi(s)*0.5*h;
    ddt21=-2*M/r21/r21/(1-2*M/r21)*dt21*dr21;
    ddr21=-M/r21/r21*(1-2*M/r21)*dt21*dt21+M/r21/r21/(1-2*M/r21)*dr21*dr21+r21*(1-2*M/r21)*dfi21*dfi21;
    ddfi21=-2/r21*dr21*dfi21;
    t22=t21+dt21*0.5*h;
    r22=r21+dr21*0.5*h;
    fi22=fi21+dfi21*0.5*h;
    dt22=dt21+ddt21*0.5*h;
    dr22=dr21+ddr21*0.5*h;
    dfi22=dfi21+ddfi21*0.5*h;
    ddt22=-2*M/r22/r22/(1-2*M/r22)*dt22*dr22;
    ddr22=-M/r22/r22*(1-2*M/r22)*dt22*dt22+M/r22/r22/(1-2*M/r22)*dr22*dr22+r22*(1-2*M/r22)*dfi22*dfi22;
    ddfi22=-2/r22*dr22*dfi22;
    
    t61=t(s)+dt(s)*h/6;
    r61=r(s)+dr(s)*h/6;
    fi61=fi(s)+dfi(s)*h/6;
    dt61=dt(s)+ddt(s)*h/6;
    dr61=dr(s)+ddr(s)*h/6;
    dfi61=dfi(s)+ddfi(s)*h/6;
    ddt61=-2*M/r61/r61/(1-2*M/r61)*dt61*dr61;
    ddr61=-M/r61/r61*(1-2*M/r61)*dt61*dt61+M/r61/r61/(1-2*M/r61)*dr61*dr61+r61*(1-2*M/r61)*dfi61*dfi61;
    ddfi61=-2/r61*dr61*dfi61;
    t62=t61+dt61*h/6;
    r62=r61+dr61*h/6;
    fi62=fi61+dfi61*h/6;
    dt62=dt61+ddt61*h/6;
    dr62=dr61+ddr61*h/6;
    dfi62=dfi61+ddfi61*h/6;
    ddt62=-2*M/r62/r62/(1-2*M/r62)*dt62*dr62;
    ddr62=-M/r62/r62*(1-2*M/r62)*dt62*dt62+M/r62/r62/(1-2*M/r62)*dr62*dr62+r62*(1-2*M/r62)*dfi62*dfi62;
    ddfi62=-2/r62*dr62*dfi62;
    
    t31=t(s)+dt(s)*h/3;
    r31=r(s)+dr(s)*h/3;
    fi31=fi(s)+dfi(s)*h/3;
    dt31=dt(s)+ddt(s)*h/3;
    dr31=dr(s)+ddr(s)*h/3;
    dfi31=dfi(s)+ddfi(s)*h/3;
    ddt31=-2*M/r31/r31/(1-2*M/r31)*dt31*dr31;
    ddr31=-M/r31/r31*(1-2*M/r31)*dt31*dt31+M/r31/r31/(1-2*M/r31)*dr31*dr31+r31*(1-2*M/r31)*dfi31*dfi31;
    ddfi31=-2/r31*dr31*dfi31;
    
    t(s+1)=t22+9*(t62-t31);
    r(s+1)=r22+9*(r62-r31);
    fi(s+1)=fi22+9*(fi62-fi31);
    dt(s+1)=dt22+9*(dt62-dt31);
    dr(s+1)=dr22+9*(dr62-dr31);
    dfi(s+1)=dfi22+9*(dfi62-dfi31);
    ddt(s+1)=ddt22+9*(ddt62-ddt31);
    ddr(s+1)=ddr22+9*(ddr62-ddr31);
    ddfi(s+1)=ddfi22+9*(ddfi62-ddfi31);
    tp(s+1)=tp(s)+h;
end
sg=length(x)-45;
figure;set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.72]);
subplot(121);hold on;grid on;
rectangle('Position',[-2*M,-2*M,4*M,4*M],'Curvature',[1,1])
axis equal;
if risk==2
    xlabel({['Trajectory'];['(risks of unstable orbit)']},'Color','m');
elseif risk==1
    xlabel({['Trajectory'];['(risks of crossing event horizon)']},'Color','r');
else
    xlabel('Trajectory');
end
plot(0,0,'k*');
temp=plot(x(1),y(1),'o','LineWidth',2,'Color',clr(2+sign(dr(1)),:));
subplot(122);hold on;grid on
title('Tangential tidal effect');
xlabel('proper time(s)');
ylabel('tidal force(m/s^2)');
for s=0:30:sg
delete(temp);
subplot(121);
CLR=clr(2+sign(dr(s+15)),:);
plot(x(s+1:s+30),y(s+1:s+30),'k.','markersize',1)
temp=plot(x(s+30),y(s+30),'o','LineWidth',2,'Color',CLR);
title({['coordinate time t=' num2str(t(s+31)) 's'];  ['proper time tp=' num2str(tp(s+31)) 's']},'Color',CLR)
subplot(122);
plot(tp(s+1:s+30),atidal(s+1:s+30),'k.');
pause(0.02);
end