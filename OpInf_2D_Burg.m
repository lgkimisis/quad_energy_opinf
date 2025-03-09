clear all
close all
run custom_cmap.m
%Import data from the simulations
Nt=401; %Number of time steps
Nti=1;
M1=Nt-3;
M=M1;
dt=.01;
xdata = importdata('.\xpos_nonl2D.txt');
xdata=xdata.data';
Nx=length(xdata);
ydata = importdata('.\ypos_nonl2D.txt');
ydata=ydata.data';
Ny=length(ydata);
[X,Y]=meshgrid(xdata,ydata);

udata = importdata('.\usol2_nonl2D.txt');
udata=udata.data;

udata=udata(:,Nti:Nt);
umean=mean(mean(udata));

NT=Nt-Nti+1;
%% Time derivatives numerical computation
ind  = 1:NT-1;
unew =(1/6*udata(:,1:end-3) - 1*udata(:,2:end-2) + 1/2*udata(:,3:end-1) + 1/3*udata(:,4:end))/dt;

%% SVD on the states (velocity)
tic
Udata = udata(:,1:end-3);

dUdata = unew;

dimvec=[5,10,15,20,25,30]; 
ndim=length(dimvec);

endtime=M1;

errU_avg=zeros(ndim,endtime);
errU_max=zeros(ndim,endtime);
errU_rel=zeros(ndim,endtime);
errUmin=zeros(ndim,endtime);

for jj=1:ndim
    [U1,S1,V] = svd(Udata,'econ');
    
    r = dimvec(jj);  
    U1 = U1(:,1:r); 
    S1 = S1(1:r,1:r);
    
    %% Make the velocity ROM
    Uold=U1.'*Udata(:,1:M);
    dUold = U1.'*dUdata(:,1:M);  %outputs
    A2=dUold;
    
    %Initialize operators
    Atilde=zeros(r,r);
    Htilde=zeros(r,r*(r+1)/2);
    
    th=0:pi/100:pi/2;
    
    k_n=r; 
    c_n=logspace(-5,3,50);
    
    vkron=[];
    for kk=1:M %timesteps (corresponds to u(2:M))
        longa=[];
        %We write this to avoid repeating terms (quadratic)
        for j=1:r
            long1=Uold(j,kk)*Uold(j:r,kk);
            longa=[longa;long1];
        end
        vkron=[vkron, longa];
    end
    
    qv=size(vkron,1);
    %
    A1=[Uold; vkron]; %removed Uin
    
    %We regularize here
    AA=A1*A1.';
    
    for ii=1:r
        
        b_res=zeros(length(c_n),1);
        x_norm=zeros(length(c_n),1);
        for f=1:length(k_n)
            for j=1:length(c_n)
                Reg=eye(r+qv,r+qv);
                Reg(1:r,1:r)=c_n(j)/k_n(f)*Reg(1:r,1:r);
                Reg(r+1:end,r+1:end)=c_n(j)*Reg(r+1:end,r+1:end);
                
                G1=(AA+Reg)\A1*A2(ii,:)'; %Transposing the eqn.
                
                b_res(j,f)=norm(A2(ii,:)-G1.'*A1,'fro'); %Frobenius norm
                
                x_norm(j,f)=norm(G1,'fro');
                
            end
        end
        
        b_resN=(b_res-min(min(b_res)))/(max(max(b_res))-min(min(b_res)));
        x_normN=(x_norm-min(min(x_norm)))/(max(max(x_norm))-min(min(x_norm)));
        
        rmin=sqrt(min(min(x_normN.^2+b_resN.^2)));
        [i_b,j_b]=find(x_normN.^2+b_resN.^2==min(min(x_normN.^2+b_resN.^2)));
        
        c=c_n(i_b); %will be kept for the next LS problem
        kn=k_n(j_b);  
        
        Reg=eye(r+qv,r+qv);
        Reg(1:r,1:r)=c/kn*Reg(1:r,1:r);
        Reg(r+1:end,r+1:end)=c*Reg(r+1:end,r+1:end);
        
        AA=A1*A1';
        
        G1=(AA+Reg)\A1*A2(ii,:)';
        
        Atilde(ii,:)=G1(1:r)'; 
        Htilde(ii,:)=G1(r+1:end)';
        
    end
    
    toc
    %% Simulation
    disp('ROM simulation')
    h1 = figure('Renderer', 'painters', 'Position', [10 10 800 260]);
    axis tight manual % this ensures that getframe() returns a consistent size
    u_old=Uold(:,1);
    tstp=endtime;
    
    dt_s=dt;
    for t=2:endtime
        
        %RK4 method
        %Step 1
        u_kron=[];
        for i=1:r
            u_kron=[u_kron; u_old(i)*u_old(i:r)];
        end
        
        du_tilde1=Atilde*u_old+Htilde*u_kron;
        u_old1=du_tilde1*dt_s/2+u_old;
        
        %Step 2
        u_kron1=[];
        for i=1:r
            u_kron1=[u_kron1; u_old1(i)*u_old1(i:r)];
        end
        
        du_tilde2=Atilde*u_old1+Htilde*u_kron1;
        u_old2=du_tilde2*dt_s/2+u_old;
        
        %Step 3
        u_kron2=[];
        for i=1:r
            u_kron2=[u_kron2; u_old2(i)*u_old2(i:r)];
        end
        
        du_tilde3=Atilde*u_old2+Htilde*u_kron2;
        u_old3=du_tilde3*dt_s+u_old;
        
        %Step 4
        u_kron3=[];
        for i=1:r
            u_kron3=[u_kron3; u_old3(i)*u_old3(i:r)];
        end
        
        du_tilde4=Atilde*u_old3+Htilde*u_kron3;
        
        %Overall
        du_tilde=(du_tilde1+2*du_tilde2+2*du_tilde3+du_tilde4)/6;
        
        u_old=du_tilde*dt_s+u_old;
        
        u_new=U1*u_old;
        
        
        if mod(t,tstp)==0
            minclim=min(udata(:,t));
            maxclim=max(udata(:,t));
            fig=figure(h1);
            set(gca, 'FontSize', 32)
            set(gca, 'FontName', 'Times New Roman')
            
            a1=subplot(1,2,1,'Parent',fig);
            contourf(a1,X,Y,reshape(udata(:,t),Nx,Ny),50,'LineColor','none')
            colormap(map)
            caxis(a1,[minclim,maxclim])
            xlabel('$x$','Interpreter','latex')
            ylabel('$y$','Interpreter','latex')
            title('Numerical simulation','Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex');

            a2= subplot(1,2,2,'Parent',fig);
          
            contourf(a2,X,Y,reshape(u_new,Nx,Ny),50,'LineColor','none')
            colormap(map)
            caxis(a2,[minclim,maxclim])
            xlabel('$x$','Interpreter','latex')
            ylabel('$y$','Interpreter','latex')
            title('\texttt{Seq\_OpInf}','Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex');
            
            set(a1,'Units','normalized', 'position', [0.1 0.2 0.32 0.6]);
            set(a2,'Units','normalized', 'position', [0.55 0.2 0.32 0.6]);
            
            h = axes(fig,'visible','off');
            c = colorbar(h,'Position',[0.9 0.168 0.022 0.7]);  % attach colorbar to h
            colormap(c,map)
            caxis(h,[minclim,maxclim]);
            set(c,'TickLabelInterpreter','latex')
            c.Title.String = 'u';
            c.Title.Interpreter = 'latex';
            drawnow
            
        end
        errU_avg(jj,t)=norm(udata(:,t)-u_new,2);
        errU_max(jj,t)=max(abs(udata(:,t)-u_new));
        errU_rel(jj,t)=mean(abs(udata(:,t)-u_new)./max(udata(:,t)));
        errUmin(jj,t)=mean(abs(udata(:,t)-U1*U1'*udata(:,t))./max(udata(:,t)));
        t
        
    end
end

%% Save error and terms energy vectors

Hfull=zeros(r,r^2);

for i=1:r
    zstart = (2*r-i)*(i-1)/2+i; %unique DOF counting (global)    
    Hfull(:,((i-1)*r+i):i*r)=Htilde(:,zstart:zstart+r-i);
end

dlmwrite('FullQuad_Mat.txt',Hfull)

T3 = table(errU_rel);
writetable(T3,'errBurg_simple.txt');

T4 = table(errUmin);
writetable(T4,'err_reproj.txt');

a=zeros(endtime,1);
for i=1:endtime
x=U1'*udata(:,i);
a(i)=x'*Hfull*kron(x,x);
end

b=zeros(endtime,1);
for i=1:endtime
x=U1'*udata(:,i);
b(i)=x'*Atilde*x;
end

dlmwrite('term_energy_simple.txt',[a, b, dt*(1:endtime)'])
