%Comparison of energy-preserving and generic sequential Operator Inference
close all
clear all
%% Comparative error plot
errE=importdata('errBurg_energy.txt');
errE=errE.data;
errS=importdata('errBurg_simple.txt');
errS=errS.data;
errP=importdata('err_reproj.txt');
errP=errP.data;


errEm=mean(errE,2);
errSm=mean(errS,2);
errPm=mean(errP,2);

rr=[5,10,15,20,25,30]';
figure(1)
plot(rr,100*errEm,'r*-')
hold on
plot(rr,100*errSm,'bo-')
plot(rr,100*errPm,'kx-')
legend({'\texttt{Seq\_OpInf\_EP} $\left( \mathbf{x}^T \mathbf{H} \left( \mathbf{x} \otimes \mathbf{x} \right) = 0 \right)$','\texttt{Seq\_OpInf}','Reprojected data'},'Interpreter','latex')
xlabel('$r$','Interpreter','latex')
ylabel('Average Error $\times \; 100 \%$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');

%% Energy contribution of each inferred term

terms_skew=importdata('term_energy_skew.txt');
terms_opinf=importdata('term_energy_simple.txt');

askew=terms_skew(:,1);
bskew=terms_skew(:,2);
time=terms_skew(:,3);

aop=terms_opinf(:,1);
bop=terms_opinf(:,2);

figure(2)
plot(time,bskew,'--r','LineWidth',3)
hold on
plot(time,askew,'--b','LineWidth',2)
plot(time,bop,'magenta','LineWidth',2)
plot(time,aop,'green','LineWidth',2)
plot(time,aop+bop,'k','LineWidth',2)
hold off
set(gca,'TickLabelInterpreter','latex');
xlabel('$t\;(s)$','Interpreter','latex')
ylabel('Energy Contribution','Interpreter','latex')
legend({'$\mathbf{x}^T \mathbf{A}_{sE} \mathbf{x} $','$\mathbf{x}^T \mathbf{H}_{sE} \left( \mathbf{x} \otimes \mathbf{x} \right) $','$\mathbf{x}^T \mathbf{A}_s \mathbf{x} $','$\mathbf{x}^T \mathbf{H}_s \left( \mathbf{x} \otimes \mathbf{x} \right)$','$\mathbf{x}^T \left( \mathbf{A}_s \mathbf{x} + \mathbf{H}_s \left( \mathbf{x} \otimes \mathbf{x} \right) \right)$'},'Interpreter','latex')
set(gca, 'FontSize', 18)