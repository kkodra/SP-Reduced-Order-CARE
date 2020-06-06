function results(output)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
set(0,'DefaultTextInterpreter','Latex',...
    'DefaultAxesFontSize',18,'defaultAxesTickLabelInterpreter','latex',...
    'defaultLineLineWidth',2.5);

P     = output.P;
P_i   = output.P_i;
Jdiff = output.Jdiff;
Pdiff = output.Pdiff;

norm_P_Pi = norm(P - P_i);
fprintf(2,'Norm of P and P_i from last iteration:\n')
disp(repmat('=',1,35))
fprintf('||P - P_i|| = %.15e\n',norm_P_Pi)

fprintf(2,'\n\nDifference of truth and obtained cost functions:\n')
disp(repmat('=',1,50))
diff_J = Jdiff

fprintf(2,'\n\nNorm of P and P_i for every iteration:\n')
disp(repmat('=',1,40))
norm_P = Pdiff

figure
plot(Jdiff,'*-'); grid on; xlabel('Iteration');
ylabel('$|J_{truth} - J_{i}|$');
title('Cost Function Difference')

figure
plot(Pdiff,'*-');grid on;xlabel('Iteration');
ylabel('$||P_{truth} - P_{i}||$')
title('Norm of CARE Solutions')
end

