% Rewrite of simulation for Automatica. A simple example used to test
% algorithm.

clear; clc
format long g

A = [0 1 0 0; 0 0 1 0; 0 0 0 100; -150 -505 -470 -125];
eps = 0.1;
B = [1 0; 2 0; 0 100; 0 100];

Q = [1 0 0.01 0.02; 0 2 0.01 0.03; 0.01 0.01 0.02 0; 0.02 0.03 0 0.02];

R = 1;

[nB,mB] = size(B);

% Reduce to order 2
red = 2;

A1 = A(1:red,1:red);
A2 = A(1:red,red+1:end);
A3 = eps*A(red+1:end,1:red);
A4 = eps*A(red+1:end,red+1:end);

redB = mB/2;

B1 = B(1:red,1:redB);
B2 = B(1:red,redB+1:end);
B3 = eps*B(red+1:end,1:redB);
B4 = eps*B(red+1:end,redB+1:end);

Q1 = Q(1:red,1:red);
Q2 = Q(1:red,red+1:end)/eps;
Q3 = Q(red+1:end,red+1:end)/eps;

Q = [Q1 eps*Q2; eps*Q2' eps*Q3];
A0 = A1 - A2/(A4)*A3;
S1 = B1/(R)*B1'; S4 = B4/(R)*B4';
S = [S1 zeros(size(A2)); zeros(size(A3)) S4/eps^2];

% Evaluate zero-order quantities P^(0)_i, i = 1,2,3
BRBT = S1 + A2/(A4)*S4/(A4)'*A2';
P30 = zeros(size(A4));
P10 = are(A0,BRBT,Q1);
P20 = -P10*A2/(A4);

% Value of original J (truth)
P = care(A,B,Q);
Aorig_tot = A - S*P;
Qorig_tot = Q + P'*S*P;
Vorig_tot = lyap(Aorig_tot',Qorig_tot);  
Jorig_tot = trace(Vorig_tot);

%% Algorithm implementation to evaluate accurate P.
% Select number of iterations (for this example converges in 10 iterations
iter = 12;
Jdiff = zeros(1,iter);
Pdiff = zeros(1,iter);

D1 = A1 - S1*P10;
D2 =  A3 - S4*P20';
As = D1 - A2/(A4)*D2;

% Initial conditions evaluated below.
E30 = lyap(A4', Q3 +P20'*A2 + A2'*P20);

H20 = P10*S1*P20 - A1'*P20;

E10 = lyap(As',(D2'*E30 + Q2 - H20)/(A4)*D2 + D2'/(A4)'*(D2'*E30 +Q2 - H20)');

E20 = -(E10*A2 + D2'*E30 +Q2 )/(A4) + H20/(A4);

P0 = [P10 eps*P20; eps*P20' eps*P30];

Ai_tot = A - S*P0;
Qi_tot = Q + P0'*S*P0;
Vi_tot = lyap(Ai_tot', Qi_tot);
Ji_tot = trace(Vi_tot);
for i = 1:iter
    H40i = P20'*S1*P20 + eps*P20'*S1*E20 + eps*E20'*S1*P20 + eps^2*E20'*S1*E20 +E30*S4*E30;
    
    % Evaluate E3 & update P3_i
    Q_Tot3 = Q3 + P20'*A2 + A2'*P20 +eps*A2'*E20 + eps*E20'*A2 - eps*H40i;
    E3i = lyap(A4',Q_Tot3);
   
    % Check if E1 is solved correctly
    E3_Check = A4'*E3i + E3i*A4 + Q_Tot3;
    
    P3i = P30 + eps*E3i;

    % Evaluate E1 & update P1_i     
    H20i = P10*S1*P20 - A1'*(P20 + eps*E20);
    H30i = P10*S1*E20 + E10*S1*P20 + eps*E10*S1*E20 + E20*S4*E3i;  % 
    H10 = E10*S1*E10 + E20*S4*E20';
    Q_E1i = (D2'*E3i + Q2 - H20i - eps*H30i)/(A4)*D2 + D2'/(A4)'*(D2'*E3i + Q2 - H20i - eps*H30i)' + eps*H10;  
    E1i = lyap(As',-Q_E1i);

    % Check if E1 is solved correctly
    E1_Check = D1'*E1i + E1i*D1 + D2'*E20' + E20*D2 - eps*H10;
    
    P1i = P10 + eps*E1i;
    
    % Evaluate E2 & update P2_i
    E2i = -(E1i*A2 + D2'*E3i +Q2 )/(A4) + (H20i + eps*H30i)/(A4);

    % Check if E2 is solved correctly
    E2_Check = E2i*A4 + E1i*A2 + D2'*E3i +Q2 - (H20i + eps*H30i);

    P2i = P20 + eps*E2i;
    
    E30 = E3i; E20 = E2i; E10 = E1i;
    
    % Form P.
    P_i = [P1i eps*P2i; eps*P2i' eps*P3i];
    E = [E1i E2i; E2i' E3i];
    Ai_tot = A - S*P_i;
    Qi_tot = Q + P_i'*S*P_i;
    Vi_tot = lyap(Ai_tot', Qi_tot);
    Ji_tot = trace(Vi_tot);
    Jdiff(i) = -(Jorig_tot - Ji_tot);
    Pdiff(i) = norm(P - P_i);
end

P_i = [P10+eps*E1i eps*(P20+eps*E2i); eps*(P20+eps*E2i)' eps*(P30+eps*E3i)];

output.P = P;
output.P_i = P_i;
output.Jdiff = Jdiff';
output.Pdiff = Pdiff';

results(output)