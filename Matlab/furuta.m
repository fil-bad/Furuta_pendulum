clear variables
close all
clc
%%

syms t th(t) a(t)

th_dot = diff(th(t),t); a_dot = diff(a(t),t);

%%
syms m_p l_p r J_arm J_p g
%%
T = (1/2)*J_arm*(th_dot)^2 + (1/2)*J_p*(a_dot)^2 + ...
    (1/2)*m_p*(-cos(th)*sin(a)*th_dot*l_p-(sin(th)*cos(a)*a_dot*l_p)-(sin(th)*th_dot*r))^2 + ...
    (1/2)*m_p*(-(sin(th)*sin(a)*th_dot*l_p)+(cos(th)*cos(a)*a_dot*l_p)+(cos(th)*th_dot*r))^2 + ...
    (1/2)*m_p*((sin(a))^2)*(a_dot^2)*(l_p^2);

V = -m_p*cos(a)*g*l_p;

L = T-V;

%%

dL_dth = diff(L, th);
dL_dthdot = diff(L, th_dot);

dL_dthdot_tDot = diff(dL_dthdot, t);

dL_da = diff(L, a);
dL_dadot = diff(L, a_dot);

dL_dadot_tDot = diff(dL_dadot, t);

%%

first_eq = dL_dthdot_tDot - dL_dth;
secnd_eq = dL_dadot_tDot - dL_da;


% % down_f_eq = subs(first_eq, [th, a], [0, 0]);

% why is the lagrangian equation different from the derived model?

%% 
clear variables
close all
clc

%% Linearization

syms r g mp lp Jp Jarm x1 x2 x3 x4 th a th_d a_d

x = [th; a; th_d; a_d];


D = [(r^2)*mp+(lp^2)*mp-(lp^2)*(cos(a)^2)*mp+Jarm, r*cos(a)*mp*lp ;
    r*cos(a)*mp*lp                              , (lp^2)*mp+Jp      ];

D_inv = simplify(D^-1);

C = [2*mp*cos(a)*a_d*(lp^2)*sin(a), -mp*sin(a)*a_d*lp*r;
    -mp*cos(a)*th_d*(lp^2)*sin(a) ,           0        ];

G = [0; mp*g*sin(a)*lp];



syms u Barm Bp Eg Kg Em Kt Km Rm 

tau_m = Eg*Kg*Em*Kt*(u-Kg*Km*th_d)/Rm;

tau = [tau_m - Barm*th_d;   -Bp*a_d];


q_ddot = D_inv*(tau - G - C*[th_d; a_d]);


de_x = jacobian(q_ddot, x);
de_u = jacobian(q_ddot, u);


down_de_x = subs(de_x, [th, a, th_d, a_d], [0, 0, 0, 0]);
down_de_u = subs(de_u, [th, a, th_d, a_d], [0, 0, 0, 0]);

up_de_x = subs(de_x, [th, a, th_d, a_d], [0, pi, 0, 0]);
up_de_u = subs(de_u, [th, a, th_d, a_d], [0, pi, 0, 0]);


%%

num_down_de_x = subs(down_de_x,...
                    [Rm Kt Em Km Kg Eg mp lp Jp Jarm Bp Barm r g],...
                    [2.6, 7.68e-3, 0.69, 7.68e-3, 70, 0.9, 0.127, 0.1556,...
                    0.0012, 0.002, 0.0024, 0.0024, 0.2159, 9.81]);

num_down_de_u = subs(down_de_u,...
                    [Rm Kt Em Km Kg Eg mp lp Jp Jarm Bp Barm r g],...
                    [2.6, 7.68e-3, 0.69, 7.68e-3, 70, 0.9, 0.127, 0.1556,...
                    0.0012, 0.002, 0.0024, 0.0024, 0.2159, 9.81]);

num_up_de_x = subs(up_de_x,...
                    [Rm Kt Em Km Kg Eg mp lp Jp Jarm Bp Barm r g],...
                    [2.6, 7.68e-3, 0.69, 7.68e-3, 70, 0.9, 0.127, 0.1556,...
                    0.0012, 0.002, 0.0024, 0.0024, 0.2159, 9.81]);

num_up_de_u = subs(up_de_u,...
                    [Rm Kt Em Km Kg Eg mp lp Jp Jarm Bp Barm r g],...
                    [2.6, 7.68e-3, 0.69, 7.68e-3, 70, 0.9, 0.127, 0.1556,...
                    0.0012, 0.002, 0.0024, 0.0024, 0.2159, 9.81]);


A_down = double([0,0,1,0;
                 0,0,0,1;
                 num_down_de_x]);        
                
A_up = double([0,0,1,0;
               0,0,0,1;
               num_up_de_x]);
           
B_down = double([0;
                 0;
                num_down_de_u]);

B_up = double([0;
               0;
               num_up_de_u]);
                              
                
save("lin_systems.mat", 'A_down', 'A_up','B_down','B_up');        

%%
clear variables
close all
clc
                
%% System analisys

load("lin_systems.mat")
                
% Stability

eig_down = eig(A_down);

eig_up = eig(A_up);

%%

R_down = ctrb(A_down,B_down);
                
R_up= ctrb(A_up,B_up);
      
                
                
                
                
                
                
                
                
                
                


