clear variables
close all
clc
%%

syms th(t) a(t)

th_dot = diff(th(t),t); a_dot = diff(a(t),t);

%%
syms m_p l_p r J_arm J_p g
%%
T = (1/2)*J_arm*(th_dot)^2 + (1/2)*J_p*(a_dot)^2 + ...
    (1/2)*m_p*(-(cos(th)*sin(a)*th_dot*l_p)-(sin(th)*cos(a)*a_dot*l_p)-(sin(th)*th_dot*r))^2 + ...
    (1/2)*m_p*(-(sin(th)*sin(a)*th_dot*l_p)-(cos(th)*cos(a)*a_dot*l_p)-(cos(th)*th_dot*r))^2 + ...
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

first_eq = dL_dthdot_tDot-dL_dth;
secnd_eq = dL_dadot_tDot-dL_da;


down_f_eq = subs(first_eq, [th, a], [0, 0]);

% why is the lagrangian equation different from the derived model?

%% Linearization

syms r g mp lp Jp Jarm x1 x2 x3 x4

x = [x1; x2; x3; x4];


D = [r^2*mp+lp^2*mp-lp^2*cos(x2)^2*mp+Jarm, r*cos(x2)*mp*lp ;
    r*cos(x2)*mp*lp                       , lp^2*mp+Jp      ];

D_inv = simplify([D(2,2), -D(1,2); -D(2,1), D(1,1)]/det(D));



G = [0; mp*g*sin(x2)*lp];

C = [2*mp*cos(x2)*x4*lp^2*sin(x2), -mp*sin(x2)*x4*lp*r;
    -mp*cos(x2)*x3*lp^2*sin(x2) ,           0        ];


syms u Barm Bp Eg Kg Em Kt Km Rm 

tau_m = Eg*Kg*Em*Kt*(u-Kg*Km*x3)/Rm;

tau = [tau_m - Barm*x3;   -Bp*x4];



q_ddot = D_inv*tau - D_inv*G - D_inv*C*[x3; x4];

de_x = jacobian(q_ddot, x);

num_de_x = subs(de_x, [x1, x2], [0, 0])


de_u = jacobian(q_ddot, u);

num_de_u = subs(de_u, [x1, x2], [0, 0])

% maybe is it necessary to assume also x3,x4 = 0 ?
%%

down_de_x = simplify(subs(num_de_x, [x3, x4], [0, 0]))


%% 





