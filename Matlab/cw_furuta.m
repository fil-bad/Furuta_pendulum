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

%% Controllability

R_down = ctrb(A_down,B_down);
                
R_up= ctrb(A_up,B_up);
      
%% Pole Placement (down)

sys_down = ss(A_down, B_down, eye(4), zeros(4,1));
damp(sys_down)

K_down = place(A_down, B_down, [-3, -18, -8, -10]);

sys_d_cl = ss(A_down-B_down*K_down, B_down, eye(4), zeros(4,1));
damp(sys_d_cl)

% figure(3)
% rlocus(sys_d_cl)

%% Pole Placement (up)

sys_up = ss(A_up, B_up, eye(4), zeros(4,1));
damp(sys_up)

K_up= place(A_up, B_up, [-3, -23, -4, -6]);


                
% sys_u_cl = ss(A_up-B_up*K_up, B_up, [1 1 1 1], 0);
sys_u_cl = ss(A_up-B_up*K_up, B_up, eye(4), zeros(4,1));
damp(sys_u_cl)

% figure(4)
% rlocus(sys_u_cl)
                
%% LQR (down)

Q_d = [0.01, zeros(1,3);
       zeros(3,1), eye(3)];
   
R_d = 10;

[K_lqr_d, P_d, CLP_d] = lqr(A_down, B_down, Q_d, R_d);

sys_lqr_d_cl = ss(A_up-B_up*K_lqr_d, B_up, eye(4), zeros(4,1));


%% LQR (up)

Q_u = Q_d;
R_u = 100;

[K_lqr_u, P_u, CLP_u] = lqr(A_up, B_up, Q_u, R_u);

sys_lqr_u_cl = ss(A_up-B_up*K_lqr_u, B_up, eye(4), zeros(4,1));

save("matrices.mat")

%% Simulations

% x0 =   [pi-0.1;      % theta
%         pi/8;     % alpha
%         0.3;      % theta_dot
%         0.1];     % alpha_dot
% 

x0 =   [0;      % theta
        0;     % alpha
        0;      % theta_dot
        0];     % alpha_dot
    
    
    
%% Plot Whole Graph

leg = legend('$\theta(t)$','$\alpha(t)$',...
                '$\dot{\theta}(t)$','$\dot{\alpha}(t)$',...
            ... %    '$\alpha_{ref}(t)$',...
                'Interpreter', 'latex');
set(leg, 'Interpreter', 'latex');

%% USE THIS FOR 2x2 LAYOUT OF SCOPE
set(gcf, 'InvertHardCopy', 'off')
saveas(gcf, '../Lyx/imgs/NL_sw_K_up_f_40a.svg','svg')


%% Plot Single evolution

leg = legend('$\alpha_{cl}(t)$', 'Interpreter', 'latex');
set(leg, 'Interpreter', 'latex');

set(gcf, 'InvertHardCopy', 'off')
saveas(gcf, '../Lyx/imgs/NL_K_down_alpha.svg','svg')






%%




%% Nonlinear evolution (ODE)

% Check stability
steps = 15; % USE AT LEAST TWO.
state_blue = [];
state_red = [];

for alpha = -pi:(2*pi/(steps-1)):pi
    for theta = -pi:(2*pi/((steps-1))):pi
        
        x_init = [theta; alpha; 0; 0];
        
        tstart = tic;
        [~,x] = ode45(@(t,x) furuta_ode(t,x,K_down),[0,10],x_init,...
                        odeset('Events',@(t,x) time_event(t,x,tstart)));
        
        error = abs( x(end,2)/(2*pi) - round(x(end,2)/(2*pi)) );
        if (error < 0.1)
            state_blue(end+1,:) = [theta, alpha];
            
        else
            state_red(end+1,:) = [theta, alpha];
        end
        theta
    end
    alpha
end

% plot graph

stab_fig = figure(2);
clf(stab_fig);
grid on
hold on

title("Stability for downward position ( $K_{lqr/pp},\ |\theta_{max}(0)| = |\alpha_{max}(0)|= \pi$)",...
    "Interpreter","latex");

xlabel('\theta')
ylabel('\alpha')
% xlim([-pi pi])
% ylim([-pi pi])

if (height(state_blue) > 0)    
    scatter(state_blue(:,1), state_blue(:,2), 'blue', 'filled');
end

if (height(state_red) > 0)    
    scatter(state_red(:,1), state_red(:,2), 'red', 'filled');
end

saveas(stab_fig, '../Lyx/imgs/Stab_K_d.svg','svg')

%%

title("Stability for downward position ( $K_{lqr},\  |\theta_{max}(0)|=\pi, |\alpha_{max}(0)|= \frac{\pi}{2}$)",...
    "Interpreter","latex");


saveas(stab_fig, '../Lyx/imgs/Stab_K_lqr_u2.svg','svg')


%% Finer interval

state_blue = [];
state_red = [];

% for alpha = 0.67:0.005:0.8 % K_up
for alpha = 0.78:0.005:0.90 % K_lqr_u
        
        theta = 0;
        
        x_init = [theta; alpha; 0; 0];
        
        tstart = tic;
        [~,x] = ode45(@(t,x) furuta_ode(t,x,K_lqr_u),[0,30],x_init,...
                        odeset('Events',@(t,x) time_event(t,x,tstart)));
        
        error = abs( x(end,2)/(2*pi) - round(x(end,2)/(2*pi)) );
        if (error < 0.1)
            state_blue(end+1,:) = [theta, alpha];
            
        else
            state_red(end+1,:) = [theta, alpha];
        end
    alpha
end



fig = figure(3);
clf(fig);
grid on
hold on

title("Finer stability, $K_{lqr}$","Interpreter","latex");

xlabel('\theta')
ylabel('\alpha')

if (height(state_blue) > 0)    
    scatter(state_blue(:,1), state_blue(:,2), 'blue', 'filled');
end

if (height(state_red) > 0)    
    scatter(state_red(:,1), state_red(:,2), 'red', 'filled');
end

saveas(fig, '../Lyx/imgs/Finer_Stab_K_lqr_u.svg','svg')








