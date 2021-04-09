clear variables
close all
clc
%%

syms th(t) a(t)

th_dot = diff(th(t),t); a_dot = diff(a(t),t);

%%
syms m_p l_p r J_arm J_p g
%%
L = (1/2)*J_arm*(th_dot)^2 + (1/2)*J_p*(a_dot)^2 + ...
    (1/2)*m_p*(-(cos(th)*sin(a)*th_dot*l_p)-(sin(th)*cos(a)*a_dot*l_p)-(sin(th)*th_dot*r))^2 + ...
    (1/2)*m_p*(-(sin(th)*sin(a)*th_dot*l_p)-(cos(th)*cos(a)*a_dot*l_p)-(cos(th)*th_dot*r))^2 + ...
    (1/2)*m_p*((sin(a))^2)*(a_dot^2)*(l_p^2) + m_p*cos(a)*g*l_p;

%%

dL_dth = diff(L, th);
dL_dthdot = diff(L, th_dot);
dL_da = diff(L, a);
dL_dadot = diff(L, a_dot);

