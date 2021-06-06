function x_dot = furuta_ode(t,x,K)
% Furuta pendulum 

% renaming state in compliance with altready written code

th = x(1); 
a = x(2); % Downward position
% a = x(2)+pi; % Upward position
th_d = x(3);
a_d = x(4);

% reference and control signal


x_ref = zeros(4,1); % origin as set point

% x_ref(2) = pi; % for Upward position

u = K*(x_ref - [th; a; th_d; a_d]);


% assigning values to symbolic values 

Rm = 2.6; Kt = 7.68e-3; Em = 0.69; Km = 7.68e-3; Kg = 70;
Eg = 0.9; mp = 0.127; lp = 0.1556; Jp = 0.0012; Jarm = 0.002;
Bp = 0.0024; Barm = 0.0024; r = 0.2159; g = 9.81;

% defining matrices

D = [(r^2)*mp+(lp^2)*mp-(lp^2)*(cos(a)^2)*mp+Jarm, r*cos(a)*mp*lp ;
    r*cos(a)*mp*lp                              , (lp^2)*mp+Jp      ];

C = [2*mp*cos(a)*a_d*(lp^2)*sin(a), -mp*sin(a)*a_d*lp*r;
    -mp*cos(a)*th_d*(lp^2)*sin(a) ,           0        ];

G = [     0; 
    mp*g*sin(a)*lp];

tau_m = Eg*Kg*Em*Kt*(u-Kg*Km*th_d)/Rm;

tau = [tau_m - Barm*th_d;
            -Bp*a_d     ];

% defining the non-linear equations
x_dot  = zeros(4,1);

x_dot(1:2) = [th_d;
               a_d]; 

x_dot(3:4) = (D^-1)*(tau - G - C*[th_d; a_d]);

end

