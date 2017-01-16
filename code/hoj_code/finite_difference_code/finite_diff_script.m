%  A resolution of Nx=12, Nt=20 seems to work.
x_min = -1.;
x_max = 1.;
Nx = 4*4;
t_min = 0.;
t_max = 1.;
Nt = 4*6;

[partial_t,partial_x] = get_pds(t_min,t_max,Nt,x_min,x_max,Nx);
L = partial_t + 0.5*partial_x;

%Boundary condition: f(t,1) = 0 for all time
delta_x_max = sparse(1,Nx);
delta_x_max(Nx) = 1.0;
Id_t = speye(Nt);
A_eq = kron(Id_t , delta_x_max );
b_eq = zeros(Nt,1);

%Boundary condition: f(t,0) = 0 for all time
%delta_x_min = sparse(1,Nx);
%delta_x_min(1) = 1.0;
%Id_t = speye(Nt);
%A_eq = cat(1,A_eq,kron(Id_t , delta_x_min ));
%b_eq = cat(1,b_eq,zeros(Nt,1));

%Constraint: L[f] <= 0
A_ub = L;
b_ub = zeros((Nx-1)*(Nt-1),1);

%Constraint: f(T,x) >= indicator(X_T) for all x
delta_T = sparse(1,Nt);
delta_T(Nt) = 1.;
A_ub = cat(1,A_ub, -kron( delta_T , speye(Nx) ) );
indicator_X_T = zeros(Nx,1);
indicator_X_T( (Nx/4+1):(3*Nx/4) ) = ones( Nx/2 ,1 );
b_ub = cat(1, b_ub , -indicator_X_T );

%Constraint: \partial_x f <= max_der
max_der = 1.0;
A_ub = cat(1,A_ub, partial_x );
b_ub = cat(1,b_ub, max_der*ones((Nx-1)*(Nt-1) , 1) );


%cost function is c[f] = int_X f(0,x) dx
delta_0 = sparse(1,Nt);
delta_0(1) = 1;
dx = 1.0 / (Nx-1);
int_X = dx*ones(1,Nx);
int_X(1) = dx*0.5;
int_X(Nx) = dx*0.5;
C = kron( delta_0 , int_X );

[f,fval,exit_flag] = linprog(C,A_ub,b_ub, A_eq, b_eq );
if(exit_flag==1)
    f = reshape(f,Nx,Nt);
    imagesc(linspace(t_min,t_max,Nt),linspace(x_min,x_max,Nx),f)
    xlabel('time');
    ylabel('space');
    title('f(t,x)');
end