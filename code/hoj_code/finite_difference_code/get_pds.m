function [partial_t,partial_x] = get_pds( t_min, t_max, Nt ,x_min, x_max, Nx )
partial_x = kron( midpoint_op(Nt) , derivative_op(x_min,x_max,Nx) );
partial_t = kron( derivative_op(t_min,t_max,Nt), midpoint_op(Nx) );
end


