#pragma once

#define V(i,j,k) (i*m_max_j*m_max_k + j*m_max_k + k) 
#define V_x(i,j,k) ((i+1)*m_max_j*m_max_k + j*m_max_k + k)
#define V_y(i,j,k) (i*m_max_j*m_max_k + (j+1)*m_max_k + k)
#define V_z(i,j,k) (i*m_max_j*m_max_k + j*m_max_k + (k+1)) 
#define alpha(phi_a,phi_b) ((phi_a)/(-phi_b+phi_a))
#define x_s(alpha,x_a,x_b) ((1.0-alpha)*x_a + alpha*x_b)

#define uint int_fast64_t // uint is an alias for unsigned int

