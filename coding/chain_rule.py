# filename: chain_rule.py

import sympy as sp

def chain_rule(f_of_u_expr, u_of_x_expr, x):
    # differentiating each function
    f_of_u_diff = sp.diff(f_of_u_expr, u)
    u_of_x_diff = sp.diff(u_of_x_expr, x)
    
    # substituting u=g(x) into f'(u)
    f_of_u_diff_subs = f_of_u_diff.subs(u, u_of_x_expr)

    # multiplying the derivatives according to the chain rule
    chain_rule_diff = f_of_u_diff_subs * u_of_x_diff

    return chain_rule_diff

x, u = sp.symbols('x u')

f_of_u = u**3  # You can replace with your function
u_of_x = sp.sin(x)  # You can replace with your function

result = chain_rule(f_of_u, u_of_x, x)

print(result)