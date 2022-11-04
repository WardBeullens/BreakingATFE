# used to generate table 2: Solving ATFE problems with auxiliary information

import time
import random
from sage.libs.giac import groebner_basis as gb_giac

n = 10
q = 524287 # 524287, 131071, 65521
rank = 6

K = GF(q)
Kn = VectorSpace(K,n)
Kn_basis = Kn.basis()

# trilinear form stuff

form_basis = [ ]
for i in range(n):
    for j in range(i+1,n):
        for k in range(j+1,n):
            form_basis.append((i,j,k))

def evaluate_form(form,u,v,w):
    eval = 0
    for x,coef in enumerate(form):
        i,j,k = form_basis[x]
        d  = u[i]*v[j]*w[k]
        d += u[k]*v[i]*w[j]
        d += u[j]*v[k]*w[i]
        d -= u[i]*v[k]*w[j]
        d -= u[j]*v[i]*w[k]
        d -= u[k]*v[j]*w[i]
        eval += d*coef
    return eval

def form_to_mat(form,u):
    return matrix(K, [ [ evaluate_form(form,Kn_basis[i],Kn_basis[j],u) for j in range(n)] for i in range(n)] )

def act_on_form(form, S):
    new_form = []
    for i,j,k in form_basis:
        new_form.append(evaluate_form(form,S.column(i),S.column(j),S.column(k)))
    return new_form

def random_invertible(K,n):
    S = random_matrix(K,n,n)
    while not S.is_invertible():
        S = random_matrix(K,n,n)
    return S

def generate_form_with_low_rank_point(rank):
    zeroes = math.comb(n-1,2) - math.comb(rank,2)

    while True:
        form = [K(0)]*zeroes + [K.random_element() for _ in range(len(form_basis)-zeroes)]

        if form_to_mat(form, Kn_basis[0]).rank() == rank:
            break

    S = random_invertible(K,n)
    S_inv = S.inverse()

    form = act_on_form(form,S)
    point = S_inv*Kn_basis[0]

    return form, point


var_names = [ x for y in zip([ "s"+str(i) for i in range(n*n-1,-1,-1) ],[ "t"+str(i) for i in range(n*n-1,-1,-1) ]) for x in y]
Solver_PR = PolynomialRing(K, var_names)

S_vars = Matrix(Solver_PR,n,n,reversed(Solver_PR.gens()[::2])).transpose()
T_vars = Matrix(Solver_PR,n,n,reversed(Solver_PR.gens()[1::2])).transpose()

def find_equivalence(form1,u1,form2,u2):

    start_time = time.time()

    M1 = form_to_mat(form1,u1)
    M2 = form_to_mat(form2,u2)

    ker1 = M1.right_kernel().basis()
    ker2 = M2.right_kernel().basis()

    if len(ker1) != len(ker2):
        print('No equivalence exists!')
        return

    ker_dim = len(ker1)

    change_of_basis1 = [u1]
    for b in ker1+Kn_basis:
        if b not in span(change_of_basis1):
            change_of_basis1.append(b)

    change_of_basis2 = [u2]
    for b in ker2 + Kn_basis: 
        if b not in span(change_of_basis2):
            change_of_basis2.append(b)

    change_of_basis1 = Matrix(K,change_of_basis1).transpose()
    change_of_basis2 = Matrix(K,change_of_basis2).transpose()
    form1 = act_on_form(form1,change_of_basis1)
    form2 = act_on_form(form2,change_of_basis2)

    S = copy(S_vars)
    T = copy(T_vars)

    for i in range(1,n):
        S[i,0] = 0
        T[i,0] = 0

    for i in range(1,ker_dim):
        for j in range(ker_dim,n):
            S[j,i] = 0
            T[j,i] = 0

    # build system of equations
    equations = []

    ST = S*T
    TS = T*S

    for i in range(n):
        for j in range(n):
            if i==j:
                equations.append(ST[i][j]-1)
                equations.append(TS[i][j]-1)
            else:
                equations.append(ST[i][j])
                equations.append(TS[i][j])

    Form2 = act_on_form(form1,S)
    for a,b in zip(Form2,form2):
        equations.append(a-b)

    Form1 = act_on_form(form2,T)
    for a,b in zip(Form1,form1):
        equations.append(a-b)

    # f2(Tx,y,z) = f1(x,Sy,Sz)
    for i in range(n):
        for j in range(n):
            for k in range(j,n):
                equations.append(evaluate_form(form2, T.column(i), Kn_basis[j], Kn_basis[k]) - evaluate_form(form1, Kn_basis[i], S.column(j), S.column(k)))
                
    # f2(Tx,Ty,z) = f1(x,y,Sz)
    for i in range(n):
        for j in range(i,n):
            for k in range(n):
                equations.append(evaluate_form(form2, T.column(i), T.column(j), Kn_basis[k]) - evaluate_form(form1, Kn_basis[i], Kn_basis[j], S.column(k)))

    print('Start GB computation')
    I = ideal(equations)
    GB = gb_giac(equations) 

    ## read solution from groebner basis
    s_zero = (-GB[0].coefficients()[1]).nth_root(3)
    S[0,0] = s_zero
    for i in range(n):
        for j in range(n):
            if (i,j) == (0,0):
                continue

            if S[i,j] == 0:
                continue

            for g in GB:
                if g.lm() == S[i,j]:
                    S[i,j] = -g.coefficients()[1]*s_zero

    S = change_of_basis1*S*change_of_basis2.inverse()

    end = time.time()
    print("Time to find the equivalence:", end - start_time)

    return S

# generate first form
form1, point1 = generate_form_with_low_rank_point(rank) 

print("Form:", form1)
print("Low rank point:", point1)
print()

# generate second form
S = random_invertible(K,n)
print("equivalence S:")
print(S)
print()

S_inv = S.inverse()
form2 = act_on_form(form1,S)
point2 = S_inv*point1*random.randrange(1,q)

# find equivalence
S_prime = find_equivalence(form1,point1,form2,point2)
print("equivalence found: S_prime")
print(S_prime)
print()

# sanity check
print("Sanity check: S_prime * S^-1 should be a 3rd root of unity:")
print(S_prime*S_inv)
print()