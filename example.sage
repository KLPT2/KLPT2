from klpt_panny import *

def _matrix_to_gens(M, B):
        """
        Converts from a matrix to generators in the quat.
        algebra
        """
        return [sum(c * g for c, g in zip(row, B)) for row in M]

def HermitianDeterminant( A ):
    a = A[0][0]
    b = A[0][1]
    c = A[1][0]
    d = A[1][1]
    return a.reduced_norm() * d.reduced_norm() + c.reduced_norm() * b.reduced_norm()\
         - (a.conjugate() * b * d.conjugate() * c).reduced_trace()

def QuaternionInverse( a ):
    return a.conjugation()/a.reduced_norm()

def MatOInverse( A ):
    B = base_ring(A)
    s = A[0][0]
    r = A[0][1]
    t = A[1][1]
    assert r == A[1][0].conjugate(), "[MatOInverse] Error: Matrix is not of the right form";
    det = s*t - r.reduced_norm();
    M = Matrix(2,2,[t,-r,-r.conjugate(),s]);
    return 1/det*M;

def MatBInverse( A ):
    B = base_ring(A)
    a = A[0][0]
    b = A[0][1]
    c = A[1][0]
    d = A[1][1]

    det = a.reduced_norm() * d.reduced_norm() + c.reduced_norm() * b.reduced_norm()\
        - (a.conjugate() * b * d.conjugate() * c).reduced_trace();
    M = Matrix(2,2,[d,-b,-c,a]);
    return 1/det*M;

def RepresentInteger(a,p):
    for z in range(1,500):
        for u in range(1,500):
            q = a - p * z ^ 2 - p * u ^ 2;
            if  q.is_prime()  and (q%4 == 1):
                try:
                    x, y = two_squares(a - p * z ^ 2 - p * u ^ 2)
                except:
                    continue;
                else:
                    return x,y,u,z;
    if verbose:
        print("Failed to represent integer");
    return 0;

def RandomSpecial(O, sbound=0.5):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];
    s = ZZ(randint(0, (p ^ sbound).floor()));
    t = ZZ(randint(0, (p ^ (3 * sbound)).floor()));
    while not (s.is_prime()):
        s = ZZ(randint(0,(p^sbound).floor()));
    while ((s * t - 1)%4) != 1:
        t = ZZ(randint(0, (p ^ (3 * sbound)).floor()));
    
    x,y,u,z = RepresentInteger(s * t - 1, p);
    r = z + u*i + x*j + y*i*j;
    return Matrix(2,2,[s,r,r.conjugate(),t]);

def RandomReduced(O, sbound=0.5):
    ctx = KLPT_Context(O.quaternion_algebra())
    s = ZZ(randint(0,(ctx.p^sbound).floor()));
    while not (s.is_prime() and s%4==1):
        s = ZZ(randint(0, (ctx.p^sbound).floor()));
    t = ZZ(randint(0, (ctx.p^(3*sbound)).floor()));
    while t%4 != 2:
        t = ZZ(randint(0, (ctx.p^(3*sbound)).floor()));
    for ind in range(1,1000):
        r = ctx.RepresentInteger(s * t - 2^ind)
        if r is None:
            continue;
        else:
            return Matrix(2,2,[s,r,r.conjugate(),t]);

def Compute_bd_LLL(a, c, O):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];
    P2.<b1,b2,b3,b4,d1,d2,d3,d4> = PolynomialRing(QQ,8);

    a1 = randint(0,p);
    a2 = randint(0,p);
    a3 = randint(0,p);
    a4 = randint(0,p);
    c1 = randint(0,p);
    c2 = randint(0,p);
    c3 = randint(0,p);
    c4 = randint(0,p);
    a  = a1 + a2*i + a3*j + a4*k;
    c  = c1 + c2*i + c3*j + c4*k;
    z1 = ZZ(a.reduced_norm());
    z2 = ZZ(c.reduced_norm());
   
    while GCD(z1,z2) != 1:
        a1 = randint(0,p);
        a2 = randint(0,p);
        a3 = randint(0,p);
        a4 = randint(0,p);
        c1 = randint(0,p);
        c2 = randint(0,p);
        c3 = randint(0,p);
        c4 = randint(0,p);
        a  = a1 + a2*i + a3*j + a4*k;
        c  = c1 + c2*i + c3*j + c4*k;
        z1 = ZZ(a.reduced_norm());
        z2 = ZZ(c.reduced_norm());
    f = a1^2*p*d3^2 + a1^2*p*d4^2 + a1^2*d1^2 + a1^2*d2^2 - 2*a1*c1*p*b3*d3 -\
        2*a1*c1*p*b4*d4 - 2*a1*c1*b1*d1 - 2*a1*c1*b2*d2 - 2*a1*c2*p*b3*d4 +\
        2*a1*c2*p*b4*d3 - 2*a1*c2*b1*d2 + 2*a1*c2*b2*d1 - 2*a1*c3*p*b1*d3 +\
        2*a1*c3*p*b2*d4 + 2*a1*c3*p*b3*d1 - 2*a1*c3*p*b4*d2 - 2*a1*c4*p*b1*d4 -\
        2*a1*c4*p*b2*d3 + 2*a1*c4*p*b3*d2 + 2*a1*c4*p*b4*d1 + a2^2*p*d3^2 +\
        a2^2*p*d4^2 + a2^2*d1^2 + a2^2*d2^2 + 2*a2*c1*p*b3*d4 - 2*a2*c1*p*b4*d3 +\
        2*a2*c1*b1*d2 - 2*a2*c1*b2*d1 - 2*a2*c2*p*b3*d3 - 2*a2*c2*p*b4*d4 -\
        2*a2*c2*b1*d1 - 2*a2*c2*b2*d2 - 2*a2*c3*p*b1*d4 - 2*a2*c3*p*b2*d3 +\
        2*a2*c3*p*b3*d2 + 2*a2*c3*p*b4*d1 + 2*a2*c4*p*b1*d3 - 2*a2*c4*p*b2*d4 -\
        2*a2*c4*p*b3*d1 + 2*a2*c4*p*b4*d2 + a3^2*p^2*d3^2 + a3^2*p^2*d4^2 +\
        a3^2*p*d1^2 + a3^2*p*d2^2 + 2*a3*c1*p*b1*d3 - 2*a3*c1*p*b2*d4 -\
        2*a3*c1*p*b3*d1 + 2*a3*c1*p*b4*d2 + 2*a3*c2*p*b1*d4 + 2*a3*c2*p*b2*d3 -\
        2*a3*c2*p*b3*d2 - 2*a3*c2*p*b4*d1 - 2*a3*c3*p^2*b3*d3 - 2*a3*c3*p^2*b4*d4 -\
        2*a3*c3*p*b1*d1 - 2*a3*c3*p*b2*d2 - 2*a3*c4*p^2*b3*d4 + 2*a3*c4*p^2*b4*d3 -\
        2*a3*c4*p*b1*d2 + 2*a3*c4*p*b2*d1 + a4^2*p^2*d3^2 + a4^2*p^2*d4^2 +\
        a4^2*p*d1^2 + a4^2*p*d2^2 + 2*a4*c1*p*b1*d4 + 2*a4*c1*p*b2*d3 -\
        2*a4*c1*p*b3*d2 - 2*a4*c1*p*b4*d1 - 2*a4*c2*p*b1*d3 + 2*a4*c2*p*b2*d4 +\
        2*a4*c2*p*b3*d1 - 2*a4*c2*p*b4*d2 + 2*a4*c3*p^2*b3*d4 - 2*a4*c3*p^2*b4*d3 +\
        2*a4*c3*p*b1*d2 - 2*a4*c3*p*b2*d1 - 2*a4*c4*p^2*b3*d3 - 2*a4*c4*p^2*b4*d4 -\
        2*a4*c4*p*b1*d1 - 2*a4*c4*p*b2*d2 + c1^2*p*b3^2 + c1^2*p*b4^2 + c1^2*b1^2 +\
        c1^2*b2^2 + c2^2*p*b3^2 + c2^2*p*b4^2 + c2^2*b1^2 + c2^2*b2^2 +\
        c3^2*p^2*b3^2 + c3^2*p^2*b4^2 + c3^2*p*b1^2 + c3^2*p*b2^2 + c4^2*p^2*b3^2 +\
        c4^2*p^2*b4^2 + c4^2*p*b1^2 + c4^2*p*b2^2;
    M = f.Gram_matrix();
    sv = M.LLL()[0];
    b = sv[0] + sv[1]*i + sv[2]*j + sv[3]*k;
    d = sv[4] + sv[5]*i + sv[6]*j + sv[7]*k;
    return b,d;      

def ConjugateTranspose( M ):
    r"""
    Computes conjugate and transpose of a given matrix M.
    Entry-wise conjugation is over the quaternion algebra.
    """
    return Matrix(2,2,[ind.conjugate() for ind in M]).transpose();

def ReducedNorm( u ):
    r"""
    Computes the reduced norm of a matrix which is defined to be
    mathcal{N}(u) = det(u*ConjugateTranspose(u)) 
                  = det(ConjugateTranspose(u)*u) 
                  = n(a)n(d) + n(b)n(c) - tr(Conjugate(a)*b*Conjugate(d)*c)
    """
    a = u[0][0];    b = u[0][1];
    c = u[1][0];    d = u[1][1];
    ## Alternatively, we can return (u*ConjugateTranspose(u)).determinant();
    return a.reduced_norm() * d.reduced_norm() + b.reduced_norm() * c.reduced_norm() - (a.conjugate() * b * d.conjugate() * c).reduced_trace();

def RandomSL2OElement( O ):
    r"""
    Returns a random element of SL(2,O), where O is a 
    maximal quaternion order specified by the user.
    """
    a = O(1);
    b = O.random_element();
    c = b.conjugate();
    d = c*b + 1;
    return Matrix(2,2,[a,b,c,d]);

def RandomConjugation( M, O ):
    r"""
    Given a matrix M and a maximal quaternion order O, this
    returns the conjugation of M by a random SL(2,O) element.
    i.e. it return u^* M u
    """
    u = RandomSL2OElement(O);
    return ConjugateTranspose(u)*M*u;

def Compute_ac_LLL( O, g, t2 ):
    r"""
    Use big matrix, and run LLL and that gives small s.
    """
    if verbose: print("Compute_ac_LLL: Setting up... ");

    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    s = ZZ(B(g[0][0]));
    r = B(g[0][1]);
    t = ZZ(B(g[1][1]));
    (r1, r2, r3, r4) = r.coefficient_tuple()
    if verbose: print("Done!\nCompute_ac_LLL: Setting up matrix... ");

    Q  = Matrix(8,8,[     s,    0,    0,    0,  r1,  -r2,-p*r3,-p*r4,
                          0,    s,    0,    0,  r2,   r1,-p*r4, p*r3,
                          0,    0,  s*p,    0,p*r3, p*r4, p*r1,-p*r2,
                          0,    0,    0,  s*p,p*r4,-p*r3, p*r2, p*r1,
                         r1,   r2, p*r3, p*r4,   t,    0,    0,    0,
                        -r2,   r1, p*r4,-p*r3,   0,    t,    0,    0,
                      -p*r3,-p*r4, p*r1, p*r2,   0,    0,  t*p,    0,
                      -p*r4, p*r3,-p*r2, p*r1,   0,    0,    0,  t*p ]);
    
    if verbose: print("Done!\nCompute_ac_LLL: LLL... ");
    RedQ, T, _ = Q.LLL();
    if verbose: print("Done!\nCompute_ac_LLL: Finding a and c... ");

    new_s = ZZ(RedQ[0][0]);
    ac = [ QQ(ind) for ind in T[1].coeffient_tuple() ];    
    ##  To check that a and c are correct, one can run this:
    ##print(Matrix(1,8,ac)*Q*Matrix(8,1,ac) eq new_s);
    a = B([ ac[0], ac[1], ac[2], ac[3] ]);
    c = B([ ac[4], ac[5], ac[6], ac[7] ]);

    gcd_ac, aprime, cprime = xgcd( ZZ(a.reduced_norm()), ZZ(c.reduced_norm()) );
    if gcd_ac != 1: return RedQ, [a,c];
    
    return RedQ, a, c;

def Compute_ac_QuadraticModuleStructure(O, g, t2):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    s = ZZ(B(g[0][0]));
    r = B(g[0][1]);
    t = ZZ(B(g[1][1]));
    (r1, r2, r3, r4) = r.coefficient_tuple()

    R = Integers(s);

    if verbose: print("Compute_ac: Computing matrix kernel... ");

    M   = Matrix(ZZ,2,4, [r1, -r2, -r3 * p, -r4 * p, r2, r1, -r4 * p, r3 * p]);
    M2  = M.transpose().kernel().basis_matrix();
    d00 = M2[0,0] + M2[0,1]*i + M2[0,2]*j + M2[0,3]*k;
    d01 = M2[1,0] + M2[1,1]*i + M2[1,2]*j + M2[1,3]*k;

    foundac = False;
    y = 1;
    if verbose: print("Done!\nCompute_ac: Starting while loop to find a and c... ");

    while not foundac:
        P2.<x> = PolynomialRing(B);
        r  = ZZ(R.random_element());
        u  = t2 - t * (x * d00 + r * d01) * (x * d00.conjugate() + r * d01.conjugate());
        u2 = R(u.monomial_coefficient(x^2));
        u1 = R(u.monomial_coefficient(x));
        u0 = R(u.constant_coefficient());
    
        P3.<l> = PolynomialRing(R);
        f = u2 * l^2 + u1 * l + u0;
        if not f.is_irreducible():
            y0 = y;
            x0 = ZZ(f.roots()[0][0]);
            w  = t2 - ZZ(t * (x0 * d00 + y0 * d01) * (x0 * d00.conjugate() + y0 * d01.conjugate()));
            w1, _ = w.quo_rem(s);
            if ZZ(w1).is_prime():
                try:
                    two_squares(w1);
                except:
                    continue;
                else:
                    foundac = true;
                    a = x0 * d00 + y0 * d01;
                    c = w1;
                    return a, c;
                if gcd(ZZ(a.reduced_norm()), ZZ(c.reduced_norm())) != 1:
                    foundac = false;
        y += 1;

def Compute_ac( O, g, L=2 ):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    s = ZZ(B(g[0][0]));
    r = B(g[0][1]);
    t = ZZ(B(g[1][1]));
    (r1, r2, r3, r4) = r.coefficient_tuple()

    ####    Solve t * p * n(r) * (c1^2 + c2^2) = L^e mod s
    K = GF(s);
    c2 = K.random_element();
    ebound = ZZ((p.log(2).n() * 7).floor())
    for e in range(ebound):
        if (L^e/K(ZZ(t*p*r.reduced_norm())) - c2^2).is_square():
            #print(f"Square at e={e}")
            c1 = (L^e/K(ZZ(t*p*r.reduced_norm())) - c2^2).sqrt();
            c = ZZ(c1) * r.conjugate() * j + ZZ(c2) * r.conjugate() * k;
            le_tnc = ZZ(L^e - t * c.reduced_norm());
            lhs, rem = le_tnc.quo_rem(s);
            if rem != 0: continue;
            else:
                #print(f"Rem=0 at e={e}")
                try:
                    a1, a2 = two_squares(lhs);
                except:
                    continue;
                else:
                    a = a1 + a2*i;
                    return a, c;

def Compute_o1o2_simple( O, alpha, a, c ):
    assert alpha.denominator() == 2, "[Compute_o1o2_simple] Using the wrong Compute_o1o2 function"
    (alpha1, alpha2, alpha3, alpha4) = alpha.coefficient_tuple();
    nc = c.reduced_norm();
    R = Integers(nc)
    gamma = c * a.conjugate();
    (gamma1, gamma2, gamma3, gamma4) = gamma.coefficient_tuple();
    M  = Matrix(4,4,[ R(gamma1), R(-gamma2), R(-gamma3 * p), R(-gamma4 * p),\
                      R(gamma2), R( gamma1), R( gamma4 * p), R(-gamma3 * p),\
                      R(gamma3), R(-gamma4), R( gamma1),     R( gamma2),\
                      R(gamma4), R( gamma3), R(-gamma2),     R( gamma1)]);
    V  = vector([R(alpha1), R(alpha2), R(alpha3), R(alpha4)]);
    o2vec = M.solve_right(V)
    o2 = ZZ(o2vec[0]) + ZZ(o2vec[1])*i + ZZ(o2vec[2])*j + ZZ(o2vec[3])*k
    oca = o2*gamma;
    (oca1, oca2, oca3, oca4) = oca.coefficient_tuple();
    beta = alpha - (1+i+j+k)/2
    (beta1, beta2, beta3, beta4) = beta.coefficient_tuple();
    o11, r1 = ZZ(beta1-oca1).quo_rem(nc)
    o12, r2 = ZZ(beta2-oca2).quo_rem(nc)
    o13, r3 = ZZ(beta3-oca3).quo_rem(nc)
    o14, r4 = ZZ(beta4-oca4).quo_rem(nc)
    assert r1==0 and r2==0 and r3==0 and r4==0, "[Compute_o1o2_simple] Error: o2 is not correct";
    o1 = o11 + o12*i + o13*j + o14*k + (1+i+j+k)/2/nc
    assert (o1 * c.reduced_norm() + o2 * c * a.conjugate()) == alpha, "[Compute_o1o2_simple] Error: output does not meet condition";

    return o1, o2;

def Compute_o1o2( O, alpha, a, c ):
    r"""
    alpha = o1 * N(c) + o2 * c * a.conjugate()
    o1 = w1 + x1*i + y1 * (i+j)/2 + z1 * (1+k)/2
    o2 = w2 + x2*i + y2 * (i+j)/2 + z2 * (1+k)/2
    and solve alpha = o1 * N(c) + o2 * c * a.conjugate()

    c * a.conjugate = gamma1 + gamma2 * i + gamma3 * j + gamma4 * k
    o1 * n(c)       = n(c) * (w1 + z1/2) + n(c) * (x1 + y1/2) * i + n(c) * y1/2 * j + n(c) * z1/2 * k
    i * c * a.conj  = -gamma2     + gamma1 * i     - gamma4 * j + gamma3 * k
    j * c * a.conj  = -gamma3 * p + gamma4 * p * i - gamma1 * j + gamma2 * k
    k * c * a.conj  = -gamma4 * p - gamma3 * p * i - gamma2 * j + gamma3 * k
    o2 * c * a.conj = [ (w2 + z2/2) * gamma1 - (x2 + y2/2) * gamma2 -    p * y2/2 * gamma3 - p * z2/2 * gamma4 ]
                + i * [ (x2 + y2/2) * gamma1 + (w2 + z2/2) * gamma2 -    p * z2/2 * gamma3 + p * y2/2 * gamma4 ]
                + j * [ (     y2/2) * gamma1 + (     z2/2) * gamma2 + (w2 + z2/2) * gamma3 - (x2 + y2/2) * gamma4 ]
                + k * [ (     z2/2) * gamma1 - (     y2/2) * gamma2 + (x2 + y2/2) * gamma3 + (w2 + z2/2) * gamma4 ]
    
    alpha1 = (w2 + z2/2) * gamma1 - (x2 + y2/2) * gamma2 -    p * y2/2 * gamma3 - p * z2/2 * gamma4    + n(c) * (w1 + z1/2)
    alpha2 = (x2 + y2/2) * gamma1 + (w2 + z2/2) * gamma2 -    p * z2/2 * gamma3 + p * y2/2 * gamma4    + n(c) * (x1 + y1/2)
    alpha3 = (     y2/2) * gamma1 + (     z2/2) * gamma2 + (w2 + z2/2) * gamma3 - (x2 + y2/2) * gamma4 + n(c) * y1/2
    alpha4 = (     z2/2) * gamma1 - (     y2/2) * gamma2 + (x2 + y2/2) * gamma3 + (w2 + z2/2) * gamma4 + n(c) * z1/2
    """

    ####    Solve for o2 first by solving it mod n(c)
    (alpha1, alpha2, alpha3, alpha4) = alpha.coefficient_tuple();
    nc = c.reduced_norm();
    # assert nc%2 == 0, "[Compute_o1o2] Odd c norm not implemented"
    if nc%2==0:
        nc2 = nc//2
    else:
        nc2 = nc
    R  = Integers(nc2);
    gamma = c * a.conjugate();
    (gamma1, gamma2, gamma3, gamma4) = gamma.coefficient_tuple();
    M  = Matrix(4,4,[ R(gamma1), R(-gamma2), R(-(gamma2/2 + p * gamma3/2)), R(gamma1/2 - p * gamma4/2),\
                      R(gamma2), R( gamma1), R( (gamma1/2 - p * gamma3/2)), R(gamma2/2 - p * gamma4/2),\
                      R(gamma3), R(-gamma4), R( (gamma1/2 - gamma4/2)),     R(gamma2/2 + gamma3/2),\
                      R(gamma4), R( gamma3), R( (gamma3/2 - gamma2/2)),     R(gamma1/2 + gamma4/2)]);
    V  = vector([R(alpha1), R(alpha2), R(alpha3), R(alpha4)]);
    o2vec = M.solve_right(V)
    o2 = ZZ(o2vec[0]) + ZZ(o2vec[1]) * i + ZZ(o2vec[2]) * (i/2+j/2) + ZZ(o2vec[3]) * (1/2+k/2)

    oca = o2*gamma;
    (oca1, oca2, oca3, oca4) = oca.coefficient_tuple();
    o11, r1 = ZZ(alpha1-oca1).quo_rem(nc)
    o12, r2 = ZZ(alpha2-oca2).quo_rem(nc)
    o13, r3 = ZZ(alpha3-oca3).quo_rem(nc)
    o14, r4 = ZZ(alpha4-oca4).quo_rem(nc)

    assert r1==0 and r2==0 and r3==0 and r4==0, "[Compute_o1o2] Error: o2 is not correct";
    o1 = (o11-o14) + (o12-o13)*i + 2*o13*(i+j)/2 + 2*o14*(1+k)/2;
    assert o1 * c.reduced_norm() + o2 * c * a.conjugate() == alpha, "[Compute_o1o2] Error: output does not meet condition";
    return o1, o2;

def Compute_bd_KLPT(O, a, c, L = 2):
    B = O.quaternion_algebra();
    if verbose:
        print("Compute_bd: KLPT Context... ");
    ctx = KLPT_Context(B)
    I = O.left_ideal( [c.reduced_norm(), c * a.conjugate()] );
    if verbose:
        print("DONE!\nCompute_bd: KLPT... ");
    alpha,J,_,_,_,_ = ctx.KLPT(I,T=2^200,returnElem=True);
    if verbose:
        print("Done!\nCompute_bd: Computing b and d... ");

    ####    Note that this is not the principle generator, so we need to find alpha first.
    ####    Do this by getting the Gram matrix of the generators and then use LLL.

    JcI = J.conjugate() * I
    JcI_basis = JcI.basis()
    JcI_gram  = JcI.gram_matrix()
    JcI_gens = _matrix_to_gens(JcI_gram.LLL_gram().transpose(), JcI_basis)
    alphavec = JcI_gens[0]
    alpha = B(alphavec)
    print(f"Checking norm of J: {J.norm()}")
    o1, o2 = Compute_o1o2(O, alpha, a, c);
    _, aa, cc = xgcd(a.reduced_norm(), c.reduced_norm());
    b =  cc * (c.reduced_norm() * o1 + a * c.conjugate() * o2);
    d = -aa * (c * a.conjugate() * o1 + a.reduced_norm() * o2);
    if verbose:
        print("Done!\n");
    return b, d;

def FindAlpha( g, O ):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    s = g[0][0];
    r = g[0][1];
    t = g[1][1];
    (r1, r2, r3, r4) = r.coefficient_tuple()

    alpha1, r1new = ZZ(r1).quo_rem(ZZ(s));
    alpha2, r2new = ZZ(r2).quo_rem(ZZ(s));
    alpha3, r3new = ZZ(r3).quo_rem(ZZ(s));
    alpha4, r4new = ZZ(r4).quo_rem(ZZ(s));

    alpha = B([alpha1, alpha2, alpha3, alpha4]);
    rnew  = B([r1new, r2new, r3new, r4new]);
    rnew_check = r - alpha*s;
    tnew  = alpha.reduced_norm() * s + (alpha.conjugate() * r).reduced_trace() + t;

    assert rnew == rnew_check, "[FindAlpha] Error: r entry has been computed wrongly";

    return Matrix(2,2,[s,rnew,rnew.conjugate(),tnew]), Matrix(2,2,[1,alpha,0,1]);

####    This is the bulk of the algorithm
def ChoosePolarisationPrimePower( g, O, t2, L = 2 ):
    if verbose: print("\nChoosePolarisationPrimePower: Setting up... ");

    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    ####    Line 2 of Algorithm 1: Computing a and c of u to control the top left entry
    ####    This is described in the first part of Section 3.3
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute first (a,c) using quadratic module structure... ");

    a, c = Compute_ac_QuadraticModuleStructure(O, g, t2);

    ####    Line 3 of Algorithm 1: Computing b and d of u to control the norm
    ####    and done according to Section 3.2 of the paper.
    ####    This is computed via the quadratic module structure and solving it using KLPT
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute first (b,d) using KLPT... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    u =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 4 for Algorithm 1: Computing alpha
    ####    This is really the intermediate step from the second part of Section 3.3.
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute alpha... ");

    gprime, alpha = FindAlpha(ConjugateTranspose(u)*g*u, O);

    ####    Line 5 for Algorithm 1: Computing a and c of uprime to control the top left entry.
    ####    Using QuadraticModuleStructure for now. 
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute second (a,c) using quadratic module structure... ");

    a, c = Compute_ac_QuadraticModuleStructure(O, gprime, t2);

    ####    Line 6 for Algorithm 1: Computing b and d of uprime to control the norm.
    ####    Using KLPT to compute this now.
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower_V2: Compute second (b,d) using KLPT... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    uprime =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 7 for Algorithm 1: Computing transformation matrix
    ui = u*alpha*uprime;

    return ui*g*ConjugateTranspose(ui), ui;

def ChoosePolarisationPrimePower_Reduced( g, O, t2, L = 2 ):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    ####    Line 5 for Algorithm 1: Computing a and c of uprime to control the top left entry.
    ####    Using QuadraticModuleStructure for now. 
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower_Reduced: Compute second (a,c) by solving norm equations... ");

    a, c = Compute_ac(O, g, t2);

    ####    Line 6 for Algorithm 1: Computing b and d of uprime to control the norm.
    ####    Using KLPT to compute this now.
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower_Reduced: Compute second (b,d) using KLPT... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    uprime =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 7 for Algorithm 1: Computing transformation matrix
    ui = uprime;

    return ConjugateTranspose(ui)*g*ui, ui;

def ConnectMatrices( D, g1, g2, O ):
    assert g1[0][1].conjugate() == g1[1][0], "[ConnectMatrices] Error: g1 not of the correct form";
    assert g2[0][1].conjugate() == g2[1][0], "[ConnectMatrices] Error: g2 not of the correct form";
    
    ##  This method uses ChoosePolarisation to find g1 and g2 that have
    ##  the same first entry which has been specially chosen.
    h1, u1 = ChoosePolarisationPrimePower( g1, O, t2, L=2 );
    h2, u2 = ChoosePolarisationPrimePower( g2, O, t2, L=2 );
    D = h1[0][0];
    assert D == h2[0][0], "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return the same first entries";
    assert D.is_power(), "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return prime power";

    ##  Constructing matrix tau, which is of the form:
    ##  [ D, r1-r2 ]
    ##  [ 0,   D   ]
    r1 = h1[0][1];
    r2 = h2[0][1];
    tau = Matrix(2,2,[D,r1-r2,0,D]);
    
    ##  Constructing matrix gamma which is given by
    ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
    gamma = ConjugateTranspose(tau)*ConjugateTranspose(u2)*Inverse(u1)*ReducedNorm(u1);

    return gamma, D^2*ReducedNorm(u1);

def ConnectMatrices_Reduced( D, g1, g2, O ):
    assert g1[0][1].conjugate() == g1[1][0], "[ConnectMatrices] Error: g1 not of the correct form";
    assert g2[0][1].conjugate() == g2[1][0], "[ConnectMatrices] Error: g2 not of the correct form";
    
    ##  This method uses ChoosePolarisation to find g1 and g2 that have
    ##  the same first entry which has been specially chosen.
    h1, u1 = ChoosePolarisationPrimePower_Reduced( g1, O, t2, L=2 );
    h2, u2 = ChoosePolarisationPrimePower_Reduced( g2, O, t2, L=2 );
    D = h1[0][0];
    assert D == h2[0][0], "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return the same first entries";
    assert D.is_power(), "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return prime power";

    ##  Constructing matrix tau, which is of the form:
    ##  [ D, r1-r2 ]
    ##  [ 0,   D   ]
    r1 = h1[0][1];
    r2 = h2[0][1];
    tau = Matrix(2,2,[D,r1-r2,0,D]);
    
    ##  Constructing matrix gamma which is given by
    ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
    gamma = ConjugateTranspose(tau)*ConjugateTranspose(u2)*Inverse(u1)*ReducedNorm(u1);

    return gamma, D^2*ReducedNorm(u1);


p = 3 * 2^43 - 1;
B.<i,j,k> = QuaternionAlgebra(-1,-p);
# O = B.maximal_order()
O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])

print(f"########    Starting example of size {p.log(2).n().floor()} bits    ########\n");
print(f"########    g1 and g2 are computed using a cheat function that      ########");
print(f"########    returns reduced matrices. However, the a and c          ########");
print(f"########    that are hardcoded are from the Compute_ac function     ########");
print(f"########    The Compute_ac function doesn't work in all cases,      ########");
print(f"########    and I have to look into it. But this example works      ########");
print(f"########    well enough to test the KLPT step.                      ########\n\n");
L = 2
g1 = Matrix(2,2,[3415397, -1428520487029 + 1525469341101*i + j, -1428520487029 - 1525469341101*i - j, 23932497436714778957]);
s = ZZ(B(g1[0][0]));
r = B(g1[0][1]);
t = ZZ(B(g1[1][1]));
a = 13372821708437065250433744034046816469522433 + 47582915196688534878375296096155410984642568*i
c = 33184923387859619749 - 41150937907654571120*i + 582423400056893113*j - 4146065789295500623*k
e = 312;
assert ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power(), "Error: a and c does not produce a top left entry which is a power of 2"

ctx = KLPT_Context(B)
I = O.left_ideal( [c.reduced_norm(), c * a.conjugate()] );
# alpha,J,_,_,_,_ = ctx.KLPT(I, T=2^400, returnElem=True);
alpha = -25689765571812964472327497785708282685072631698985562624892464707272357710068872183841/2 + 9076128410896141384668373733345378105307090893606155360536604946024341999140081116661/2*i - 388467772179038320933806742074405547675085872567584835324855186801495527503819/2*j - 4120936610677708668968160584231689139108083183242852177033477324650799901962015/2*k;
J = O.left_ideal([645562469521727147413979793000752968582426448207305878207664839135161905504210298657411338320034457858975792993186873344,\
322781234760863573706989896500376484291213224103652939103832419567580952752105149328705669160017228929487896496593436672 + 322781234760863573706989896500376484291213224103652939103832419567580952752105149328705669160017228929487896496593436672*i,\
603091228545843506322930216376321629697778935433741436249788558490402611995564783041562066413922508557243964607691497113 + 219897151776359039437561934344299669110461485943847084307376098563527890230602333158560967616206176097020949573180912396*i + j,\
1028756546291211614299348075032774929169743897697200230150077299062036627269172748540412437117750790319198808027697458061/2 + 177425910800475398346512357719868330225813973170282642349499817918768596721956817542711695710094226795289121187685536165/2*i + 1/2*j + 1/2*k])
print(factor(J.norm()))
JcI = J.conjugate() * I
Ialpha = O.left_ideal(alpha)

o1, o2 = Compute_o1o2_simple(O, alpha, a, c);
_, aa, cc = xgcd(a.reduced_norm(), c.reduced_norm());
b =  cc * (c.reduced_norm() * o1 + a * c.conjugate() * o2);
d = -aa * (c * a.conjugate() * o1 + a.reduced_norm() * o2);

uprime =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);
u1 = uprime;
h1 = ConjugateTranspose(u1)*g1*u1;
print(ZZ(h1[0][0]).is_prime_power())
na = a.reduced_norm()
nb = b.reduced_norm()
nc = c.reduced_norm()
nd = d.reduced_norm()
Tabdc = (a.conjugate()*b*d.conjugate()*c).reduced_trace()
print(ZZ(na*nd + nc*nb - Tabdc).is_prime_power())

# b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
# u1 =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);
# h1 = ConjugateTranspose(u1)*g1*u1;
# 
# 
# g2 = Matrix(2,2,[2439293, -13806033694566 - 7520372325196*i + j - 8*k, -13806033694566 + 7520372325196*i - j + 8*k, 133044213525032504367]);
# s = ZZ(B(g2[0][0]));
# r = B(g2[0][1]);
# t = ZZ(B(g2[1][1]));
# a = 35750393660878712955326929919887252786631683 + 46286600504200301725430876607031515132494632*i
# c = -124091621182608522944 - 69704401659973054893*i - 8147814418960481716*j - 6650759827868279666*k
# e  = 312
# b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
# u2 =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);
# h2 = ConjugateTranspose(u2)*g2*u2;
# 
# D = h1[0][0];
# assert D == h2[0][0], "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return the same first entries";
# assert ZZ(D).is_prime_power(), "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return prime power";
# 
# ##  Constructing matrix tau, which is of the form:
# ##  [ D, r1-r2 ]
# ##  [ 0,   D   ]
# r1 = h1[0][1];
# r2 = h2[0][1];
# tau = Matrix(2,2,[D,r1-r2,0,D]);
# 
# ##  Constructing matrix gamma which is given by
# ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
# gamma = ConjugateTranspose(tau)*ConjugateTranspose(u2)*MatBInverse(u1)*ReducedNorm(u1);
