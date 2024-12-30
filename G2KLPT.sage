from klpt_panny import *

def _matrix_to_gens(M, B):
    """
    This piece of code was taken from the LearningToSQI repo.
    Converts from a matrix to generators in the quat.
    algebra
    """
    return [sum(c * g for c, g in zip(row, B)) for row in M]

def QuaternionInverse( a ):
    """
    Returns the inverse of a quaternion element
    """
    return a.conjugation()/a.reduced_norm()

def MatInverse( u ):
    """
    Returns the inverse of a matrix with quaternion elements
    """
    a = u[0][0];  b = u[0][1];
    c = u[1][0];  d = u[1][1];

    unorm = a.reduced_norm() * d.reduced_norm() + b.reduced_norm() * c.reduced_norm() - (a.conjugate()*b*d.conjugate()*c).reduced_trace()

    M = u * ConjugateTranspose(u)
    x = M[0][0];
    y = M[0][1];
    z = M[1][1];
    
    uinv = ConjugateTranspose(u) * Matrix(2,2,[z, -y, -(y.conjugate()), x]) / unorm
    return uinv

def RandomPolarisation(O, sbound=20):
    """
    This function returns a random polarisation matrix.
    """
    ctx = KLPT_Context(O.quaternion_algebra())
    s = ZZ(randint(0, (p ^ (sbound)).floor()));
    t = ZZ(randint(0, (p ^ (3 * sbound)).floor()));
    
    r = ctx.RepresentInteger(s * t - 1);
    return Matrix(2,2,[s,r,r.conjugate(),t]);

def RandomReduced(O, sbound=0.5):
    """
    This function returns a random reduced matrix where the definition
    of a reduced matrix is given in Definition 3.11 of the paper.
    """
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
    """
    This is the function that computes b and d using LLL.
    It takes as input a and c.
    """
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
    This is the reduced norm as given by Lemma 2.8 in the paper.
    mathcal{N}(u) = det(u*ConjugateTranspose(u)) 
                  = det(ConjugateTranspose(u)*u) 
                  = n(a)n(d) + n(b)n(c) - tr(Conjugate(a)*b*Conjugate(d)*c)
    """
    a = u[0][0];    b = u[0][1];
    c = u[1][0];    d = u[1][1];
    return a.reduced_norm() * d.reduced_norm() + b.reduced_norm() * c.reduced_norm() - (a.conjugate() * b * d.conjugate() * c).reduced_trace();

def RandomSL2OElement( O ):
    """
    Returns a random element of SL(2,O), where O is a 
    maximal quaternion order specified by the user.
    """
    a = O(1);
    b = O.random_element();
    c = b.conjugate();
    d = c*b + 1;
    return Matrix(2,2,[a,b,c,d]);

def RandomConjugation( M, O ):
    """
    Given a matrix M and a maximal quaternion order O, this
    returns the conjugation of M by a random SL(2,O) element.
    i.e. it return u^* M u
    """
    u = RandomSL2OElement(O);
    return ConjugateTranspose(u)*M*u;

def Compute_ac_LLL( O, g ):
    """
    This computes a and c using the LLL method.
    Given a polarisation matrix, this function works with the
    Compute_bd_KLPT function to produce a transformation matrix 
    U = [a b]  such that U^* g U = [s' *].
        [c d]                      [ * *]
    From the paper, we see that the target s' is given by
    s' = K((a, c)) := s * n(a) + t * n(c) + tr(cbar * rbar * a)
    which only depends on the a and c values of U, and also values
    from the initial polarisation matrix.
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

    Bgram  = Matrix(8,8,[     s,    0,    0,    0,  r1,  -r2,-p*r3,-p*r4,
                              0,    s,    0,    0,  r2,   r1,-p*r4, p*r3,
                              0,    0,  s*p,    0,p*r3, p*r4, p*r1,-p*r2,
                              0,    0,    0,  s*p,p*r4,-p*r3, p*r2, p*r1,
                             r1,   r2, p*r3, p*r4,   t,    0,    0,    0,
                            -r2,   r1, p*r4,-p*r3,   0,    t,    0,    0,
                          -p*r3,-p*r4, p*r1, p*r2,   0,    0,  t*p,    0,
                          -p*r4, p*r3,-p*r2, p*r1,   0,    0,    0,  t*p ]);
    
    if verbose: print("Done!\nCompute_ac_LLL: LLL... ");
    Q = Bgram.LLL_gram()

    a1 = B(Q.transpose()[0][0:4])
    c1 = B(Q.transpose()[0][4:8])
    a2 = B(Q.transpose()[1][0:4])
    c2 = B(Q.transpose()[1][4:8])
    a3 = B(Q.transpose()[2][0:4])
    c3 = B(Q.transpose()[2][4:8])
    a4 = B(Q.transpose()[3][0:4])
    c4 = B(Q.transpose()[3][4:8])
    
    
    if verbose:
        print(f"[Debug] Compute_ac_LLL: 1st short vector gives s' = {QuadraticForm_TopLeftEntry(a1,c1,g1)}")
        print(f"[Debug] Compute_ac_LLL: 2nd short vector gives s' = {QuadraticForm_TopLeftEntry(a2,c2,g1)}")
        print(f"[Debug] Compute_ac_LLL: 3rd short vector gives s' = {QuadraticForm_TopLeftEntry(a3,c3,g1)}")
        print(f"[Debug] Compute_ac_LLL: 4th short vector gives s' = {QuadraticForm_TopLeftEntry(a4,c4,g1)}")
    
    min_a = B(0)
    min_c = B(0)
    min_s = +Infinity
    
    if verbose: print("Finding smallest prime...")
    for ind in range(-10,10):
        for jnd in range(-10,10):
            for knd in range(-10,10):
                for lnd in range(-10,10):
                    if ind == jnd and jnd == knd and knd == lnd and lnd == 0:
                        continue;
                    a = ind*a1 + jnd*a2 + knd*a3 + lnd*a4
                    c = ind*c1 + jnd*c2 + knd*c3 + lnd*c4
                    na = a.reduced_norm()
                    nc = c.reduced_norm()
                    s = QuadraticForm_TopLeftEntry(a,c,g1)
                    if gcd(na,nc) == 1:
                        if ZZ(s).is_prime():
                            if ZZ(s) < min_s:
                                min_s = ZZ(s)
                                min_a = a
                                min_c = c
    if verbose: print("DONE!")
    assert min_s != 0, "Can't find a prime";
    
    return min_a, min_c;

def Compute_ac_QuadraticModuleStructure(O, g):
    """
    This is the old version of compute a and c.
    Keeping it here for now.
    """
    t2 = 2^(Ceiling(p.log(2).n()*10))
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
    """
    Given a polarisation matrix g, this function returns a and c
    such that 
     - n(a) and n(c) are co-prime
     - s * n(a) + t * n(c) + tr(c.conj * r.conj * a) = ell^e2
    Note that the output a and c will form the transformation 
    matrix u such that u^* g u is a matrix whose top-left entry
    is a power of L. Also, u = [a *] where the missing entries
                               [c *]
    will be found in another function called Compute_bd_KLPT.
    """
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

def Compute_o1o2( O, alpha, a, c ):
    """
    alpha is the generator of the principal ideal J * I.conjugate,
    and is obtained from the KLPT step. a and c are from the function
    Compute_ac.
    The purpose of this function is to return o1 and o2 such that 
    (o1 * c.reduced_norm() + o2 * c * a.conjugate()) == alpha
    This function does so by randomising o2 and computing o1 by
    making o1 the subject in the above equation.
    """
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];
    nc = c.reduced_norm()

    o2 = O.random_element();
    o1nc = alpha - o2*c*a.conjugate();
    return o1nc/nc, o2;
 
def Compute_o1o2_Integer( O, alpha, a, c ):
    """
    alpha is the generator of the principal ideal J * I.conjugate,
    and is obtained from the KLPT step. a and c are from the function
    Compute_ac.
    The purpose of this function is to return o1 and o2 such that 
    (o1 * c.reduced_norm() + o2 * c * a.conjugate()) == alpha
    This version of the function computes o1 and o2 with integer 
    coefficients in the basis <1,i,j,k>.
    """
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

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
    if alpha.denominator() == 2:
        beta = alpha - (1+i+j+k)/2
    elif alpha.denominator() == 1:
        beta = alpha
    (beta1, beta2, beta3, beta4) = beta.coefficient_tuple();
    o11, r1 = ZZ(beta1-oca1).quo_rem(nc)
    o12, r2 = ZZ(beta2-oca2).quo_rem(nc)
    o13, r3 = ZZ(beta3-oca3).quo_rem(nc)
    o14, r4 = ZZ(beta4-oca4).quo_rem(nc)
    assert r1==0 and r2==0 and r3==0 and r4==0, "[Compute_o1o2_simple] Error: o2 is not correct";
    if alpha.denominator() == 2:
        o1 = o11 + o12*i + o13*j + o14*k + (1+i+j+k)/2/nc
    elif alpha.denominator() == 1:
        o1 = o11 + o12*i + o13*j + o14*k
    assert (o1 * c.reduced_norm() + o2 * c * a.conjugate()) == alpha, "[Compute_o1o2_simple] Error: output does not meet condition";

    return o1, o2;
       
def Compute_o1o2_Random_o1( O, alpha, a, c ):
    """
    alpha is the generator of the principal ideal J * I.conjugate,
    and is obtained from the KLPT step. a and c are from the function
    Compute_ac.
    The purpose of this function is to return o1 and o2 such that 
    (o1 * c.reduced_norm() + o2 * c * a.conjugate()) == alpha
    This version of the function randomises o1 and computes o2
    by making o2 the subject.
    """
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];
    nc = c.reduced_norm()

    o1 = O.random_element();
    o2ca = alpha - o1*c.reduced_norm()
    return o1, o2ca*(c*a.conjugate()).inverse();

def Compute_bd_KLPT(O, a, c, L = 2):
    """
    This function is a continuation from the Compute_ac
    function and is used to obtain a matrix u = [a b] where 
                                                [c d]
    u^* g u is a matrix whose top-left entry is a power of L.
    Furthermore, this function ensures that the ReducedNorm
    of u is also a power of L.
    """
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
    Ialpha = O.left_ideal(alpha)
    
    o1, o2 = Compute_o1o2(O, alpha, a, c);
    _, aa, cc = xgcd(a.reduced_norm(), c.reduced_norm());
    b =  cc * ( c.reduced_norm() * o1.conjugate() + a * c.conjugate() * o2.conjugate());
    d = -aa * (c * a.conjugate() * o1.conjugate() +  a.reduced_norm() * o2.conjugate());
    if verbose:
        print("Done!\n");
    return b, d;

def FindAlpha( g, O ):
    """
    This function is used for the reduction step in the middle of 
    Section 3.3. This is to find the alpha in the [1 alpha]
                                                  [0     1]
    matrix in the paper.
    """
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

def ChoosePolarisationPrimePower( g, O, L = 2 ):
    """
    This is almost the full algorithm KLPT^2 as described in Algorithm 1.
    Given a matrix g, the KLPT^2 algorithm outputs the transformed
    matrix h, and also the transformation matrix u.
    This function covers lines 2-7 of the algorithm.
    """
    if verbose: print("\nChoosePolarisationPrimePower: Setting up... ");

    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    ####    Line 2 of Algorithm 1: Computing a and c of u to control the top left entry
    ####    This is described in the first part of Section 3.3
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute first (a,c) using LLL for to get a reduced matrix... ");

    a, c = Compute_ac_LLL(O, g);

    ####    Line 3 of Algorithm 1: Computing b and d of u to control the norm
    ####    and done according to Section 3.2 of the paper.
    ####    This is computed via the quadratic module structure and solving it using KLPT
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute first (b,d) using KLPT to get reduced matrix... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    u =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 4 for Algorithm 1: Computing alpha
    ####    This is really the intermediate step from the second part of Section 3.3.
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute alpha... ");

    gprime, alpha = FindAlpha(ConjugateTranspose(u)*g*u, O);

    v = Matrix(2, 2, [B(1), alpha, B(0), B(1)])

    g_red = ConjugateTranspose(v) * ConjugateTranspose(u) * g * u * v
    assert ZZ(g_red[0][0]).is_prime(), "Error: ChoosePolarisation failed to produce reduced matrix: s is not prime"
    assert ZZ(ReducedNorm(g_red)).is_prime_power(), "Error: ChoosePolarisation failed to produce reduced matrix: det is not a prime power"
    assert ZZ(r.reduced_norm())%2 == 1, "Error: ChoosePolarisation failed to produce reduced matrix: ell divides n(r)"

    ####    Line 5 for Algorithm 1: Computing a and c of uprime to control the top left entry.
    ####    Using QuadraticModuleStructure for now. 
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute second (a,c)... ");

    a, c = Compute_ac(O, gprime, t2);

    ####    Line 6 for Algorithm 1: Computing b and d of uprime to control the norm.
    ####    Using KLPT to compute this now.
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute second (b,d) using KLPT... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    uprime =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 7 for Algorithm 1: Computing transformation matrix
    ui = u*alpha*uprime;

    return ConjugateTranspose(ui)*g*ui, ui;

def ChoosePolarisationPrimePower_Reduced( g, O, t2, L = 2 ):
    """
    This is the reduced KLPT^2 algorithm which starts from a g
    which has been reduced as per Definition 3.11.
    Given such a matrix g, the reduced KLPT^2 algorithm outputs the 
    transformed matrix h, and also the transformation matrix u.
    This function covers lines 5-6 of the algorithm.
    """
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    ####    Line 5 for Algorithm 1: Computing a and c of uprime to control the top left entry.
    ####    Using QuadraticModuleStructure for now. 
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower_Reduced: Compute second (a,c) by solving norm equations... ");

    a, c = Compute_ac(O, g);

    ####    Line 6 for Algorithm 1: Computing b and d of uprime to control the norm.
    ####    Using KLPT to compute this now.
    ####    Can also use LLL.
    if verbose: print("DONE!\nChoosePolarisationPrimePower_Reduced: Compute second (b,d) using KLPT... ");

    b, d = Compute_bd_KLPT(O, B(a), B(c), L = 2 );
    uprime =  Matrix(2, 2, [B(a), B(b), B(c), B(d)]);

    ####    Line 7 for Algorithm 1: Computing transformation matrix
    u = uprime;

    return ConjugateTranspose(u)*g*u, u;

def ConnectMatrices( D, g1, g2, O ):
    """
    This is the full KLPT^2 algorithm.
    The bulk of the implementation of the algorithm is done in the above
    function. The implementation of Line 8 of Algorithm 1 in the paper
    is found here.
    It returns the transformation matrix gamma, and also l^e such that
    gamma^* g2 gamma = l^e g1
    """
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
    ##  u_2 * tau * u_1^{-1} * mathcal{N}(u_1)
    gamma = u2 * tau * MatInverse(u1) * ReducedNorm(u1)

    return gamma, D^2*ReducedNorm(u1);

def ConnectMatrices_Reduced( D, g1, g2, O ):
    """
    This is the reduced KLPT^2 algorithm.
    The bulk of the implementation of the algorithm is done in the above
    function. The implementation of Line 8 of Algorithm 1 in the paper
    is found here.
    It returns the transformation matrix gamma, and also l^e such that
    gamma^* g2 gamma = l^e g1
    """
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
    ##  u_2 * tau * u_1^{-1} * mathcal{N}(u_1)
    gamma = u2 * tau * MatInverse(u1) * ReducedNorm(u1)

    return gamma, D^2*ReducedNorm(u1);

def QuadraticForm_TopLeftEntry( a, c, G ):
    s = G[0][0];
    r = G[0][1];
    t = G[1][1];
    return s*a.reduced_norm() + t*c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()

def test_reduced_matrix():
    p = 3 * 2^43 - 1;
    B.<i,j,k> = QuaternionAlgebra(-1,-p);
    O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])
    
    print("Computing random polarisation with sbound = 10...");
    g1 = RandomPolarisation(O, sbound=10);
    a, c = Compute_ac_LLL(O, g1);
    b, d = Compute_bd_KLPT(O, a, c);
    
    print(f"n(a)*n(d) + n(c)*n(b) - tr(abar*b*dbar*c) is prime power: {ZZ(a.reduced_norm()*d.reduced_norm() + c.reduced_norm()*b.reduced_norm() - (a.conjugate()*b*d.conjugate()*c).reduced_trace()).is_prime_power()}")
    
    U = Matrix(2,2,[a,b,c,d])
    g = ConjugateTranspose(U)*g1*U;
    r = g[0][1]
    print(f" Top left is prime: {ZZ(g[0][0]).is_prime()}")
    print(f"det is prime power: {ZZ(ReducedNorm(g)).is_prime_power()}")
    print(f"     n(r) not even: {ZZ(r.reduced_norm())%2==1}")

