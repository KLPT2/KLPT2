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
    ## Alternatively, we can return (u*ConjugateTranspose(u)).determinant();
    return a.reduced_norm() * d.reduced_norm() + b.reduced_norm() * c.reduced_norm() - (a.conjugate() * b * d.conjugate() * c).reduced_trace();

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
    if verbose: print("DONE!\nChoosePolarisationPrimePower: Compute first (a,c) using quadratic module structure... ");

    a, c = Compute_ac_QuadraticModuleStructure(O, g);

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
    ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
    gamma = ConjugateTranspose(tau)*ConjugateTranspose(u2)*Inverse(u1)*ReducedNorm(u1);

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
    ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
    gamma = ConjugateTranspose(tau)*ConjugateTranspose(u2)*Inverse(u1)*ReducedNorm(u1);

    return gamma, D^2*ReducedNorm(u1);

def FastExample():
    """
    This is a hardcoded example to show off the correctness
    of KLPT^2.

    g1 and g2 are computed using a cheat function that  
    returns reduced matrices. However, the a and c      
    that are hardcoded are from the Compute_ac function 
    Also, alpha from the KLPT algorithm has been        
    hardcoded as well. Note that these values can be    
    computed using Compute_ac and Compute_bd_KLPT too.  
    This just provides a proof-of-concept for the algo. 
    """
    p = 3 * 2^43 - 1;
    B.<i,j,k> = QuaternionAlgebra(-1,-p);
    O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])
    
    print(f"########    Starting example of size {p.log(2).n().floor()} bits    ########\n");
    print(f"########    g1 and g2 are computed using a cheat function that      ########");
    print(f"########    returns reduced matrices. However, the a and c          ########");
    print(f"########    that are hardcoded are from the Compute_ac function     ########");
    print(f"########    Also, alpha from the KLPT algorithm has been            ########");
    print(f"########    hardcoded as well. Note that these values can be        ########");
    print(f"########    computed using Compute_ac and Compute_bd_KLPT too.      ########");
    print(f"########    This just provides a proof-of-concept for the algo.     ########\n\n");

    L = 2
    g1 = Matrix(2,2,[50890016657, -863067848249479309881 + 1137349014214088243244*i + j - k,\
     -863067848249479309881 - 1137349014214088243244*i - j + k, 40082714730313402201920563137231]);
    s1 = ZZ(B(g1[0][0]));
    r1 = B(g1[0][1]);
    t1 = ZZ(B(g1[1][1]));
    a1 = 4300011028323335800307687665176258924009304382753121012 + 5679056451800766578285204460258507613777496939217999113*i
    c1 = -935637750462041355697723 - 1633022010301668477335021*i + 43950814457374960853024343929313*j - 57034733971562311049700603963456*k
    e1 = 400
    assert ZZ(s1 * a1.reduced_norm() + t1 * c1.reduced_norm() + (c1.conjugate()*r1.conjugate()*a1).reduced_trace()).is_prime_power(), "Error: a and c does not produce a top left entry which is a power of 2"
    
    ctx = KLPT_Context(B)
    I1 = O.left_ideal( [c1.reduced_norm(), c1 * a1.conjugate()] );
    # alpha1,J1,_,_,_,_ = ctx.KLPT(I1, T=2^400, returnElem=True);
    alpha1 = 21608308590885486802594499433628246787938189956653589712495470864358 + 85084780364216044308509191983820367957747057892469323651827799995351*i - 14714197354587346781157831555787798425708757042345800555829446*j + 3488217817734115868441001642408756739946031136488909232125207*k
    J1 = O.left_ideal([100433627766186892221372630771322662657637687111424552206336,\
     50216813883093446110686315385661331328818843555712276103168 + 50216813883093446110686315385661331328818843555712276103168*i,\
     8419462552381568417906684228168707061460691806934826892801 + 5031292796058782773067849868175706277206876879226405761504*i + j,\
     103821797522509677866211465131315663441891502039132973337633/2 + 13450755348440351190974534096344413338667568686161232654305/2*i + 1/2*j + 1/2*k])
    print(factor(J1.norm()))
    
    o1, o2 = Compute_o1o2(O, alpha1, a1, c1);
    _, aa, cc = xgcd(a1.reduced_norm(), c1.reduced_norm());
    b1 =  cc * ( c1.reduced_norm() * o1.conjugate() + a1 * c1.conjugate() * o2.conjugate());
    d1 = -aa * (c1 * a1.conjugate() * o1.conjugate() +  a1.reduced_norm() * o2.conjugate());
    
    u1 =  Matrix(2, 2, [B(a1), B(b1), B(c1), B(d1)]);
    h1 = ConjugateTranspose(u1)*g1*u1;
    print(f"Top left entry is a prime power: {ZZ(h1[0][0]).is_prime_power()}.")
    print(f"               Factorisation is: {factor(ZZ(h1[0][0]))}")
    print(f"     Norm of u is a prime power: {ZZ(ReducedNorm(u1)).is_prime_power()}.")
    print(f"               Factorisation is: {factor(ZZ(ReducedNorm(u1)))}")
    
    g2 = Matrix(2,2,[9649191169, -923954084130295985514 + 159836740039816258651*i + j - 15*k, -923954084130295985514 - 159836740039816258651*i - j + 15*k, 91261541728430194367609532073051]);
    s2 = ZZ(B(g2[0][0]));
    r2 = B(g2[0][1]);
    t2 = ZZ(B(g2[1][1]));
    a2 = 3198778801702244695722908677977591305436716001631634795 + 16043095261242602481185019370886986035806314135305502608*i
    c2 = -3691039507708977979394517 - 1511178102039228580798355*i - 1415642428282538826335037147749*j - 9320423832680276825962762278126*k
    e2 = 400
    assert ZZ(s2 * a2.reduced_norm() + t2 * c2.reduced_norm() + (c2.conjugate()*r2.conjugate()*a2).reduced_trace()).is_prime_power(), "Error: a and c does not produce a top left entry which is a power of 2"
    
    ctx = KLPT_Context(B)
    I2 = O.left_ideal( [c2.reduced_norm(), c2 * a2.conjugate()] );
    # alpha2,J2,_,_,_,_ = ctx.KLPT(I2, T=2^400, returnElem=True);
    alpha2 = 2618926897906107395758138216691929455610496108274702823676141167459 - 11317065132259730066657577240251116838894678255075433770558568122998*i - 1898697672799189830631347837377801503830527091671220146188867*j - 455505479499459925751606151599113528683863138415457638980806*k
    J2 = O.left_ideal([100433627766186892221372630771322662657637687111424552206336,\
     50216813883093446110686315385661331328818843555712276103168 + 50216813883093446110686315385661331328818843555712276103168*i,\
     5013472415338319236751372538437525161606967286536142727041 + 15543029969841702579936484545020812654490884448072322648464*i + j,\
     89904070211683508878187518764739375164753769949888372284913/2 + 20556502385180021816687857083458337816097851734608465375505/2*i + 1/2*j + 1/2*k])
    print(factor(J2.norm()))
    
    o1, o2 = Compute_o1o2(O, alpha2, a2, c2);
    _, aa, cc = xgcd(a2.reduced_norm(), c2.reduced_norm());
    b2 =  cc * ( c2.reduced_norm() * o1.conjugate() + a2 * c2.conjugate() * o2.conjugate());
    d2 = -aa * (c2 * a2.conjugate() * o1.conjugate() +  a2.reduced_norm() * o2.conjugate());
    
    u2 =  Matrix(2, 2, [B(a2), B(b2), B(c2), B(d2)]);
    h2 = ConjugateTranspose(u2)*g2*u2;
    print(f"Top left entry is a prime power: {ZZ(h2[0][0]).is_prime_power()}.")
    print(f"               Factorisation is: {factor(ZZ(h2[0][0]))}")
    print(f"     Norm of u is a prime power: {ZZ(ReducedNorm(u2)).is_prime_power()}.")
    print(f"               Factorisation is: {factor(ZZ(ReducedNorm(u2)))}")
    
    D = h1[0][0];
    assert D == h2[0][0], "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return the same first entries";
    assert ZZ(D).is_prime_power(), "[ConnectMatrices] Error: ChoosePolarisationPrimePower did not return prime power";
    
    ##  Constructing matrix tau, which is of the form:
    ##  [ D, r1-r2 ]
    ##  [ 0,   D   ]
    r1 = h1[0][1];
    t1 = h1[1][1];
    
    r2 = h2[0][1];
    t2 = h2[1][1];
    tau = Matrix(2,2,[D,r1-r2,0,D]);
    
    ##  Constructing matrix gamma which is given by
    ##  tau^* * u_2^* * u_1^{-1} * mathcal{N}(u_1)
    gamma = u2 * tau * MatInverse(u1) * ReducedNorm(u1)
    
    LHS = ConjugateTranspose(gamma)*g2*gamma;
    RHS = 2^1192 * g1;
    print(f"\n\n   gamma^* g2 gamma == ell^e g1: {LHS == RHS}")
    print(f"       e = 1192, L = 2")

FastExample()
