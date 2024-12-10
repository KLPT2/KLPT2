from klpt_panny import *

def twoadic(N):
    M = N
    while (N % 2) == 0:
        N = N//2;
    return (M // N);

def Compute_ac( O, g, L=2, e=None ):
    temp = O.quaternion_algebra();
    B.<i,j,k> = temp;
    p = B.ramified_primes()[0];

    s = ZZ(B(g[0][0]));
    r = B(g[0][1]);
    t = ZZ(B(g[1][1]));
    (r1, r2, r3, r4) = r.coefficient_tuple()

    ####    Solve t * p * n(r) * (c1^2 + c2^2) = L^e mod s
    K = GF(s);
    if e is None:
        e= ZZ((p.log(2).n() * 7).floor())
    print(f" e = {e}");
    Le = 2^e
    foundsol = False
    u = K(ZZ(t*p*r.reduced_norm()))
    while not foundsol:
        c20 = K.random_element();
        if (Le * u^(-1) - c20^2).is_square():
            c10 = (Le * u^(-1) - c20^2).sqrt();
            c2 = ZZ(c20);
            c1 = ZZ(c10);
            if ((c1-c2) % 2) == 1:
                c  = c1 * r.conjugate() * j + c2 * r.conjugate() * k;
                A0 = ZZ(Le-t*c.reduced_norm());
                A1 = A0 // s;
                A  = A1 // twoadic(A1);
                if A.is_prime() and A%4 == 1:
                    a1, a2 = two_squares(A)
                    foundsol = True
    a = a1 + a2*i;
    return a, c, e

def RandomReduced(O, sbound=0.5):
    ctx = KLPT_Context(O.quaternion_algebra())
    s = ZZ(randint(0,(ctx.p^sbound).floor()));
    while not (s.is_prime() and s%4==1):
        s = ZZ(randint(0, (ctx.p^sbound).floor()));
    t = ZZ(randint(0, (ctx.p^(3*sbound)).floor()));
    while not t%2==1:
        t = ZZ(randint(0, (ctx.p^(3*sbound)).floor()));
    for ind in range(1,1000):
        r = ctx.RepresentInteger(s * t - 2^ind)
        if r is None:
            continue;
        else:
            return Matrix(2,2,[s,r,r.conjugate(),t]);

def MoreRandomCheat(e=None):
    p = 3 * 2^43 - 1;
    B.<i,j,k> = QuaternionAlgebra(-1,-p);
    O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])

    ####    Solve t * p * n(r) * (c1^2 + c2^2) = L^e mod s
    if e is None:
        e = ZZ((p.log(2).n() * 7).floor())
    print(f" e = {e}");
    Le = 2^e
    foundsol = False
    while not foundsol:
        g = RandomReduced(O, sbound=0.5);
        s = ZZ(B(g[0][0]));
        r = B(g[0][1]);
        t = ZZ(B(g[1][1]));
        K = GF(s);
        u = K(ZZ(t*p*r.reduced_norm()))

        c20 = K.random_element();
        if (Le * u^(-1) - c20^2).is_square():
            c10 = (Le * u^(-1) - c20^2).sqrt();
            c2 = ZZ(c20);
            c1 = ZZ(c10);
            if ((c1-c2) % 2) == 1:
                c  = c1 * r.conjugate() * j + c2 * r.conjugate() * k;
                A0 = ZZ(Le-t*c.reduced_norm());
                A1 = A0 // s;
                A  = A1 // twoadic(A1);
                if A.is_prime() and A%4 == 1:
                    a1, a2 = two_squares(A)
                    foundsol = True
    a = a1 + a2*i;
    print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
    print(f"           g   = [{g[0,0]}, {g[0,1]}]")
    print(f"                 [{g[1,0]}, {g[1,1]}]")
    print(f"           a   = {a}")
    print(f"           c   = {c}")
    print(f"           e   = {e}")
    print(f"       s mod 4 = {s%4}")
    print(f"       t mod 4 = {t%4}")
    print(f"top left cond  = {ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power()}");

def GenerateCheatPairs(e=None):
    p = 3 * 2^43 - 1;
    B.<i,j,k> = QuaternionAlgebra(-1,-p);
    O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])

    if e is None:
        e = ZZ((p.log(2).n() * 7).floor())

    g1 = RandomReduced(O, sbound=0.5);
    s = ZZ(g1[0][0]);
    r = g1[0][1];
    t = ZZ(g1[1][1]);
    K = GF(s);
    a, c, e = Compute_ac( O, g1, e=e )
    print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
    print(f"           g1  = [{g1[0,0]}, {g1[0,1]}]")
    print(f"                 [{g1[1,0]}, {g1[1,1]}]")
    print(f"           a   = {a}")
    print(f"           c   = {c}")
    print(f"           e   = {e}")
    print(f"top left cond  = {ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power()}");

    g2 = RandomReduced(O, sbound=0.5);
    s = ZZ(g1[0][0]);
    r = g1[0][1];
    t = ZZ(g1[1][1]);
    K = GF(s);
    a, c, e = Compute_ac( O, g2, e=e )
    print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
    print(f"           g2  = [{g2[0,0]}, {g2[0,1]}]")
    print(f"                 [{g2[1,0]}, {g2[1,1]}]")
    print(f"           a   = {a}")
    print(f"           c   = {c}")
    print(f"           e   = {e}")
    print(f"top left cond  = {ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power()}");

GenerateCheatPairs(e=1000)
#MoreRandomCheat(e=1000)



#Le = 2^1000;
#p  = 3 * 2^43 - 1;
#B.<i,j,k> = QuaternionAlgebra(-1,-p);
#O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])
#
#r = 183473551550885215239093361195645100494348123432123+103489357288336036570742620116745037321594270065985*i+5125*j+6646*k;
#s = 5421341;
#t = 8184799884502414753566140534100785841846441788692629529322485064161801624924877657812355294081;
#R = ZZ(r.reduced_norm())
#print(s*t-R);
#K = GF(s)
#u = K(t*p*R);
#foundsol = False
#while not foundsol:
#    c20 = K.random_element()
#    if (Le*(u^(-1))-c20^2).is_square():
#        c10 = (Le*(u^(-1))-c20^2).sqrt();
#        c2  = ZZ(c20);
#        c1  = ZZ(c10);
#        if ((c1-c2) % 2) == 1:
#            c = c1*r.conjugate()*j + c2*r.conjugate()*k;
#            A0 = ZZ(Le-t*c.reduced_norm());
#            A1 = A0//s;
#            A  = A1//twoadic(A1);
#            if A.is_prime() and A%4==1:
#                a1, a2 = two_squares(A)
#                foundsol = True
#a = a1 + a2*i
#print(s%4)
#print(t%4)
#print(ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power())
