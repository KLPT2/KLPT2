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
    for ind in range(86,1000):
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

    print("Generating g1...")
    g1 = RandomReduced(O, sbound=0.5);
    print("DONE!\nGenerating a and c...")
    s = ZZ(g1[0][0]);
    r = g1[0][1];
    t = ZZ(g1[1][1]);
    K = GF(s);
    a, c, e = Compute_ac( O, g1, e=e )
    print("DONE!\nParameters generated:\n")
    print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
    print(f"           g1  = [{g1[0,0]}, {g1[0,1]}]")
    print(f"                 [{g1[1,0]}, {g1[1,1]}]")
    print(f"           a   = {a}")
    print(f"           c   = {c}")
    print(f"           e   = {e}")
    print(f"       s*t - R = {s*t-r.reduced_norm()}")
    print(f"       s mod 4 = {s%4}")
    print(f"       t mod 4 = {t%4}")
    print(f"top left cond  = {ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power()}");

    print("\nGenerating g2...")
    g2 = RandomReduced(O, sbound=0.5);
    print("DONE!\nGenerating a and c...")
    s = ZZ(g2[0][0]);
    r = g2[0][1];
    t = ZZ(g2[1][1]);
    K = GF(s);
    a, c, e = Compute_ac( O, g2, e=e )
    print("DONE!\nParameters generated:\n")
    print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
    print(f"           g2  = [{g2[0,0]}, {g2[0,1]}]")
    print(f"                 [{g2[1,0]}, {g2[1,1]}]")
    print(f"           a   = {a}")
    print(f"           c   = {c}")
    print(f"           e   = {e}")
    print(f"       s*t - R = {s*t-r.reduced_norm()}")
    print(f"       s mod 4 = {s%4}")
    print(f"       t mod 4 = {t%4}")
    print(f"top left cond  = {ZZ(s * a.reduced_norm() + t * c.reduced_norm() + (c.conjugate()*r.conjugate()*a).reduced_trace()).is_prime_power()}");

GenerateCheatPairs()
