from klpt_panny import *

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
    ebound = ZZ((p.log(2).n() * 7).floor())
    for e in range(ebound):
        c2 = K.random_element();
        if (L^e/K(ZZ(t*p*r.reduced_norm())) - c2^2).is_square():
            #print(f"Square at e={e}")
            c1 = (L^e/K(ZZ(t*p*r.reduced_norm())) - c2^2).sqrt();
            c = ZZ(c1) * r.conjugate() * j + ZZ(c2) * r.conjugate() * k;
            le_tnc = ZZ(L^e - t * c.reduced_norm());
            lhs, rem = le_tnc.quo_rem(s);
            if rem != 0: continue;
            #print(f"Rem=0 at e={e}")
            if gcd(c.reduced_norm(),lhs) != 1: continue;
            try:
                a1, a2 = two_squares(lhs);
            except:
                continue;
            else:
                a = a1 + a2*i;
                return [a, c];
    return None

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

def GenerateCheatPairs():
    p = 3 * 2^43 - 1;
    B.<i,j,k> = QuaternionAlgebra(-1,-p);
    O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])

    L = 2
    e = 350;
    print(f"Generating g1...")
    foundsol = False
    while not foundsol:
        g1 = RandomReduced(O, sbound=0.5);
        s = ZZ(g1[0][0]);
        r = g1[0][1];
        t = ZZ(g1[1][1]);
        K = GF(s);
        ac = Compute_ac( O, g1 )
        if ac is None: continue;
        else:
            print(f"DONE!\n")
            print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
            print(f"           g1  = [{g1[0,0]}, {g1[0,1]}]")
            print(f"                 [{g1[1,0]}, {g1[1,1]}]")
            print(f"           a   = {a}")
            print(f"           c   = {c}")
            print(f"           e   = {e}")
            print(f"\nGenerating g2...");
            foundsol = True
    foundsol = False
    while not foundsol:
        g2 = RandomReduced(O, sbound=0.5);
        s = ZZ(g1[0][0]);
        r = g1[0][1];
        t = ZZ(g1[1][1]);
        K = GF(s);
        ac = Compute_ac( O, g2 )
        if ac is None: continue;
        else:
            print(f"DONE!\n\n")
            print(f"gcd(n(a),n(c)) = {gcd(a.reduced_norm(),c.reduced_norm())}")
            print(f"           g2  = [{g2[0,0]}, {g2[0,1]}]")
            print(f"                 [{g2[1,0]}, {g2[1,1]}]")
            print(f"           a   = {a}")
            print(f"           c   = {c}")
            print(f"           e   = {e}")
            foundsol = True

#GenerateCheatPairs()

def twoadic(N):
    M = N
    while (N % 2) == 0:
        N = N//2;
    return (M // N);

p = 3 * 2^43 - 1;
Le = 2^1000;
B.<i,j,k> = QuaternionAlgebra(-1,-p);
O = B.quaternion_order([1,i,i/2+j/2,1/2+k/2])
r = 183473551550885215239093361195645100494348123432123 + 103489357288336036570742620116745037321594270065985*i+5125*j+6646*k;
s = 5421341;
t = 8184799884502414753566140534100785841846441788692629529322485064161801624924877657812355294081;
R = ZZ(r.reduced_norm())
print(s*t-R);
F = GF(s);
u = F(t*p*R);
for ind in range(2000):
    c20 = F.random_element();
    if (Le * u^(-1) - c20^2).is_square():
        c10 = (Le * u^(-1) - c20^2).sqrt();
        c2 = ZZ(c20);
        c1 = ZZ(c10);
        if ((c1-c2) % 2) == 1:
            c  = c1 * r.conjugate() * j + c2 * r.conjugate() * k;
            A0 = ZZ(Le-t*c.reduced_norm());
            A1 = A0 // s;
            A  = A1 // twoadic(A1);
            if A.is_prime():
                a1, a2 = two_squares(A)

##  twoadic:= function(N)
##      M:=N;
##      while (N mod 2) eq 0 do
##      N:=N div 2;
##      end while;
##      return (M div N);
##  end function;
##  p:= 3 * 2^43 - 1;
##  Le:=2^1000;
##  B<i,j,k>:= QuaternionAlgebra<Rationals()|-1,-p>;
##  O:=QuaternionOrder([1,i,i/2+j/2,1/2+k/2]);
##  r:=183473551550885215239093361195645100494348123432123+
##  103489357288336036570742620116745037321594270065985*i+5125*j+6646*k;
##  s:=5421341;
##  t:=8184799884502414753566140534100785841846441788692629529322485064161801624924877\
##  657812355294081;
##  R:=Integers()!Norm(r);
##  s*t-Norm(r);
##  F:=FiniteField(s);
##  u:=F!t*p*R;
##  for i in [1..2000] do
##      c20:=Random(F);
##      if IsSquare(Le*(u^(-1))-c20^2) eq true then
##          c10:=Sqrt(Le*(u^(-1))-c20^2);
##          c2:=Integers()!c20;
##          c1:=Integers()!c10;
##          if ((c1-c2) mod 2) eq 1 then
##              c:=c1*Conjugate(r)*j+c2*Conjugate(r)*k;
##              A0:=Integers()!(Le-t*Norm(c));
##              A1:=A0 div s;
##              A:=A1 div twoadic(A1);
##              if IsPrime(A) then
##                  NormEquation(1,A);
##              end if;
##          end if;
##      end if;
##  end for;
