from flint import arb, ctx
from fractions import Fraction
print("prec", ctx.prec)
x = arb(0.1); print("arb(0.1) =", x.str(30, radius=True), "rad", x.rad(), "exact?", x.rad()==0, Fraction(0.1) == Fraction(x.mid().str(50,radius=False)) if False else "")
print("arb(0.1).mid() man_exp:", x.mid().man_exp(), Fraction(0.1).as_integer_ratio())
m,e = x.mid().man_exp(); print("matches float exactly:", Fraction(int(m))*Fraction(2)**int(e) == Fraction(0.1))
b = arb(1, 1)   # [0,2]
print("ball [0,2] > 0 :", b > 0, " >= 0:", b >= 0, " < 3:", b < 3)
c = arb(0, 1)
print("arb(0,1):", c, " lower", c.lower(), " upper", c.upper(), " >0", c>0, " >=0", c>=0, "<0", c<0)
d = arb(1, 0.5); print("[0.5,1.5] > 0:", d > 0, "lower:", d.lower(), type(d.lower()), "lower>=0", d.lower() >= 0)
e = arb(0, 0); print("exact 0: >=0", e >= 0, ">0", e > 0)
f = arb(-1e-300, 1e-300); print("[-2e-300,0] >= 0:", f>=0, " lower>=0", f.lower() >= 0, "lower", f.lower())
g = arb(3, 2**-60) ; print("rad of arb(3, 2^-60):", g.rad(), type(g.rad()), "rad == 2^-60?", g.rad() >= arb(2)**-60)
# rad() returns an upper bound for the radius: test with rounding
h = arb(1)/3; print("1/3:", h.str(radius=True), "rad", h.rad(), "mid", h.mid(), "contains 1/3 check: |mid-1/3|<=rad:", abs(Fraction(*h.mid().man_exp()) if False else 0))
mm, ee = h.mid().man_exp(); rm, re_ = h.rad().man_exp()
mid = Fraction(int(mm))*Fraction(2)**int(ee); rad = Fraction(int(rm))*Fraction(2)**int(re_)
print("|mid-1/3| <= rad:", abs(mid-Fraction(1,3)) <= rad)
# nan / inf
try:
    n = arb(float('nan')); print("arb(nan):", n, "n>0", n>0, "n>=0", n>=0, "n.lower()>=0", n.lower() >= 0)
except Exception as ex: print("nan exc", ex)
try:
    i = arb(float('inf')); print("arb(inf):", i, "i>0", i>0, "lower", i.lower(), "lower>=0", i.lower()>=0)
except Exception as ex: print("inf exc", ex)
# big int exactness
bi = 3**80 + 1; a = arb(bi); print("arb(3^80+1) rad:", a.rad(), " exact?", a.rad()==0)
print("arb * big int:", (arb(1)*bi).rad())
# mid of a ball with nonexact mid
print("mid().rad():", h.mid().rad())
# upper().sqrt().upper()
r2 = arb(2); s = r2.upper().sqrt().upper(); print("sqrt2 upper:", s, s.rad(), s >= arb(2).sqrt())
sm, se = s.man_exp(); print("upper^2 >= 2 exactly:", (Fraction(int(sm))*Fraction(2)**int(se))**2 >= 2)
print("comparisons of overlapping balls:", arb(1,1) > arb(0.5), arb(1,1) < arb(0.5), arb(1,1)==arb(1,1))
print("bool(arb comparison) type:", type(arb(1) > 0))
