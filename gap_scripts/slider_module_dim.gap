# Compute the dimension of the SLIDER MODULE M_p = R[B_3] / (R[B_3] p(x2) + R[B_3] z)
# (cyclic left B_3-module) via GBNP's module Grobner basis (SGrobnerModule / BaseQM).
# Validation target: linear p = x2 - t  must give dim 5.
LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetInfoLevel(InfoGBNPTime, 1);
SetRecursionTrapInterval(200000);

p_prime := 2^31-1;
R := GF(p_prime);
q := 5*One(R);
t := 7*One(R);

A := FreeAssociativeAlgebraWithOne(R, "x1","x2","y1","y2");
g := GeneratorsOfAlgebra(A);
e := g[1]; x1:=g[2]; x2:=g[3]; y1:=g[4]; y2:=g[5];

# two-sided group B_3 relations (matches zipper.gap, n=3)
twosid := [
  x1*y1 - e, y1*x1 - e, x2*y2 - e, y2*x2 - e,
  x1*x2*x1 - x2*x1*x2,
  y1*y2*y1 - y2*y1*y2,
  x1*x2*y1 - y2*x1*x2,
  y1*y2*x1 - x2*y1*y2,
  x1*y2*y1 - y2*y1*x2,
  y1*x2*x1 - x2*x1*y2
];

# validated zipper element z:  Theta + q x1 + (1-q) x2 x1 - x2 x2 x1 - q,  Theta = x1 x2 x2 x1
zpoly := x1*x2*x2*x1 + q*x1 + (1-q)*x2*x1 - x2*x2*x1 - q*e;

RunCase := function(name, ppoly, maxno)
  local pref, GBR, B, t0, dimstr;
  Print("\n================ CASE: ", name, " ================\n");
  pref := [ ppoly, zpoly ];   # module relations: p(x2) u = 0 and z u = 0
  t0 := Runtime();
  GBR := SGrobnerModule(GP2NPList(pref), GP2NPList(twosid));
  Print("[", name, "] SGrobnerModule done in ", Runtime()-t0, " ms; #ts=",
        Length(GBR.ts), " #prefix=", Length(GBR.p), "\n");
  B := BaseQM(GBR, 4, 1, maxno);
  if maxno > 0 and Length(B) >= maxno then
    dimstr := Concatenation(">= ", String(maxno), " (cap hit)");
  else
    dimstr := String(Length(B));
  fi;
  Print("[", name, "] MODULE DIMENSION = ", dimstr, "\n");
  return Length(B);
end;

# 1) validation: linear, full basis (maxno=0)
RunCase("linear x2 - t", x2 - t*e, 0);

QUIT_GAP();
