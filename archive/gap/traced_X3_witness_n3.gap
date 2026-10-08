LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetRecursionTrapInterval(50000);

# Hard-code n=3 traced run + X_3 cofactor extraction.
ReadParam := function(name, default_)
  if IsBound(GAPInfo.SystemEnvironment.(name))
     and GAPInfo.SystemEnvironment.(name) <> "" then
    return GAPInfo.SystemEnvironment.(name);
  fi;
  return default_;
end;

n := 3;
pStr := ReadParam("PSTR", "2^31-1");
q_int := Int(ReadParam("Q", "5"));
t_int := Int(ReadParam("T", "7"));

p := EvalString(pStr);
R := GF(p);
q := q_int * One(R);
t := t_int * One(R);

generator_names := Concatenation(
    List([1..n-1], i -> Concatenation("x", String(i))),
    List([1..n-1], i -> Concatenation("y", String(i)))
);
F := FreeAssociativeAlgebraWithOne(R, generator_names);
gens := GeneratorsOfAlgebra(F);
one := gens[1];
xs := gens{[2..n]};
ys := gens{[n+1..2*n-1]};

# Build relations (mirrors zipper.gap, untwist=bigelow).
build_relations := function(one_, xs_, ys_, q_, t_)
  local n_, inv, braid, zipper, untwisting, S, i, j, k,
        yProd, xProd, yProd_after, xProd_after, yProd_full, xProd_full,
        xProd_reverse, yProd_reverse;
  n_ := Length(xs_) + 1;
  inv := []; braid := []; zipper := []; untwisting := [];
  for i in [1..n_-1] do
    Add(inv, xs_[i]*ys_[i] - one_);
    Add(inv, ys_[i]*xs_[i] - one_);
  od;
  for i in [1..n_-2] do
    for j in [i+1..n_-1] do
      if j - i >= 2 then
        Add(braid, xs_[i]*xs_[j] - xs_[j]*xs_[i]);
        Add(braid, ys_[i]*ys_[j] - ys_[j]*ys_[i]);
        Add(braid, ys_[i]*xs_[j] - xs_[j]*ys_[i]);
        Add(braid, xs_[i]*ys_[j] - ys_[j]*xs_[i]);
      else
        Add(braid, xs_[i]*xs_[j]*xs_[i] - xs_[j]*xs_[i]*xs_[j]);
        Add(braid, ys_[i]*ys_[j]*ys_[i] - ys_[j]*ys_[i]*ys_[j]);
      fi;
    od;
  od;
  for i in [1..n_-2] do
    Add(braid, xs_[i]*xs_[i+1]*ys_[i]   - ys_[i+1]*xs_[i]*xs_[i+1]);
    Add(braid, ys_[i]*ys_[i+1]*xs_[i]   - xs_[i+1]*ys_[i]*ys_[i+1]);
    Add(braid, xs_[i]*ys_[i+1]*ys_[i]   - ys_[i+1]*ys_[i]*xs_[i+1]);
    Add(braid, ys_[i]*xs_[i+1]*xs_[i]   - xs_[i+1]*xs_[i]*ys_[i+1]);
  od;
  S := List([1..n_], i -> []);
  for j in [1..n_-1] do
    S[2][j] := q_*ys_[j] + (1-q_)*one_ - xs_[j];
  od;
  for i in [3..n_] do
    for j in [1..n_-i+1] do
      yProd := ys_[j];
      for k in [j+1..j+i-2] do yProd := yProd * ys_[k]; od;
      xProd := xs_[j];
      for k in [j+1..j+i-2] do xProd := xProd * xs_[k]; od;
      S[i][j] := (q_^(i-1)*yProd - xProd) * S[i-1][j];
    od;
  od;
  for j in [1..n_-2] do
    Add(zipper,
        (q_*ys_[j+1] + (1-q_)*one_ - xs_[j+1]
         - q_*ys_[j]*ys_[j+1] + xs_[j]*xs_[j+1]) * S[2][j]);
  od;
  for i in [3..n_-1] do
    for j in [1..n_-i] do
      yProd_after := ys_[j+1];
      for k in [j+2..j+i-1] do yProd_after := yProd_after * ys_[k]; od;
      yProd_full := ys_[j] * yProd_after;
      xProd_after := xs_[j+1];
      for k in [j+2..j+i-1] do xProd_after := xProd_after * xs_[k]; od;
      xProd_full := xs_[j] * xProd_after;
      Add(zipper,
          (q_^(i-1)*yProd_after - xProd_after
           - q_^(i-1)*yProd_full + xProd_full) * S[i][j]);
    od;
  od;
  for j in [1..n_-1] do
    Add(untwisting, (xs_[j] - t_*one_) * S[2][j]);
  od;
  for i in [3..n_] do
    for j in [1..n_-i+1] do
      xProd_reverse := xs_[j+i-2];
      yProd_reverse := ys_[j+i-2];
      for k in [2..i-1] do
        xProd_reverse := xProd_reverse * xs_[j+i-1-k];
        yProd_reverse := yProd_reverse * ys_[j+i-1-k];
      od;
      Add(untwisting, (xProd_reverse - q_^(i-1)*yProd_reverse) * S[i][j]);
    od;
  od;
  return rec(inv := inv, braid := braid, zipper := zipper,
             untwisting := untwisting, S := S);
end;

rels := build_relations(one, xs, ys, q, t);
inv_GB := GP2NPList(rels.inv);
zipper_simplified := Filtered(List(GP2NPList(rels.zipper),
                              pp -> StrongNormalFormNP(pp, inv_GB)),
                              pp -> pp <> []);
untwisting_simplified := Filtered(List(GP2NPList(rels.untwisting),
                                  pp -> StrongNormalFormNP(pp, inv_GB)),
                                  pp -> pp <> []);
I := Concatenation(inv_GB,
                   GP2NPList(rels.braid),
                   zipper_simplified,
                   untwisting_simplified);
Print("Input relations: ", Length(I), "\n");

Print("Running SGrobnerTrace at n=3 ...\n");
GBT := SGrobnerTrace(I);
B := List(GBT, r -> r.pol);
Print("Traced GB size: ", Length(B), "\n");

# X_3 = S[3][1]
X3 := rels.S[3][1];
X3np := GP2NPList([X3])[1];
Print("X_3 NP term count = ", Length(X3np[1]), "\n");

# StrongNormalFormTraceDiff: returns rec(pol := nf, trace := cofactor info)
# expressing X3 - sum(coeffs * basis * coeffs) = nf.
Print("Running StrongNormalFormTraceDiff ...\n");
nfd := StrongNormalFormTraceDiff(X3np, GBT);
Print("Reduction record fields = ", RecNames(nfd), "\n");

isZero := (nfd.pol = []) or
          (IsList(nfd.pol) and Length(nfd.pol) >= 1 and nfd.pol[1] = []);
if isZero then
  Print("RESULT n=3: NF(X_3, B) = 0  --> X_3 lies in ideal at n=3.\n");
  Print("Cofactor trace length = ", Length(nfd.trace), "\n");
  # The trace expresses X3 = sum_k coef_k * left_k * I[idx_k] * right_k, in
  # terms of the *input* relations I. That's an explicit witness from GBNP.
else
  Print("RESULT n=3: NF != 0; remainder term count = ", Length(nfd.pol[1]), "\n");
fi;

QUIT_GAP();
