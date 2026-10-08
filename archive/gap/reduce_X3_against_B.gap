LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetInfoLevel(InfoGBNPTime, 1);
SetRecursionTrapInterval(50000);

# Parameters (override via env). Defaults: n=4, q=5, t=7.
ReadParam := function(name, default_)
  if IsBound(GAPInfo.SystemEnvironment.(name))
     and GAPInfo.SystemEnvironment.(name) <> "" then
    return GAPInfo.SystemEnvironment.(name);
  fi;
  return default_;
end;

n := Int(ReadParam("N", "4"));
pStr := ReadParam("PSTR", "2^31-1");
q_int := Int(ReadParam("Q", "5"));
t_int := Int(ReadParam("T", "7"));

p := EvalString(pStr);
R := GF(p);
q := q_int * One(R);
t := t_int * One(R);

# Load the saved (untraced) Gröbner basis B for this (n, p, q, t).
gbFile := Concatenation("logs_and_traces/gap_traces/grobner_zipper_bigelow-untraced-",
                        String(n), "-(", pStr, ")-",
                        String(q_int), "-", String(t_int), ".gap");
Print("Loading GB from: ", gbFile, "\n");
Read(gbFile);
Print("Loaded basis B with ", Length(B), " elements.\n");

# Rebuild the free algebra and compute X_3 = S[3][1] symbolically over GF(p).
generator_names := Concatenation(
    List([1..n-1], i -> Concatenation("x", String(i))),
    List([1..n-1], i -> Concatenation("y", String(i)))
);
F := FreeAssociativeAlgebraWithOne(R, generator_names);
generators := GeneratorsOfAlgebra(F);
one := generators[1];
xs := generators{[2..n]};
ys := generators{[n+1..2*n-1]};

# Construct sliders S[i][j] just as zipper.gap does.
S := List([1..n], i -> []);
for j in [1..n-1] do
  S[2][j] := q*ys[j] + (1-q)*one - xs[j];
od;
for i in [3..n] do
  for j in [1..n-i+1] do
    yProd := ys[j];
    for k in [j+1..j+i-2] do yProd := yProd * ys[k]; od;
    xProd := xs[j];
    for k in [j+1..j+i-2] do xProd := xProd * xs[k]; od;
    S[i][j] := (q^(i-1)*yProd - xProd) * S[i-1][j];
  od;
od;

X3 := S[3][1];
Print("X_3 = S[3][1] constructed.\n");

# Convert X_3 to NP form and reduce against B.
X3np := GP2NPList([X3])[1];
Print("X_3 as NP: leading monomial = ", X3np[1][1],
      ", number of terms = ", Length(X3np[1]), "\n");

Print("Reducing X_3 against B via StrongNormalFormNP ...\n");
nf := StrongNormalFormNP(X3np, B);

Print("Raw NF result = ", nf, "\n");
# An NP polynomial is zero iff its monomial list (nf[1]) is empty (or nf=[]).
isZero := (nf = []) or (IsList(nf) and Length(nf) >= 1 and nf[1] = []);
if isZero then
  Print("RESULT: StrongNormalFormNP(X_3, B) = 0  (X_3 IS in the ideal)\n");
  Print("VERDICT: GBNP independently confirms X_3 = 0 at n=", n,
        ", (q,t)=(", q_int, ",", t_int, "), p=", pStr, ".\n");
  # Cross-check via NP->GAP polynomial: should equal Zero(F)
  Print("Cross-check: NP2GP(nf, F) = Zero(F) ? ",
        NP2GP(nf, F) = Zero(F), "\n");
else
  Print("RESULT: StrongNormalFormNP(X_3, B) != 0 (non-zero remainder)\n");
  Print("Remainder term count = ", Length(nf[1]), "\n");
fi;

QUIT_GAP();
