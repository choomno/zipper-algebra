# EXACT over Q(q,t): confirm min poly of Theta_3 on X_3 = p*(lambda) in
# Z_3/(zipper + two-strand p_2), for p_2 in {linear, (sigma-t)^2, sigma^2-t}.
# n=3, FunctionField(Rationals,[q,t]); validate presentation via dim.
LoadPackage("gbnp");
SetInfoLevel(InfoGBNP,0); SetInfoLevel(InfoGBNPTime,0);
F := FunctionField(Rationals, ["q","t"]);
ind := IndeterminatesOfFunctionField(F); q := ind[1]; t := ind[2];
A := FreeAssociativeAlgebraWithOne(F, "x1","x2","y1","y2");
gg := GeneratorsOfAlgebra(A); e:=gg[1]; x1:=gg[2]; x2:=gg[3]; y1:=gg[4]; y2:=gg[5];
xs:=[x1,x2]; ys:=[y1,y2];

baseRels := function()
  local r;
  r := [];
  Add(r, x1*y1-e); Add(r, y1*x1-e); Add(r, x2*y2-e); Add(r, y2*x2-e);
  Add(r, x1*x2*x1-x2*x1*x2); Add(r, y1*y2*y1-y2*y1*y2);
  Add(r, x1*x2*y1-y2*x1*x2); Add(r, y1*y2*x1-x2*y1*y2);
  Add(r, x1*y2*y1-y2*y1*x2); Add(r, y1*x2*x1-x2*x1*y2);
  return r;
end;

S2 := [ q*y1+(1-q)*e-x1, q*y2+(1-q)*e-x2 ];
S31 := (q^2*y1*y2 - x1*x2)*S2[1];
Theta3 := x1*x2*x2*x1;
zipperZ2 := (q*y2+(1-q)*e-x2 - q*y1*y2 + x1*x2)*S2[1];

NFof := function(GB, elt)
  return StrongNormalFormNP(GP2NPList([elt])[1], GB);
end;
IsZ := function(nf) return nf=[] or (IsList(nf) and Length(nf)>=1 and nf[1]=[]); end;

runcase := function(name, untwPolys, checks)
  local rels, j, GB, d, c;
  rels := baseRels();
  Add(rels, zipperZ2);
  for j in [1,2] do Add(rels, untwPolys[j]); od;
  GB := SGrobner(GP2NPList(rels));
  d := DimQA(GB, 4);
  Print("=== ", name, " ===\n");
  Print("  dim = ", d, "   |GB| = ", Length(GB), "\n");
  Print("  X_3 = 0? ", IsZ(NFof(GB, S31)), "\n");
  for c in checks do
    Print("  ", c[1], " reduces to 0? ", IsZ(NFof(GB, c[2])), "\n");
  od;
end;

# linear p_2 = x - t; p* = lam + t
runcase("linear (sigma-t)", [ (x1-t*e)*S2[1], (x2-t*e)*S2[2] ],
  [ ["(Theta_3 + t) X_3   [=p*]", Theta3*S31 + t*S31] ]);

# square-zero p_2 = (x-t)^2; p* = (lam+t)^2
runcase("square-zero (sigma-t)^2", [ (x1*x1-2*t*x1+t^2*e)*S2[1], (x2*x2-2*t*x2+t^2*e)*S2[2] ],
  [ ["(Theta_3 + t) X_3        [proper divisor, expect NONzero]", Theta3*S31 + t*S31],
    ["(Theta_3+t)^2 X_3   [=p*]", Theta3*(Theta3*S31) + 2*t*(Theta3*S31) + t^2*S31] ]);

# framing p_2 = x^2 - t; p* = lam^2 - t
runcase("framing (sigma^2-t)", [ (x1*x1-t*e)*S2[1], (x2*x2-t*e)*S2[2] ],
  [ ["(Theta_3 + t) X_3        [expect NONzero]", Theta3*S31 + t*S31],
    ["(Theta_3^2 - t) X_3 [=p*]", Theta3*(Theta3*S31) - t*S31] ]);

QUIT_GAP();
