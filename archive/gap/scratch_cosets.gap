# Explore B_3 / <s2> coset structure via a finite quotient or free reduction.
# B_3 = < s1, s2 | s1 s2 s1 = s2 s1 s2 >.
F := FreeGroup("a","b");;
a := F.1;; b := F.2;;
B3 := F / [ a*b*a*(b*a*b)^-1 ];;
s1 := B3.1;; s2 := B3.2;;
Print("Order of B3 (expect infinity / fail): \n");
