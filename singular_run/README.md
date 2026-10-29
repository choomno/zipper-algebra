# Independent Singular cross-check

`x3-membership-Q.sing` performs a separate specialized calculation with
[Singular](https://www.singular.uni-kl.de/)'s `freegb.lib` and Letterplace
implementation. It uses the four-strand presentation with 31 defining relations
in `x1,x2,x3,y1,y2,y3`, specialized to `q=5`, `t=7` over the rationals.
The target `X_3` is not adjoined as a relation.

From this directory, with Singular installed, run:

```sh
Singular -q x3-membership-Q.sing
```

The script computes `J = twostd(I)` with Letterplace degree bound 12, then prints
`NF(X3, J)`. The expected remainder is `0`. The degree bound and specialization
are explicit parts of this calculation; the script does not assert that an
incomplete computation supplies a completed Gröbner basis.

This is a specialized membership cross-check. It does not read or verify a
supplied cofactor representation, and it does not establish the generic identity
over `Q(q,t)`. The [certificate reproduction guide](../docs/X3_REPRODUCTION.md)
describes production and exact verification of that generic identity.

The script is covered by the repository's [MIT license](../LICENSE). Singular
and its libraries retain their own licenses. This independent solver run is
separate from the default Python test suite.
