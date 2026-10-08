"""
zip_lib.py -- exact arithmetic in the braid group algebra R B_n, R = Z[q^{+-1}] (coeffs as Fractions),
for auditing the zipping algorithm of thesis.tex (section:algorithm) and the diagrammatic identities.

Braid words: tuples of nonzero ints; +i = sigma_i, -i = sigma_i^{-1}.
Multiplication convention (thesis): composition; in a product a*b the RIGHT factor b is applied first.
Word concatenation a+b corresponds to the product a*b.

Canonical form: image tuple of the free-group generators x_1..x_n under the Artin action
(sigma_k: x_k -> x_k x_{k+1} x_k^{-1}, x_{k+1} -> x_k), composed letter by letter.
Faithfulness of the Artin representation B_n -> Aut(F_n) is a standard theorem (Artin), used as a black box.
Whether the map word -> automorphism is a homomorphism or an antihomomorphism is irrelevant:
it respects concatenation and is injective on B_n, so two words are equal in B_n iff their keys agree.

Coefficients: Laurent polynomials in q with Fraction coefficients, stored as dict {exponent: Fraction}.
"""
from fractions import Fraction


# ---------- Laurent polynomials in q ----------
class LP:
    __slots__ = ('d',)

    def __init__(self, d=None):
        if d is None:
            self.d = {}
        elif isinstance(d, dict):
            self.d = {e: Fraction(c) for e, c in d.items() if c != 0}
        else:  # scalar
            c = Fraction(d)
            self.d = {0: c} if c != 0 else {}

    @staticmethod
    def q(e=1, c=1):
        return LP({e: Fraction(c)})

    def __add__(self, o):
        o = o if isinstance(o, LP) else LP(o)
        d = dict(self.d)
        for e, c in o.d.items():
            d[e] = d.get(e, 0) + c
        return LP(d)
    __radd__ = __add__

    def __neg__(self):
        return LP({e: -c for e, c in self.d.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, LP) else LP(o)))

    def __rsub__(self, o):
        return LP(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, LP) else LP(o)
        d = {}
        for e1, c1 in self.d.items():
            for e2, c2 in o.d.items():
                d[e1 + e2] = d.get(e1 + e2, 0) + c1 * c2
        return LP(d)
    __rmul__ = __mul__

    def iszero(self):
        return not self.d

    def __eq__(self, o):
        o = o if isinstance(o, LP) else LP(o)
        return self.d == o.d

    def __hash__(self):
        return hash(tuple(sorted(self.d.items())))

    def __repr__(self):
        if not self.d:
            return '0'
        parts = []
        for e in sorted(self.d):
            c = self.d[e]
            cs = str(c) if c.denominator != 1 else str(c.numerator)
            if e == 0:
                parts.append(cs)
            elif e == 1:
                parts.append(cs + '*q')
            else:
                parts.append(cs + '*q^' + str(e))
        return ' + '.join(parts)


ONE = LP(1)
ZERO = LP()


def qp(e, c=1):
    return LP.q(e, c)


# ---------- free group words ----------
def freduce(w):
    out = []
    for a in w:
        if out and out[-1] == -a:
            out.pop()
        else:
            out.append(a)
    return tuple(out)


def finv(w):
    return tuple(-a for a in reversed(w))


def subst(w, images):
    """substitute x_j -> images[j-1] in free word w."""
    out = []
    for a in w:
        img = images[a - 1] if a > 0 else finv(images[-a - 1])
        for b in img:
            if out and out[-1] == -b:
                out.pop()
            else:
                out.append(b)
    return tuple(out)


# ---------- Artin action ----------
_ID_CACHE = {}


def identity_key(n):
    if n not in _ID_CACHE:
        _ID_CACHE[n] = tuple((i,) for i in range(1, n + 1))
    return _ID_CACHE[n]


def apply_letter(key, a):
    """key = tuple of images; apply sigma_a (a>0) or sigma_{-a}^{-1} (a<0) as in verify_central.py."""
    w = list(key)
    k = abs(a)
    i = k - 1
    A, B = w[i], w[i + 1]
    if a > 0:
        w[i] = freduce(A + B + finv(A))
        w[i + 1] = A
    else:
        w[i] = B
        w[i + 1] = freduce(finv(B) + A + B)
    return tuple(w)


def key_of_word(word, n):
    key = identity_key(n)
    for a in word:
        assert 1 <= abs(a) <= n - 1, 'letter %d out of range for n=%d' % (a, n)
        key = apply_letter(key, a)
    return key


def compose_keys(A, B):
    """key of (word_a + word_b) given keys A = key(word_a), B = key(word_b)."""
    return tuple(subst(b, A) for b in B)


# ---------- elements of R B_n ----------
class Elt:
    """dict key -> [word, coeff]; word is a representative braid word."""
    __slots__ = ('n', 'terms')

    def __init__(self, n, terms=None):
        self.n = n
        self.terms = {} if terms is None else terms

    @staticmethod
    def word(n, w, c=ONE):
        w = tuple(w)
        c = c if isinstance(c, LP) else LP(c)
        if c.iszero():
            return Elt(n)
        return Elt(n, {key_of_word(w, n): [w, c]})

    @staticmethod
    def one(n):
        return Elt.word(n, ())

    @staticmethod
    def gen(n, i, sign=1):
        return Elt.word(n, (i if sign > 0 else -i,))

    def copy(self):
        return Elt(self.n, {k: [w, c] for k, (w, c) in self.terms.items()})

    def __add__(self, o):
        if not isinstance(o, Elt):
            o = Elt.one(self.n) * o
        assert self.n == o.n
        t = {k: [w, c] for k, (w, c) in self.terms.items()}
        for k, (w, c) in o.terms.items():
            if k in t:
                nc = t[k][1] + c
                if nc.iszero():
                    del t[k]
                else:
                    t[k][1] = nc
            else:
                t[k] = [w, c]
        return Elt(self.n, t)
    __radd__ = __add__

    def __neg__(self):
        return Elt(self.n, {k: [w, -c] for k, (w, c) in self.terms.items()})

    def __sub__(self, o):
        if not isinstance(o, Elt):
            o = Elt.one(self.n) * o
        return self + (-o)

    def __rsub__(self, o):
        return (Elt.one(self.n) * o) - self

    def scale(self, c):
        c = c if isinstance(c, LP) else LP(c)
        if c.iszero():
            return Elt(self.n)
        return Elt(self.n, {k: [w, cc * c] for k, (w, cc) in self.terms.items()})

    def __mul__(self, o):
        if not isinstance(o, Elt):
            return self.scale(o)
        assert self.n == o.n
        t = {}
        for k1, (w1, c1) in self.terms.items():
            for k2, (w2, c2) in o.terms.items():
                k = compose_keys(k1, k2)
                c = c1 * c2
                if k in t:
                    nc = t[k][1] + c
                    if nc.iszero():
                        del t[k]
                    else:
                        t[k][1] = nc
                else:
                    t[k] = [w1 + w2, c]
        return Elt(self.n, t)

    def __rmul__(self, o):
        return self.scale(o)

    def iszero(self):
        return not self.terms

    def __eq__(self, o):
        return (self - o).iszero()

    def nterms(self):
        return len(self.terms)

    def inv(self):
        """inverse of a single group element (monomial with coefficient 1)."""
        assert len(self.terms) == 1
        (k, (w, c)), = self.terms.items()
        assert c == ONE
        return Elt.word(self.n, finv(w))

    def __repr__(self):
        if not self.terms:
            return '0'
        items = sorted(self.terms.values(), key=lambda wc: (len(wc[0]), wc[0]))
        return ' + '.join('(%s)*%s' % (c, wstr(w)) for w, c in items)


def wstr(w):
    if not w:
        return '1'
    return ''.join((('s%d' % a) if a > 0 else ('S%d' % -a)) for a in w)


# ---------- braid words of the thesis ----------
def w_v(i, j):      # varsigma_{i,j} = s_j ... s_{i+j-1}
    return tuple(range(j, i + j))


def w_vb(i, j):     # varsigma-bar_{i,j} = s_j^{-1} ... s_{i+j-1}^{-1}
    return tuple(-m for m in range(j, i + j))


def w_vinv(i, j):   # varsigma_{i,j}^{-1} = s_{i+j-1}^{-1} ... s_j^{-1}
    return finv(w_v(i, j))


def w_vbinv(i, j):  # varsigma-bar_{i,j}^{-1} = s_{i+j-1} ... s_j
    return finv(w_vb(i, j))


def w_vk(i, j, k, order='asc'):
    """varsigma_{i,j,k}: thesis prints varsigma_{i,j} varsigma_{i,j+1} ... varsigma_{i,j+k-1} ('asc').
    'desc' = varsigma_{i,j+k-1} ... varsigma_{i,j}."""
    fs = [w_v(i, j + m) for m in range(k)]
    if order == 'desc':
        fs = fs[::-1]
    return tuple(a for f in fs for a in f)


def w_vbk(i, j, k, order='asc'):
    fs = [w_vb(i, j + m) for m in range(k)]
    if order == 'desc':
        fs = fs[::-1]
    return tuple(a for f in fs for a in f)


def w_Delta(i, j):
    """half twist on strands j..j+i-1: Delta_{i,j} = varsigma_{i-1,j} Delta_{i-1,j} (thesis mirrored-slider def.)"""
    if i <= 1:
        return ()
    return w_v(i - 1, j) + w_Delta(i - 1, j)


def w_Theta(i, j):  # lasso full twist Theta_{i,j} = varsigma_{i,j} varsigma-bar_{i,j}^{-1}
    return w_v(i, j) + w_vbinv(i, j)


def w_ThetaK(k):    # Theta_2 = s_1; Theta_k = Theta_{k-1,1}
    return (1,) if k == 2 else w_Theta(k - 1, 1)


def E(n, w, c=ONE):
    return Elt.word(n, w, c)


def W(n, *words):
    """product of braid words (each a tuple) as an Elt."""
    w = tuple(a for ww in words for a in ww)
    return E(n, w)


_XCACHE = {}


def X(n, i, j):
    """slider X_{i,j} (Definition def:shifted_sliders)."""
    assert i >= 2 and j >= 1 and i + j - 1 <= n, (n, i, j)
    key = (n, i, j)
    if key in _XCACHE:
        return _XCACHE[key]
    x = E(n, (-j,), qp(1)) + E(n, (), qp(0, 1) - qp(1)) - E(n, (j,))
    for m in range(3, i + 1):
        x = (E(n, w_vb(m - 1, j), qp(m - 1)) - E(n, w_v(m - 1, j))) * x
    _XCACHE[key] = x
    return x


_MXCACHE = {}


def mX(n, i, j):
    """mirrored slider _{i,j}X by the recursive definition:
    _{2,j}X = X_{2,j}; _{i,j}X = (q^{i-1} varsigma_{i-1,j}^{-1} - varsigma-bar_{i-1,j}^{-1}) _{i-1,j+1}X."""
    key = (n, i, j)
    if key in _MXCACHE:
        return _MXCACHE[key]
    if i == 2:
        r = X(n, 2, j)
    else:
        r = (E(n, w_vinv(i - 1, j), qp(i - 1)) - E(n, w_vbinv(i - 1, j))) * mX(n, i - 1, j + 1)
    _MXCACHE[key] = r
    return r


def c_(n, i, j):
    """c_{i,j} = q^i vbar_{i,j} - v_{i,j}   (so X_{i,j} = c_{i-1,j} X_{i-1,j})."""
    return E(n, w_vb(i, j), qp(i)) - E(n, w_v(i, j))


def zip_relator(n, i, j):
    """Shifted zipper relator Z_i at position j (lem:zipshifted), as LHS - RHS:
    i=2:  X_{2,j+1} X_{2,j} - (q vbar_{2,j} - v_{2,j}) X_{2,j}
    i>=3: (q^{i-1} vbar_{i-1,j+1} - v_{i-1,j+1}) X_{i,j} - (q^{i-1} vbar_{i,j} - v_{i,j}) X_{i,j}."""
    if i == 2:
        rhs = (E(n, w_vb(2, j), qp(1)) - E(n, w_v(2, j))) * X(n, 2, j)
        return X(n, 2, j + 1) * X(n, 2, j) - rhs
    lhs = (E(n, w_vb(i - 1, j + 1), qp(i - 1)) - E(n, w_v(i - 1, j + 1))) * X(n, i, j)
    rhs = (E(n, w_vb(i, j), qp(i - 1)) - E(n, w_v(i, j))) * X(n, i, j)
    return lhs - rhs


def F(n, i, j, k, order='asc'):
    """katamari friends F_{i,j,k} = varsigma_{i-1,j+1,k} (thesis); order 'asc' as printed, 'desc' reversed."""
    return E(n, w_vk(i - 1, j + 1, k, order))


def phi_word(w, n):
    """antiautomorphism phi: s_i -> s_{n-i}, reversing word order."""
    return tuple((n - a) if a > 0 else -(n + a) for a in reversed(w))


def phi(e):
    n = e.n
    out = Elt(n)
    for k, (w, c) in e.terms.items():
        out = out + Elt.word(n, phi_word(w, n), c)
    return out
