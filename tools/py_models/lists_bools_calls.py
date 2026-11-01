"""bool variables, lists and calls between the functions of one file (see tools/pyast2lean.py).

bool is Lean `Bool` (a comparison is `decide`d when stored in one); list[int|float|bool] is `Array`. A list index
or a call to a function with a precondition gives a `_safe` theorem: 0 <= i < len(a), resp. the callee's
precondition at that call. Lean arrays are values, so a list argument is never changed in place (`b = list(a)`),
and `for x in a` has the hidden index `idx_x`, which the asserts before the loop may mention."""


def is_small(x: int) -> bool:
    return x < 10


def sq(n: Nat) -> int:
    assert n >= 0
    return n * n


def sum_sq(n: Nat) -> int:
    assert n >= 0
    t = 0
    i = 0
    assert 0 <= t and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        t += sq(i)
        i += 1
    assert t >= 0
    return t


def count_small(n: Nat) -> int:
    assert n >= 0
    t = 0
    i = 0
    assert 0 <= t and t <= i and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        if is_small(i):
            t += 1
        i += 1
    assert t <= n
    return t


def seen_big(n: Nat) -> bool:
    assert n >= 0
    seen = False
    i = 0
    assert (seen == False or i > 3) and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        if i > 3:
            seen = True
        i += 1
    assert (not seen) or n > 3
    return seen


def pick(flag: bool, x: int, y: int) -> int:
    if flag:
        return x
    return y


def fill_squares(n: Nat) -> list[int]:
    assert n >= 0
    a = [0] * n
    i = 0
    assert 0 <= i and i <= n and len(a) == n
    decreases(n - i)
    while i < n:
        a[i] = i * i
        i += 1
    assert len(a) == n
    return a


def evens(n: Nat) -> list[int]:
    assert n >= 0
    r: list[int] = []
    i = 0
    assert 2 * len(r) <= i + 1 and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        if i % 2 == 0:
            r.append(i)
        i += 1
    return r


def swap_first(a: list[int]) -> list[int]:
    assert len(a) >= 2
    b = list(a)
    b[0], b[1] = b[1], b[0]
    assert b[0] == a[1] and b[1] == a[0]
    return b


def total(a: list[int]) -> int:
    s = 0
    assert 0 <= idx_x and idx_x <= len(a)
    for x in a:
        s += x
    return s
