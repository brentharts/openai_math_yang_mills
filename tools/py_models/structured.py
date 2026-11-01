"""Nested loops, several loops in a row, break and continue (see tools/pyast2lean.py).

Every loop has its own invariant: the asserts right before it. For nested loops the inner invariant is stated
in closed form (it names `i`, `j`, `s` directly): a fact about a variable the inner loop changes is not carried
into the inner loop's conditions, a fact about the other variables is. Asserts after a loop and before the next
statement are that loop's exit condition; with `break` they must hold on every way out."""


def nested(n: Nat) -> int:
    assert n >= 0
    s = 0
    i = 0
    assert s == i * n and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        j = 0
        assert s == i * n + j and 0 <= j and j <= n and i < n
        decreases(n - j)
        while j < n:
            s += 1
            j += 1
        i += 1
    assert s == n * n
    return s


def two_phase(n: Nat) -> int:
    assert n >= 0
    a = 0
    i = 0
    assert a == i and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        a += 1
        i += 1
    assert a == n
    b = 0
    k = 0
    assert b == 2 * k and 0 <= k and k <= a
    decreases(a - k)
    while k < a:
        b += 2
        k += 1
    assert b == 2 * n
    return b


def first_neg(a: list[int]) -> int:
    k = 0
    r = -1
    assert 0 <= k and k <= len(a) and r == -1
    decreases(len(a) - k)
    while k < len(a):
        if a[k] < 0:
            r = k
            break
        k += 1
    assert r == -1 or (0 <= r and r < len(a))
    return r


def count_odd(n: Nat) -> int:
    assert n >= 0
    c = 0
    i = 0
    assert 0 <= c and c <= i and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        i += 1
        if i % 2 == 0:
            continue
        c += 1
    assert c <= n
    return c


def tri(n: Nat) -> int:
    assert n >= 0
    t = 0
    assert 2 * t == i * (i - 1) and 0 <= i and i <= n
    for i in range(0, n):
        t += i
    assert 2 * t == n * (n - 1)
    return t
