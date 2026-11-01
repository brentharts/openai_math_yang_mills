"""Loops with their invariants, written as ordinary Python (see tools/pyast2lean.py).
Asserts at the top are the precondition; asserts right before a loop are its invariant;
asserts after the loop are the postcondition; decreases(...) names the variant."""


def sum_to(n: Nat) -> int:
    assert n >= 0
    s = 0
    i = 0
    assert 2 * s == i * (i - 1) and 0 <= i and i <= n
    decreases(n - i)
    while i < n:
        s = s + i
        i = i + 1
    assert 2 * s == n * (n - 1)
    return s


def fib(n: Nat) -> int:
    assert n >= 0
    a, b = 0, 1
    k = 0
    assert 0 <= a and a <= b and k <= n
    decreases(n - k)
    while k < n:
        a, b = b, a + b
        k += 1
    assert 0 <= a
    return a


def power(x: int, n: Nat) -> int:
    assert n >= 0
    r = 1
    e = 0
    assert r == x ** e and 0 <= e and e <= n
    decreases(n - e)
    while e < n:
        r = r * x
        e += 1
    assert r == x ** n
    return r


def gcd(a: Nat, b: Nat) -> int:
    assert a >= 0 and b >= 0
    x, y = a, b
    assert x >= 0 and y >= 0
    decreases(y)
    while y > 0:
        x, y = y, x % y
    assert x >= 0
    return x
