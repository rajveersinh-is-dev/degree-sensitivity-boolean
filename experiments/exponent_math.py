from decimal import Decimal, getcontext

getcontext().prec = 80
C_CRIT = Decimal(6).ln() / Decimal(3).ln()

def pow_gt(n, d):
    """exact comparison  n > d^c  via high-precision decimal"""
    return Decimal(n) > Decimal(d) ** C_CRIT

def nmin(d):
    """The least n such that n > d^c"""
    x = Decimal(d) ** C_CRIT
    n = int(x) + 1
    while Decimal(n) <= x:
        n += 1
    while Decimal(n - 1) > x:
        n -= 1
    return n, x
