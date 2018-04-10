import random as rd
from sage.all import *

def check_point_on_curve(cur, x, y):
    return (x**2 + y**2) % cur['p'] == (1 + cur['d'] * x**2 * y**2) % cur['p']

def inv_mod(a ,n):
    return pow(a, n-2, n)

def get_isomorphic_curve(cur):
    d = cur['d']
    p = cur['p']
    return ( ((-(1 + 14*d + d**2 ) % p) * inv_mod(48, p)) % p, ((-(1 - 33*d - 33 * d**2 + d**3) % p) * inv_mod(864, p)) % p )

def check_mov(cur):
    n = cur['n']
    p = cur['p']

    for k in range(2, 33):
        if (p**k - 1) % n == 0:
            return False
    return True

def jakobi_symbol(n, p):
    return pow(n, (p-1)/2, p)

K = 10


def extended_euclid(a, b):
    if (b == 0):
        return a, 1, 0
    d, x, y = extended_euclid(b, a % b)
    return d, y, x - a // b * y


def pre_division_test(n):
    if n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0 or n % 11 == 0 or n % 13 == 0:
        return False
    else:
        return True


def miller_rabine_test(n):
    d = n - 1
    s = 0
    k = 0
    while d % 2 == 0:
        d = d // 2
        s += 1
    while k <= K:
        x = rd.randint(2, n - 1)
        if extended_euclid(x, n)[0] > 1:
            return False
        x_d = pow(x, d, n)
        if x_d == 1 or x_d == n - 1:
            k += 1
        else:
            fl = False
            for r in range(1, s):
                x_r = pow(x, d * 2**r, n)
                if x_r == n - 1:
                    fl = True
                    k += 1
                    break
                elif x_r == 1:
                    return False
                else:
                    continue
            if not fl:
                return False
    return True

def main():
    cur192 = {
        'p': 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF,
        'x': 0x44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B,
        'y': 0x15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029,
        'n': 0x3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455,
        'd': 0x6DBA6A,
    }
    

    # Find 521bit curve
    d = 363
    cur521 = {
        'p': 2**521-1,
    }
    print "finding curve... begin with d =", d
    while True:
        if jakobi_symbol(d, cur521['p']) != 1:
            cur521['d'] = d
            a,b = get_isomorphic_curve(cur521)
            E = EllipticCurve(GF(2**521-1), [0,0,0,a,b])

            n = E.cardinality()
            t = (cur521['p']+1) - n
            n1 = n//4
            n2 = ((cur521['p']+1)+t)//4
            if miller_rabine_test(n1):
                cur521['n'] = n1
                print "found curve!"
                if check_mov(cur521):
                    print "mov checked, its our curve!"
                    print cur521['d'], cur521['n']
                    break
            elif miller_rabine_test(n2):
                cur521['n'] = n2
                cur521['d'] = inv_mod(cur521['d'], cur521['p'])
                print "found curve! (twisted)"
                if check_mov(cur521):
                    print "mov checked, its our curve!"
                    print cur521['d'], cur521['n']
                    break
            else:
                print "curve is not appropriate, skip", d
            d += 1
                
        else:
            d += 1

if __name__ == '__main__':
    main()