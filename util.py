import random as rd
from sage.all import *

def check_point_on_curve(cur, x, y):
    return (x**2 + y**2) % cur['p'] == (1 + cur['d'] * x**2 * y**2) % cur['p']

def inv_mod(a ,n):
    return pow(a, n-2, n)

def sqrt_mod(a, n):
    k = (n-3)//4 # for n=4k+3 only!
    return pow(a, k+1, n)

def jakobi_symbol(n, p):
    return pow(n, (p-1)/2, p)

def get_isomorphic_curve(cur):
    d = cur['d']
    p = cur['p']
    return ( ((-(1 + 14*d + d**2 ) % p) * inv_mod(48, p)) % p, ((-(1 - 33*d - 33 * d**2 + d**3) % p) * inv_mod(864, p)) % p )

def get_base_point(cur):
    # x**2 + y**2 = 1 + d * x**2 * y**2
    while True:
        x = rd.randint(1, cur['p'])
        y_2 = (1-x**2 % cur['p']) * inv_mod(1- cur['d']*x**2, cur['p']) % cur['p']
        if jakobi_symbol(y_2, cur['p']) == 1:
            return x, sqrt_mod(y_2, cur['p'])


def check_mov(cur):
    n = cur['n']
    p = cur['p']

    for k in range(2, 33):
        if (p**k - 1) % n == 0:
            return False
    return True



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

def find521curve(d):
    # Find 521bit curve
    
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
                    print get_base_point(cur521)
                    break
            elif miller_rabine_test(n2):
                cur521['n'] = n2
                cur521['d'] = inv_mod(cur521['d'], cur521['p'])
                print "found curve! (twisted)"
                if check_mov(cur521):
                    print "mov checked, its our curve!"
                    print cur521['d'], cur521['n']
                    print get_base_point(cur521)
                    break
            else:
                print "curve is not appropriate, skip", d
            d += 1
                
        else:
            d += 1

def main():
    cur521 = {
        'd': 0x16a,
        'n': 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff46087bc4a294fcc80b3f45d8cebfb21479ba651ba07de913ad1d8392de3ff8af,
        'p': 2**521 - 1
    }
    x, y = get_base_point(cur521)
    if check_point_on_curve(cur521, x, y):
        print hex(x), hex(y)
    else:
        print "smth gone wrong"
    
    find521curve(750)

    

if __name__ == '__main__':
    main()