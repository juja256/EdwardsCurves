import random as rd
from sage.all import *

def check_point_on_curve(cur, x, y):
    if not 'a' in cur:
        cur['a'] = 1
    return (x**2 + cur['a']*y**2) % cur['p'] == (1 + cur['d'] * x**2 * y**2) % cur['p']

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
    a = 1
    if a in cur:
        a = cur['a']
    
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

def edwards_j(d,a,p):
    return ((16* pow((pow(a, 2, p) + pow(d, 2, p) + 14 * a * d) % p, 3, p)) * inv_mod(a * d * pow(a-d % p, 4, p), p)) % p

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


    curs = {
        'UA_384_2': {
            'p': 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF115,
            'n': 0x400000000000000000000000000000000000000000000000231C4F8D91A6B595EAA0789F9CCFF3C7FA9F05F6C028EACF,
            'd': 0xA8,
            'a': 2,
            'y': 0xA819A5ADFFD39342A4AF74E1A1B5EEDEDA54DD1BD1E5EA0E871F325F4702736E3336D2048AA8BB89DDA7DA6B3C969AEF,
            'x': 0x9B8C366A8BE70A632BBEB36B1C164998B949953EDBEC19CFD7065C8DD197823285C8254CAEE2015CB36597827F6A37E2
            #'y': 0x89FD5124BBCFC2FBFA908A0A2F8D46E9A443EA0D34A8101CC28EA068C32EEA2A49466ECEDD25E2DABECBDE0B016C8ACD
            #'x': 0x45DF45F8020CA03A0DE417410BDA8BFB795A653450104321344EFB9A1C5086B25A970EC79ED51DCA9D362AFF9A86F528
        },
        'UA_384_1': {
            'p': 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF79D,
            'n': 0x4000000000000000000000000000000000000000000000000AF905B73674AC7D4AF38C53331DC208A517DCB3F340EECF,
            'd': 0x214,
            'a': 2,
            'y': 0xCA2C0B2E41084153BDBB86F443452DAB96813E26194F8428EFAFE1F84B0B7D1810EBE419B398C6E6B9160187BA7017E0,
            'x': 0x3EC279CBB1230F32704B13BC0BF90C86518AD689A3BCC59701B5D041D7A97A9AE838B359AD9BD3E9AD756A99724582A3
            #'y': 0x70A53749F80DD864332CB6E9672B15E1BBCA81BF473F3EF45A4963EF3782AA779DC83176F0E9088CA1C9E1375335E73B,
            #'x': 0xA7F96A06EC6BB1EB6210E910D185AE7E95642587736E7934ED2FBEB5142F11038E8A7E222474AA4CDD82205009DEC375
        },
        'UA_256_5': {
            'p': 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF9FD,
            'n': 0x40000000000000000000000000000000412E994309D80E83AFA4642085233F05,
            'd': 0x170,
            'a': 2,
            'y': 0xCC4A085C63F3F5A59306AF8B234F1DFA9A177A62922967D7250B99F797F4FAF9,
            'x': 0xD9D411C794BC01F7A07DEC0ECD314F8DE16720623DF9BA1AE4E08339E067E746
        }
    }

    for name, cur in curs.items():
        print("Verifying " + name)
        print("Base Point Valid: " + str(check_point_on_curve(cur, cur['x'], cur['y'])))
        jinv = edwards_j(cur['d'], cur['a'], cur['p'])
        print("J invariant: " + hex(jinv))
        E = EllipticCurve(GF(cur['p']), j=jinv)
        n = E.cardinality()
        if n//4 == cur['n']:
            print(name + " curve order is correct")
        else:
            print(name + " curve order is invalid, valid = " + hex(n//4).upper())
    

    
    #find521curve(750)
    

if __name__ == '__main__':
    main()