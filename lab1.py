import math
import random
import numpy as np
from itertools import chain, combinations
import time

def jacobi(a: int, n: int) -> int:
    if math.gcd(a, n) != 1:
        return 0
    res = 1
    while a != 1:
        if a < 0:
            a = -a
            res *= (-1) ** ((n - 1) // 2)
        while a % 2 == 0:
            a //= 2
            res *= (-1) ** ((n ** 2 - 1) // 8)
        if a == 1:
            break
        if a < n:
            temp = a
            a = n
            n = temp
            res *= (-1) ** ((n - 1) // 2 * ((a - 1) // 2))
        if a > n:
            a %= n
    if a == 1:
        return 1 * res

def pow_mod(t: int, k: int, mod: int) -> int:
    res = 1
    while k:
        if k & 1:
            res *= t
            res %= mod
        k //= 2
        if k == 0:
            break
        t = t ** 2
        t %= mod
    res %= mod
    return res

def pseudo_Euler(a: int, p: int) -> bool:
    if math.gcd(a, p) != 1:
        return False
    else:
        jac = jacobi(a, p)
        if jac == -1:
            jac = p - 1
        deg = int((p - 1) / 2)
        pow = pow_mod(a, deg, p)
        if jac == pow % p:
            return True
        else:
            return False

def solovei_shtrassen(p: int, k: int) -> bool:
    i = 0
    if(p==1 or p==2):
        return True
    while i < k:
        x = random.randint(2, p - 1)
        if math.gcd(x, p) == 1:
            if not pseudo_Euler(a=x, p=p):
                return False
            else:
                i += 1
        else:
            return False
    return True

primes = [2, 3, 5, 7, 11]

def trial_divisions(num):
    if num < 1:
        return -1

    for d in primes:
        if(num%d==0):
            return d
    return 1

def f(x, num):
    return (x ** 2 + 1) % num

def rho_pollard(num: int, x_0: int = 2) -> int:
    d=1
    x = f(x_0, num)
    y = f(f(x_0, num), num)

    while x != y:
        x = f(x, num)
        y = f(f(y, num), num)

        if x == y:
            return 1

        d = math.gcd(num, (x-y) % num)
        if d != 1:
            return d

    return 1

def calculate_a(n, sqrt_n, V, Al, A, U):
    V.append((n - U[-1] ** 2) / V[-1])
    Al.append((sqrt_n + U[-1]) / V[-1])
    A.append(math.floor(Al[-1]))
    U.append(A[-1] * V[-1] - U[-1])

def calculate_b(n, B, B2, A):
    B.append((B[-1] * A[-1] + B[-2]) % n)
    b2 = pow(B[-1], 2, n)
    B2.append(b2 if b2 < n / 2 else b2 - n)

def legendre_symbol(num, prime):
    if num % prime == 0:
        return 0
    return pow(num, int((prime - 1) / 2), prime)

#перевірка що символ Лежандра =1
def check_divisor(n, La, p):
    ls = legendre_symbol(n, p)
    return ls == 1 and p < La

def factorBase(num, n, La):
    P = {}
    if num < 0:
        P[-1] = 1
        num *= -1
    p = 2
    while pow(p, 2) <= num:
        if num % p == 0:
            P[p] = 1
            if not check_divisor(n, La, p):
                return False
            num /= p
            while num % p == 0:
                P[p] = 1 - P[p]
                num /= p
        p += 1
    if num > 1:
        num = int(num)
        P[num] = 1
        if not check_divisor(n, La, p):
            return False
    Fb = {p: 1 for p in P}
    return True

def build_vectors(Fb, Vec_b, vec_len, P, B, B2):
    append_size = len(Fb) - vec_len
    for vec in Vec_b.values():
        for i in range(append_size):
            vec.append(0)
    vec = []
    for i in Fb:
        if i in P:
            vec.append(P[i])
        else:
            vec.append(0)
    Vec_b[(B[-1], B2[-1])] = vec

def xor(vector_a, vector_b):
    vector_c = vector_a[:]
    for i in range(len(vector_a)):
        vector_c[i] = (vector_c[i] + vector_b[i]) % 2
    return vector_c

def key(dictionary, value):
    for key, val in dictionary.items():
        if val == value:
            return key
    return None

def solve(Vec_b_dict, Fb_dict, n):
    Vec_power = list(chain.from_iterable(combinations(Vec_b_dict.values(), r) for r in range(1, len(Vec_b_dict.values()) + 1)))
    for Vecs in Vec_power:
        sum_vector = [0] * len(Fb_dict)
        for vector in Vecs:
            sum_vector = xor(sum_vector, vector)
        if sum_vector == [0] * len(Fb_dict):
            X = 1
            Y = 1
            for vector in Vecs:
                b, b2 = key(Vec_b_dict, vector)
                X = (X * b) % n
                Y = (Y * b2) % (n ** 2)
            Y = int(Y ** 0.5)
            if X != Y and X != n - Y:
                d1 = math.gcd(X + Y, n)
                d2 = math.gcd(X - Y, n)
                if d1 != 1 and d2 != 1:
                    return (d1, d2)
                global alpha
                alpha += 1
    return (-1, -1)



def Brilhart_Morrison(n):
    sqrt_n = n ** 0.5
    global alpha 
    alpha = (1 / 2) ** 0.5
    La = math.e ** (alpha * (math.log2(n) * math.log2(math.log2(n))) ** 0.5)

    V = [1]
    Al = [sqrt_n]
    A = [math.floor(Al[0])]
    U = [A[0]]
    B = [0, 1]
    B2 = [0, 1]
    P = {}
    Fb = {}
    Vec_b = {}
    for i in range(1, 100):
        vec_len = len(Fb)
        if i != 1:
            calculate_a(n, sqrt_n, V, Al, A, U)
        calculate_b(n, B, B2, A)
        i += 1
        if not factorBase(B2[-1], n, La):
            continue
        build_vectors(Fb, Vec_b, vec_len, P, B, B2)
        res = solve(Vec_b, Fb, n)
        if res != (-1, -1):
            return res[0]
        print()
    return n

def check_prime(a):
    if a < 2:
        return False
    for i in range(2, int(math.sqrt(a)) + 1):
        if a % i == 0:
            return False
    return True

def canonical_distribution(n):
    result=[]
    k = random.randrange(3, max(n + 1,4))
    if solovei_shtrassen(n,k):
        return [n]
    while True:
        if(n<=47):
            d=trial_divisions(n)
            if d==1: 
                break
            #print("trial, division", d)
            result.append(d)
            n=n//d
        else:
            break
    d=rho_pollard(n)
    if(d!=1):
        #print("rho_pollard",d)
        temp=d
        if not check_prime(d):
            temp=canonical_distribution(d)
        result.append(temp)
        n//=d
    while True:
        k = random.randrange(3, max(n + 1,4))
        if solovei_shtrassen(n,k):
            result.append(n)
            return result
        d1=Brilhart_Morrison(n)
        if d1==1:
            result.append(n)
            print("I can't find canonical schedule of numbers :(")
            return result
        n//=d1
        #print("Brilhart_Morrison", d1)
        result.append(d1)
        


print(canonical_distribution(901667173167834173))
