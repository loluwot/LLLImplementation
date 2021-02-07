from vector import *
import itertools
from decimal import *
from Crypto.Util.number import *

getcontext().prec = 2000
def gram_schmidt(basis):
    new_basis = [basis[0]]
    for i in range(1, len(basis)):
        adjust = Vector([0 for _ in range(len(basis[0].coords))])
        for j in range(len(new_basis)):
            adjust += new_basis[j]*Decimal((basis[i].dot(new_basis[j])))/Decimal(new_basis[j].sqmag())

        new_basis.append(basis[i] - adjust)
        
    return new_basis

ans = gram_schmidt([Vector([4,1,3,-1]), Vector([2,1,-3,4]), Vector([1,0,-2,7]), Vector([6,2,9,-5])])

def mucalc(grambasis, basis):
    #print('MU')
    mu = [[0 for _ in range(len(basis))] for __ in range(len(basis))]
    for i in range(len(grambasis)):
        for j in range(i):
            # if not i == j:
                #print(i, j)
                #print(str(grambasis[j]))
            mu[i][j] = Decimal(basis[i].dot(grambasis[j]))/Decimal(grambasis[j].sqmag())
    return mu

def lll(basis, delta=Decimal(0.5)):
    print(len(basis))
    grambasis = gram_schmidt(basis)
    newbasis = basis
    mu = mucalc(grambasis, basis)
    k = 1
    while k < len(grambasis):
        for j in range(k-1, -1, -1):
            #print(float(mu[k][j]))
            if abs(mu[k][j]) > 0.5:
                #print('CONDITION 1')
                newbasis[k] -= newbasis[j]*int(round(mu[k][j]))
                #print(newbasis[k])
                grambasis = gram_schmidt(newbasis)
                mu = mucalc(grambasis, newbasis)
        if grambasis[k].sqmag() >= (delta - mu[k][k-1]**2)*grambasis[k-1].sqmag():
            k += 1
        else:
            #print('CONDITION2 ')
            temp = newbasis[k]
            newbasis[k] = newbasis[k-1]
            newbasis[k-1] = temp
            grambasis = gram_schmidt(newbasis)
            mu = mucalc(grambasis, newbasis)
            #newbasis[k-1] = temp
            k = max(k-1, 1)
    
    return newbasis

def constants_approx(num, constants_list=[(Decimal(2).sqrt(), 'sqrt(2)'), (Decimal(3).sqrt(), 'sqrt(3)')], max_degree=3):
    ACCURACY = 10**9
    constants, names = list(zip(*constants_list))
    constants = list(constants)
    constants.append(1)
    terms = set()
    for tup in itertools.product([i for i in range(len(constants))], repeat=max_degree):
        prod = 1
        count = [0 for _ in range(len(constants)-1)]
        for idx in tup:
            if idx == len(constants) - 1:
                continue
            prod *= constants[idx]
            count[idx] += 1
        count = tuple(count)
        if prod == 1:
            continue
        prod = math.floor(prod*ACCURACY)   
        terms.add((prod, count))
    basis = []
    terms = list(terms)
    print(terms)
    for i, tup in enumerate(terms):
        term, name = tup
        coords = [0 for _ in range(len(terms) + 2)]
        coords[i] = 1
        coords[-1] = term
        print(term)
        basis.append(Vector(coords))
    
    coords = [0 for _ in range(len(terms) + 2)]
    coords[-2] = 1
    coords[-1] = -math.floor(num*ACCURACY)
    basis.append(Vector(coords))
    for b in basis:
        print(str(b))
    rbasis = lll(basis)
    for b in rbasis:
        print(str(b))
    smallest = rbasis[-1]
    #print(smallest.coords[-1])
    s = []
    su = 0
    for i, v in enumerate(smallest.coords[:-2]):
        if not v == 0:
            temp = terms[i][1]
            st = []
            value = 1
            for ii, vv in enumerate(temp):
                if not vv == 0:
                    if vv == 1:
                        st.append('{}'.format(names[ii]))
                    else:
                        st.append('{}**{}'.format(names[ii],vv))
                    value *= constants[ii] ** vv
            st = '*'.join(st)
            if len(st) == 0:
                s.append('{}'.format(v))
            else:
                s.append('{}*{}'.format(v, st))
            su += Decimal(float(v))*value
    return '+'.join(s)


def algebraic_approx(num, max_degree=3, accuracy=10):
    ACCURACY = 10**accuracy
    basis = []
    for i in range(max_degree):
        coords = [0 for _ in range(max_degree + 1)]
        coords[i] = 1
        coords[-1] = int(ACCURACY*Decimal(num)**Decimal(i))
        basis.append(Vector(coords))
    rbasis = lll(basis)
    smallest = rbasis[0].coords[:-1]
    # terms = []
    # for i, v in enumerate(smallest[:-1]):
    #     if not v == 0:
    #         if v == 1:
    #             terms.append('x**{}'.format(i))
    #         else:
    #             terms.append('{}*x**{}'.format(v, i))
    # return terms
    return smallest
    
def quadratic_approx(num, accuracy=10):
    c, b, a = algebraic_approx(num, max_degree=3, accuracy=accuracy)
    print(c, b, a)
    s = (-b + Decimal(b**2 - 4*a*c))/2/a
    d = (-b - Decimal(b**2 - 4*a*c))/2/a
    pm = '+' if abs(s - num) < abs(d - num) else '-'
    
    return '({}{}sqrt({}))/{}'.format(-b, pm, b**2 - 4*a*c, 2*a)
    
#constants_list = [(Decimal(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706), 'pi'), (Decimal(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274), 'e')]
#print(quadratic_approx(Decimal(9.81)))
#print(algebraic_approx(Decimal(3.14159)))
#print(constants_approx(Decimal(63710088)/10**4, constants_list=constants_list))
#print(constants_approx(Decimal()))
def isint(val):
    return val-int(val) == 0

#Coppersmith's attack

p = getPrime(512)
dp = p - (p % 2**86)
#print(p%2**86)
q = getPrime(512)
N = p*q

X = 2**86
#print(X)

basis = [Vector([X**2, 2*X*dp, dp**2]), Vector([0, X, dp]), Vector([0, 0, N])]
rbasis = lll(basis)
a, b, c = (rbasis[0]).coords
a = a // X**2
b = b // X
#print(a, b, c)
res1 = (-b - int(Decimal(b**2 - 4*a*c).sqrt()))/(2*a)
res2 = (-b + int(Decimal(b**2 - 4*a*c).sqrt()))/(2*a)
if isint(res1):
    print((dp + res1) == p)
else:
    print((dp + res2) == p)
#print((a + res) % p)
