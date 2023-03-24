# attach('~/overpartitions/eta-quotients.sage') # can be obtained from https://github.com/nsirolli/eta-quotients
# etaq = EtaQuotient({2:1,1:-2}) # needs the above attachment to work

nmin = 2000 # for smaller n the naive method must be used

from os.path import exists
oppath= os.path.expanduser('~')+'/overpartitions/'


def nksupport(n,k,supp=None,zeros=None):
    '''
        Returns a 3-tuple consisting of:
        - a dictionary whose keys are the pairs needed to compute Atilde(n,k),
          according to Theorem 6.1, and whose values at each key are -1 or
          theta, according to Theorem 6.2.
        - a list of pairs (n1,k1) at which Atilde vanishes.
        - the number of cosines to be multiplied.
    '''
    if not supp:
        supp = {}
    if not zeros:
        zeros = set()
    ek = 0
    res = supp, zeros, ek
    if (n,k) in zeros:
        return res
    n3 = n
    k2 = k
    for p,e in factor(k):
        k1 = p**e
        k2 = k2//k1
        n1 = n3 / k2**2 % k1
        n2 = n3 / k1**2 % k2
        n3 = n2
        if (n1,k1) in zeros or (n2,k2) in zeros:
            zeros.add((n,k))
            return res
        else:
            if n1 % p == 0:
                if e > 1:
                    zeros.add((n,k))
                    zeros.add((n1,k1))
                    return res
                else:
                    resk1 = -1 # meaning here there is no factor actually
            else:
                n1m = mod(n1,k1)
                if is_square(-n1m):
                    th = sqrt(-n1m/16)
                    resk1 = th.lift()
                    ek += 1
                else:
                    zeros.add((n,k))
                    zeros.add((n1,k1))
                    return res
            supp[(n1,k1)] = resk1
        if k2 == 1:
            break
    res = supp, zeros, ek
    return res

def nsupport(n,supp=None,zeros=None,N=None):
    '''
        Gathers the output of nksupport(n,k) for k up to N.
    '''
    if not supp:
        supp = {}
    if not zeros:
        zeros = set()
    ekdict = {}
    if not N:
        N = ceil(sqrt(n))
    for k in xsrange(3,N+1,2):
        suppk, zerosk, ek = nksupport(n,k,supp,zeros)
        supp |= suppk
        zeros = zeros.union(zerosk)
        if ek != 0:
            ekdict[(n,k)] = ek
    res = supp, zeros, ekdict
    return res

def Atilde(n,k,supp=None,zeros=None,prec=None,by_def=False):
    '''
        Computes Atilde(n,k)/sqrt(k) from the data obtained from nksupport(n,k).
        '''
    if k == 1:
        return 1
    if by_def:
        sums = [exp(pi*I*(2*dedekind_sum(h,k)-dedekind_sum(2*h,k)-2*n*h/k)) \
                for h in xsrange(k) if gcd(k,h) == 1]
        res = sum(sums) / sqrt(k)
        return res.real()
    if not supp or not zeros:
        supp, zeros, _ = nksupport(n,k)
    if (n,k) in zeros:
        return 0
    if prec:
        RR = RealField(prec)
    n3 = n
    k2 = k
    res = 1
    for p,e in factor(k):  # Algorithm 8.1 in the paper
        k1 = p**e
        k2 = k2//k1
        n1 = n3 / k2**2 % k1
        n2 = n3 / k1**2 % k2
        if (n1,k1) in zeros or (n2,k2) in zeros:
            return 0
        th = supp[(n1,k1)]
        if th != -1:
            x = pi*4*th/k1
            if prec:
                x = RR(x)
            res *= 2*cos(x)
        if k2 == 1:
            break
        n3 = n2
    return res

def nkprec(n,k,ek):
    '''
        Returns the precision that is needed for computing tk, according to
        Proposition 8.2.
        Here ek is the number of cosines that are to be multiplied.
    '''
    cn = pi*sqrt(n)
    x = cn/k
    bits = ceil(max([
        -1/2*log(n,2) + x*log(e,2) + ek + \
                log(10*x + 7*(2*ek-4) + 22,2), \
        1/2*log(n,2) + 5,\
        11
        ]))
    bits += 10 # for safety
    return bits

def pbar(n,suppz=None,ver=0,l=None):
    '''
        different versions:
        0: naive, using q-expansions
        1: computes the HRR sum symbolically
        2: following Algorithm 8.2
        3: modulo l, following Algorithm 8.3
    '''
    if ver == 0 or n < nmin:
        # the remaining versions only are known to work for n >= nmin
        res = qexpansion(n+1,l)[n]
    else:
        if not suppz:
            supp, zeros, ekdict = nsupport(n)
        N = ceil(sqrt(n))
        cn = pi*sqrt(n)
        ent2 = ceil(log(1/8/n*exp(cn),2)) - 1 # for the integer part
        dec2 = ceil(log(4*sqrt(n), 2)) + 1 # for the fractional part
        ext2 = 10 # for safety
        dig2 = ent2 + dec2 + ext2 # precision needed to store the result
        RR = RealField(dig2)

        if ver == 1:
            sums = [Atilde(n,k,supp,zeros) * U(cn/k) \
                    for k in xsrange(1,N+1,2)] # symbolic
            res = 1/4/n * sum(sums)
            res = RR(res)
            res = round(res)

        if ver == 2:
            res = 0
            for k in xsrange(1,N+1,2):
                if (n,k) in zeros:
                    continue
                try:
                    ek = ekdict[(n,k)]
                except:
                    ek = 0
                bits = nkprec(n,k,ek)
                tk = RR(1)
                if ek > 0:
                    tk *= RR(Atilde(n,k,supp,zeros,prec=bits))
                tk *= RR(U(cn/k,prec=bits))
                res += tk
            res *= 1/4/n
            res = round(res)

        if ver == 3:
            r = 0
            s = 0 # will belong to [0,1)
            parity = 0
            for k in xsrange(1,N+1,2):
                if (n,k) in zeros:
                    continue
                try:
                    ek = ekdict[(n,k)]
                except:
                    ek = 0
                bits = nkprec(n,k,ek)
                bits += 10
                tk = RR(1)
                if ek > 0:
                    tk *= RR(Atilde(n,k,supp,zeros,prec=bits))
                tk *= RR(U(cn/k,prec=bits)/4/n)
                rk = floor(tk) # we could also round!
                sk = tk - rk
                # for safety
                assert sk >= 0
                assert sk < 1, f'for k = {k} we have that sk = {sk}'
                s += sk
                # for safety
                assert s < 2
                if s >= 1:
                    rk += 1
                    s -= 1
                # for safety
                assert s < 1
                assert s >= 0
                parity += rk % 2
                parity = parity % 2
                r += rk % l
                r = r % l
            res = (r+parity) % l

    return res

def U(x,prec=None):
    if prec:
        RR = RealField(prec)
        x = RR(x)
    return cosh(x) - sinh(x)/x

def nsfile(l,Q,path=oppath+'ns/'):
    '''
        Creates a file with the needed n's where pbar is to be evaluated when
        running Algortihm 9.1.
        Each line of the file contains a 3-tuple (n,n*Q**2,n/Q**2).
    '''
    p = prime_divisors(l)[0]
    j = valuation(l,p)
    try:
        eps,nmax,_ = ldata[l]
    except:
        cc = etaq.conditionC(p)
        assert cc[0]
        eps = cc[1]
        nmax = etaq.sturm_bound(p,j)
    ns = []
    for n in xsrange(1,nmax+1):
        if kronecker(n,p) == -eps:
            if n % (Q**2) == 0:
                ndiv = n // (Q**2)
            else:
                ndiv = 0
            ns.append((n,n*Q**2,ndiv))
    file = path + f'ns-l{l}Q{Q}'
    with open(file,'w') as f:
        for n in ns:
            f.write(str(n) + '\n')

# given a pair (l,Q), the idea for deciding efficiently if Algorithm 9.1 returns
# True is to split the file given by nsfile(l,Q) into chunks, and run in
# parallel the following command over the different chunks.

def congruences(l,Q,ns,path=''):
    '''
        Given a file ns containing 3-tuples (n,n*Q**2,n/Q**2), verifies the
        condition in line 12 of Algorithm 9.1 for those n, by computing the
        values of pbar at the components of these tuples; unless there
        exists a file unint claiming that Q does not yield congruences.
        If for some line in ns the condition does not hold, it writes the
        relevant information in a file unint.
        Otherwise, it writes a file with the computed values.
    '''
    unint = oppath+f'congruences/u-l{l}Q{Q}'
    pbarns = path+f'pbar-l{l}Q{Q}ns{os.path.basename(ns)}' 
    eps, _, kappa = ldata[l]
    sign = (-1)**((kappa-1)//2)
    Q2 = Q**(kappa-2) % l
    Q3 = Q**((kappa-3)//2) % l

    with open(ns,'r') as ns:
        with open(pbarns,'w') as ff:
            for npd in ns:
                if not exists(unint):
                    n, npor, ndiv = sage_eval(npd)
                    pbn = pbar(n,ver=3,l=l)
                    ff.write(str((n,pbn)) + '\n')
                    pbnpor = pbar(npor,ver=3,l=l)
                    ff.write(str((npor,pbnpor)) + '\n')
                    if ndiv != 0:
                        pbndiv = pbar(ndiv,ver=3,l=l)
                        ff.write(str((ndiv,pbndiv)) + '\n')
                    else:
                        pbndiv = 0
                    cond = pbnpor \
                            + kronecker_symbol(sign *n,Q) * Q3 * pbn \
                            + Q2 * pbndiv
                    cond =  cond % l
                    if cond != 0:
                        with open(unint,'a+') as f:
                            f.write(f'{n}\n')
                            f.write(f'{(n,pbn)}\n')
                            f.write(f'{(npor,pbnpor)}\n')
                            f.write(f'{(ndiv,pbndiv)}\n')
                            f.write(f'{cond}\n')
                            f.write('\n')

# code for computing q-expansions, borrowed from
# https://github.com/nsirolli/eta-quotients; included here to avoid the need to
# import it
from sage.modular.etaproducts import qexp_eta

def qexp_eta_scaled(ps_ring, prec, d):
    """
        computes qexp_eta(R,prec)(x^d)
    """
    prec = Integer(prec)
    if not prec > 0:
        raise ValueError("prec must be a positive integer")
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    try:
        n = 1
        while True:
            pm = -pm
            v[d*n*(3*n-1)//2] = pm
            v[d*n*(3*n+1)//2] = pm
            n += 1
    except IndexError:
        pass
    return ps_ring(v, prec=prec)

def qexpansion(prec,l=None):
    if not l:
        coeff_ring = ZZ
    else:
        coeff_ring = IntegerModRing(l)
    assert prec >= 0
    v_int = 0 
    rdict = {2:1,1:-2}
    ring = coeff_ring[['q']]
    qw = ring.0
    Nprec = prec - v_int
    eta_q = 1
    if prec <= v_int:
        return eta_q * O(qw**prec)
    for delta in rdict:
        rd = rdict[delta]
        ff = qexp_eta_scaled(ring, Nprec, delta) 
        eta_q = eta_q*ff**rd
    return qw**v_int * eta_q

# here we store the needed computations from
# https://github.com/nsirolli/eta-quotients: for each prime l, the 3-tuple
# consisting of eps_l, nmax and kappa.
ldata = {\
        3:(-1,852,71),\
        9:(-1,852,71),\
        27:(-1,2580,215),\
        81:(-1, 7764, 647),\
        5:(1,3570,119),\
        25:(1,3570,119),\
        125:(1,17970,599),\
        7:(-1,2632,47),\
        49:(-1, 18760, 335),\
        11:(-1,15708,119),\
        13:(1,30394,167),\
        17:(1,87822, 287),\
        19:(-1,136420, 359),\
        23:(-1, 290904, 527)\
        }
