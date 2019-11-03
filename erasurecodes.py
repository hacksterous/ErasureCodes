#GF routines from:
#https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
#Does not use polynomial division

GF_EXP = [0] * 512 # Create list of 512 elements. In Python 2.6+, consider using bytearray
GF_LOG = [0] * 256

def init_tables(prim=0x11d):
    '''Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial.'''
    # prim is the primitive (binary) polynomial. Since it's a polynomial in the binary sense,
    # it's only in fact a single galois field value between 0 and 255, and not a list of gf values.
    global GF_EXP, GF_LOG
    GF_EXP = [0] * 512 # anti-log (exponential) table
    GF_LOG = [0] * 256 # log table
    # For each possible value in the galois field 2^8, we will pre-compute the logarithm and anti-logarithm (exponential) of this value
    x = 1
    for i in range(0, 255):
        GF_EXP[i] = x # compute anti-log for this value and store it in a table
        GF_LOG[x] = i # compute log at the same time
        x = gf_mult_noLUT(x, 2, prim)

        # If you use only generator==2 or a power of 2, you can use the following which is faster than gf_mult_noLUT():
        #x <<= 1 # multiply by 2 (change 1 by another number y to multiply by a power of 2^y)
        #if x & 0x100: # similar to x >= 256, but a lot faster (because 0x100 == 256)
            #x ^= prim # substract the primary polynomial to the current value (instead of 255, so that we get a unique set made of coprime numbers), this is the core of the tables generation

    # Optimization: double the size of the anti-log table so that we don't need to mod 255 to
    # stay inside the bounds (because we will mainly use this table for the multiplication of two GF numbers, no more).
    for i in range(255, 512):
        GF_EXP[i] = GF_EXP[i - 255]

def gf_mult_noLUT(x, y, prim=0, field_charac_full=256, carryless=True):
    '''Galois Field integer multiplication using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction).'''
    r = 0
    while y: # while y is above 0
        if y & 1: r = r ^ x if carryless else r + x # y is odd, then add the corresponding x to r (the sum of all x's corresponding to odd y's will give the final product). Note that since we're in GF(2), the addition is in fact an XOR (very important because in GF(2) the multiplication and additions are carry-less, thus it changes the result!).
        y = y >> 1 # equivalent to y // 2
        x = x << 1 # equivalent to x*2
        if prim > 0 and x & field_charac_full: x = x ^ prim # GF modulo: if x >= 256 then apply modular reduction using the primitive polynomial (we just subtract, but since the primitive number can be above 256 then we directly XOR).

    return r

def gf_mul(x,y):
    if x==0 or y==0:
        return 0
    i = GF_LOG[x] + GF_LOG[y]
    v = GF_EXP[i]
    #print ("gf_mul: gf_exp value is ", v)
    return GF_EXP[i] 
    # should be GF_EXP[(GF_LOG[x]+GF_LOG[y])%255] if GF_EXP wasn't oversized

def gf_exp(x, power):
    return GF_EXP[(GF_LOG[x] * power) % 255]

def gf_inverse(x):
    return GF_EXP[255 - GF_LOG[x]] # gf_inverse(x) == gf_div(1, x)

def gf_div(x,y):
    if y==0:
        return 0
    if x==0:
        return 0
    return GF_EXP[(GF_LOG[x] + 255 - GF_LOG[y]) % 255]

def gf_add(x, y):
    return x ^ y

def gf_sub(x, y):
    return x ^ y # in binary galois field, subtraction is just the same as addition (since we mod 2)

def mat_submatrix (m, sr, sc, lr, lc):
    #extract submatrix from matrix m
    #pass start/end rows and columns
    #rows/columns are counted starting from 1
    nr = lr-sr+1
    nc = lc-sc+1
    
    if nr > len(m) or len(m) == 0:
        return []
    if nc > len(m[0]): 
        return []
    mat = [[0 for j in range(nc)] for i in range(nr)]
    try:
        for i in range(nr):
            for j in range(nc):
                mat[i][j] = m[sr + i -1][sc + j -1]
    except IndexError:
        return []
    return mat

def mat_submatrix (m, sr, sc, lr, lc):
    #extract submatrix from matrix m
    #pass start/end rows and columns
    #rows/columns are counted starting from 1
    nr = lr-sr+1
    nc = lc-sc+1
    
    if nr > len(m) or len(m) == 0:
        return []
    if nc > len(m[0]): 
        return []
    mat = [[0 for j in range(nc)] for i in range(nr)]
    try:
        for i in range(nr):
            for j in range(nc):
                mat[i][j] = m[sr + i -1][sc + j -1]
    except IndexError:
        return []
    return mat

def mat_submat (m, sr, sc, lr, lc):
    #extract submatrix from matrix m
    #pass start/end rows and columns
    #rows/columns are counted starting from 1
    try:
        return [[m[i][j] for j in range(sc-1, lc)] for i in range(sr-1, lr)]
    except IndexError:
        return []

def mat_unity (n):
    mat = [[(1 if i == j else 0) for j in range(n)] for i in range(n)]
    return mat

def mat_vandermonde (nr, nc):
    # Create a Vandermonde matrix, which is guaranteed to have the
    # property that any subset of rows that forms a square matrix
    # is invertible.
    mat = [[0 for j in range(nc)] for i in range(nr)]
    for i in range (nr):
        for j in range(nc):
            mat[i][j] = gf_exp(i, j)
    return mat

def mat_print (mat, h=True):
    nr = len(mat)
    nc = 0
    if nr > 0:
        nc = len(mat[0])
    if nc == 0 or nr == 0:
        print ('[]')
        return

    for i in range (nr):
        for j in range(nc):
            #print (hex(mat[i*nc + j]), ",", end='')
            if h == True:
                print (hex(int(mat[i][j])), '\t', end='')
            else:
                print (mat[i][j], '\t', end='')
        print ('')


def mat_mul (m1, m2, mulfunc=None, addfunc=None):
    #m1 is nr1 x nc1
    #m2 is nc2 x nr2
    #nc1 must be equal to nr2
    #result will be nr1 x nc2
    try:
        nr1 = len(m1)
        #print ("mat_mul: nr1 is ", nr1)
        nc1 = 0
        if nr1 > 0:
            nc1 = len(m1[0])
        #print ("mat_mul: nc1 is ", nc1)
        if nc1 == 0 or nr1 == 0:
            return []
        
        nr2 = len(m2)
        nc2 = 0
        if nr2 > 0:
            nc2 = len(m2[0])
        #print ("mat_mul: nr2 is ", nr2)
        #print ("mat_mul: nc2 is ", nc2)
        if nc2 == 0 or nr2 == 0:
            return []
        
        if nc1 != nr2:
            return []
    except (IndexError, ValueError, TypeError):
        return []

    mat = [[0 for j in range(nc2)] for i in range(nr1)]
    for i in range (nr1):
        for j in range(nc2):
            mat[i][j] = 0
            for u in range(nr2):
                if mulfunc == None and addfunc == None:
                    mat[i][j] += m1[i][u] * m2[u][j] 
                else:
                    #print ("mat_mul: passing to gf_mul:", m1[i][u], m2[u][j])
                    mat[i][j] = addfunc(mat[i][j], mulfunc(m1[i][u], m2[u][j]))
            #print ("mat_mul: row", i, "col", j, "is", mat[i][j])
    return mat

def mat_transpose (m):
    try:
        nr = len(m)
        nc = 0
        if nr > 0:
            nc = len(m[0])
        if nc == 0 or nr == 0:
            return []
    except (IndexError, ValueError, TypeError):
        print ("transpose error")
        return []

    return [[m[i][j] for i in range(nr)] for j in range(nc)]
    

def mat_round (m, rnd=1000):
    nr = len(m)
    nc = 0
    if nr > 0:
        nc = len(m[0])
    if nc == 0 or nr == 0:
        return []

    return [[nround(m[i][j], rnd) for j in range(nc)] for i in range(nr)]
    
def mat_invert (m, rnd=3):
    #NOTE: will change the value of m to unity matrix
    #if matrix inversion is successful. Else value
    #might remain unchanged or partially changed.
    n = len(m)
    if n != len(m[0]):
        return []

    mat = mat_unity(n)
    for i in range(n):
        #iterate through each row
        #choose each diagonal element
        if m[i][i] == 0:
            return []
        
        #scale all elements on this row i by the diagonal element
        V = m[i][i]
        for j in range(n):
            m[i][j] /= V
            mat[i][j] /= V

        for k in range(n):
            #iterate through all the other rows
            #subtract from all other rows of this column V*m[i][i]
            V = m[i][i]
            U = m[k][i]/V #use i for column number
            for j in range(n):
                #zero out elements of column j that are not a diagonal
                #element 
                if k != i:
                    m[k][j] -= U*m[i][j]
                    mat[k][j] -= U*mat[i][j]

    return mat

def mat_gf_invert (m, rnd=3):
    #NOTE: will change the value of m to unity matrix
    #if matrix inversion is successful. Else value
    #might remain unchanged or partially changed.
    n = len(m)
    if n != len(m[0]):
        return []

    mat = mat_unity(n)
    for i in range(n):
        #iterate through each row
        #choose each diagonal element
        if m[i][i] == 0:
            return []
        
        #scale all elements on this row i by the diagonal element
        V = m[i][i]
        for j in range(n):
            m[i][j] = gf_div(m[i][j], V)
            mat[i][j] = gf_div(mat[i][j], V)

        for k in range(n):
            #iterate through all the other rows
            #subtract from all other rows of this column V*m[i][i]
            V = m[i][i]
            U = gf_div(m[k][i], V) #use i for column number
            for j in range(n):
                #zero out elements of column j that are not a diagonal
                #element 
                if k != i:
                    m[k][j] = gf_sub(m[k][j], gf_mul(U, m[i][j]))
                    mat[k][j] = gf_sub(mat[k][j], gf_mul(U, mat[i][j]))

    return mat

def mat_cauchy (kset, mset):
    nr = len(mset)
    nc = len(kset)
    mat = [[0 for j in range(nc)] for i in range(nr)]
    for i in range(nr):
        for j in range(nc):
            mat[i][j] = gf_inverse(gf_add(mset[i], kset[j]))
    return mat

def mat_catrows (m1, m2):
    #cat m2 matrix at the bottom of m1
    nr1 = len(m1)
    nr2 = len(m2)
    if nr1 > 0 and nr2 > 0:
        nc = len(m2[0])
    else:
        return []
    if len (m1[0]) != nc:
        return []
    mat = [[(m1[i][j] if i < nr1 else m2[i-nr1][j]) for j in range(nc)] for i in range(nr1+nr2)]
    return mat

def mat_catcols (m1, m2):
    #cat m2 matrix at the bottom of m1
    nr = len(m1)
    if nr == len(m2):
        nc1 = len(m1[0])
        nc2 = len(m2[0])
    else:
        return []
    mat = [[(m1[i][j] if j < nc1 else m2[i][j-nc1]) for j in range(nc1+nc2)] for i in range(nr)]
    return mat

def mat_delrows (m, sr, n):
    nr = len(m)
    if nr <= 0:
        return []
    nc = len(m[0])
    if nc <= 0:
        return []

    try:
        mat = [[m[i][j] for j in range(nc)] for i in range(sr-1)]
        mat += [[m[i][j] for j in range(nc)] for i in range(sr+n-1, nr)]
    except IndexError:
        return []
    return mat

def nround(x, fac):
    return round(x*fac)/fac

# Configuration of the parameters and input message
prim = 0x11d
n = 12      
k = 8       #no. of data disks
m = n - k   #no. of parity disks

# Initializing the log/antilog tables
init_tables(prim)
#print ("finally, GF_EXP is", GF_EXP)

mat = mat_vandermonde (n, k)
#mat_print (mat)
#print ("Extracted generator matrix")
mat = mat_submatrix(mat, k+1, 1, n, k)
#mat_print (mat)

m1 = [[1,2,3], [4,5,6]]
m2 = [[1,2], [3, 4], [5,6]]

#mat = mat_mul (m1, m2, mulfunc=gf_mul, addfunc=gf_add)
#mat_print (mat, h=False)
#mat_print (mat)

#Cauchy matrix
Kset = [251, 1, 2, 3, 45, 155, 6, 7]
Mset = [8, 9, 110, 11]
mat = mat_cauchy (Kset, Mset)
#mat_print (mat_unity(8))
#print ("Cauchy matrix")
#mat_print (mat)
#print ("mat inversion")
A = [[1,2,3], [4, 5, 6], [7, 8, 9]]
A = [[5,3,1], [3, 9, 4], [1, 3, 5]]
mat = mat_invert (A)
#mat_print (mat_round(mat, 1000), h=False)
#print ("mat orig")
A = [[1,2,3], [4, 5, 6], [7, 8, 9]]
A = [[5,3,1], [3, 9, 4], [1, 3, 5]]
#mat_print (A, h=False)
#print ("mat check")
#mat = mat_mul (A, mat, mulfunc=None)
#mat_print (mat_round(mat, 1000), h=False)
A = [[5,3,1], [3, 9, 4], [1, 3, 5]]
mat = mat_invert (A)
A = [[5,3,1], [3, 9, 4], [1, 3, 5], [6, 7, 8]]
#print ("mat catrows")
mat2 = mat_catrows(A, mat)
#mat_print (mat_round(mat2, 1000), h=False)
A = [[5,3,1], [3, 9, 4], [1, 3, 5]]
B = [[1,2,3], [4, 5, 6], [7, 8, 9]]
#print ("mat catcols")
mat2 = mat_catcols(A, B)
#mat_print (mat_round(mat2, 1000), h=False)

Kset = [251, 1, 2, 3, 45, 155, 6, 7]
Mset = [8, 9, 110, 11]
Kdata = [[10, 11, 12, 13, 14, 15, 216, 217]]
print ("data")
mat_print (Kdata, h=True)
matpose = mat_transpose(Kdata)

mat = mat_cauchy (Kset, Mset)
mat2 = mat_catrows(mat_unity(len(Kset)), mat)
print ("encodingmatrix")
mat_print (mat2, h=True)
matcode = mat_mul (mat2, matpose, mulfunc=gf_mul, addfunc=gf_add)
print ("output code")
mat_print (matcode, h=True)
badmatcode = mat_delrows(matcode, 7, 4)
print ("errored code")
mat_print (badmatcode, h=True)
mat3 = mat_delrows(mat2, 7, 4)
print ("encoding matrix after deletion")
mat_print (mat3, h=True)
mat4inv = mat_gf_invert (mat3)
print ("recovery matrix after inversion")
mat_print (mat4inv, h=True)
print ("recovered data after decode")
matdecode = mat_mul (mat4inv, badmatcode, mulfunc=gf_mul, addfunc=gf_add)
matpose = mat_transpose(matdecode)
mat_print (matpose, h=True)
