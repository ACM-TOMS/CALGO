import numpy as np
import timeit
import time
import csv,os
import string
import sys
import copy
from six import print_ 
#from collections import Counter

#TODO
#some form of calculation count
#use arbitrary hall basis to see which will lead to the fewest calculations being necessary
#order the words in Polynomial according to length
#retain ratios of integers instead of numbers
#profile this - somewhat done
#template over float vs double to reduce file length, or just use floats.
#don't template over d,m - can just overload
#Write a function to convert a log signature into a signature by expanding the lie brackets
# (e.g. 37[1,2] becomes 37(12)-37(21) ) and exponentiating. and vice versa 
#split bch.cpp into several files to exploit make -j
#try in python 3
#optim: pull lookups in loops into locals outside 

#current calculation times for d=2 only, no float, fn=2
#// 1 Wed Mar 25 17:36:27 2015
#// 2 Wed Mar 25 17:36:27 2015
#// 3 Wed Mar 25 17:36:27 2015
#// 4 Wed Mar 25 17:36:27 2015
#// 5 Wed Mar 25 17:36:27 2015
#// 6 Wed Mar 25 17:36:27 2015
#// 7 Wed Mar 25 17:36:28 2015
#// 8 Wed Mar 25 17:36:28 2015
#// 9 Wed Mar 25 17:36:28 2015
#// 10 Wed Mar 25 17:36:30 2015
#// 11 Wed Mar 25 17:36:42 2015//compiling d=2 m=11 takes ages and much much ram, even w/o -O
#// 12 Wed Mar 25 17:37:46 2015
#// 13 Wed Mar 25 17:43:56 2015
#// 14 Wed Mar 25 18:25:51 2015


#SECTION 1: NUMBERS
# n_add n_mult n_multScalar n_neg n_iszero makenumber n_difference n_one
# are defined which can either be actual floating point numbers (not in use)
# or an algebra of variables

UseNumbers=False
UseSymbols=True

if not UseNumbers:
    if not UseSymbols:
        class Number:
            def __init__(self, a=0):
                self.numb = a
            def __str__(self):
                return "N"+str(self.numb)

        n_one = Number(1)
        def n_add(a,b):
            return Number(a.numb+b.numb)
        def n_mult(a,b):
            return Number(a.numb*b.numb)
        def n_multScalar(a,b):
            return Number(a.numb*b)
        def n_neg(a):
            return Number(-a.numb)
        def n_difference(a,b): #just for reporting
            return a.numb-b.numb
        def n_iszero(a):
            return 0==a.numb
        def makenumber(a): #temporary
            return Number(a)
    else:
        #this class is treated as immutable so that N_polynomial can rely on shallow copies
        class N_monomial:
            def __init__(self, a):#a is a sorted array of strings
                self.strings = a
            def __repr__(self):
                return "M"+str(self.strings)
            def __eq__(self,other):
                return self.strings==other.strings
            def __hash__(self):
                return hash(tuple(self.strings)) #Perhaps we should make a be a tuple instead 
    
        def prodMonomials(a,b):
            ar=sorted(a.strings+b.strings)
            return N_monomial(ar)     
        class N_constantNum:
            def __init__(self,a):
                self.num_=a
        class N_polynomial:
            def __init__(self, **kwargs):
                if "n_word" in kwargs:
                    a=kwargs["n_word"]
                    if isinstance(a,N_monomial):
                        self.terms_ = {a:1}
                    else:
                        raise Exception( "bad initialisation of N_polynomial ", kwargs)
                elif "n_ar" in kwargs:
                    self.terms_ = kwargs["n_ar"]
                elif kwargs:
                    raise Exception( "bad initialisation of N_polynomial ", kwargs)
                else:
                    self.terms_ = {}
            def copy(self):
                #N_monomial is basically immutable so this is safe
                return N_polynomial(n_ar=self.terms_.copy())
            def sumWith(self, b):
                a=b.terms_
                for key in a:
                    if key in self.terms_:
                        f=a[key]+self.terms_[key]
                        if 0==f:
                            del self.terms_[key]
                        else:
                            self.terms_[key]=f
                    else:
                        self.terms_[key]=a[key]
            def scalarMultiply(self, b):
                for key in self.terms_:
                    self.terms_[key]=b*self.terms_[key]
            def __str__(self):
                return str(self.terms_)
            
        n_one=N_constantNum(1)
        def n_add(a,b):
            if isinstance(a,N_polynomial) and isinstance(b,N_polynomial):
                out=a.copy()
                out.sumWith(b)
                return out
            if isinstance(a,N_constantNum) and isinstance(b,N_constantNum):
                return N_constantNum(a.num_+b.num_)
            raise NotImplementedError()
        def n_neg(a):
            if isinstance(a,N_polynomial):
                out=copy.deepcopy(a)
                out.scalarMultiply(-1)
                return out
            if isinstance(a,N_constantNum):
                return N_constantNum(-a.num_)
            raise NotImplementedError()            
        def n_mult(a,b):
            if isinstance(a,N_polynomial) and isinstance(b,N_polynomial):
                total={}
                for keya in a.terms_:
                    for keyb in b.terms_:
                        m = prodMonomials(keya,keyb)
                        coeff=a.terms_[keya]*b.terms_[keyb]
                        if m in total:
                            newc = coeff + total[m]
                            if newc==0:
                                del total[m]
                            else:
                                total[m]=newc
                        else:
                            total[m]=coeff
                return N_polynomial(n_ar=total)
            if isinstance(a,N_constantNum) and isinstance(b,N_polynomial):
#                print "multiplyBy ", a.a
                if(a.num_==1):
                    return b
                out=b.copy()
                out.scalarMultiply(a.num_)
                return out
            if isinstance(a,N_constantNum) and isinstance(b,N_constantNum):
                return N_constantNum(a.num_*b.num_)
            print_(a, b)
            raise NotImplementedError()
        def n_multScalar(a,b):
            if isinstance(a,N_polynomial):
                out=a.copy()
                out.scalarMultiply(b)
                return out
            raise NotImplementedError()            
        def n_iszero(a):
            if isinstance(a,N_polynomial):
                return not bool(a.terms_)
            if isinstance(a,N_constantNum):
                return a.num_==0
            raise NotImplementedError()            
        def n_difference(a,b):
            return n_add(a,n_neg(b))
        g_counter=0
        def makenumber(a):#this ignores a
            global g_counter
            g_counter=g_counter+1
            mon=N_monomial([str(g_counter)])
            return N_polynomial(n_word=mon)
    

else:
    n_one=1
    def n_add(a,b):
        return a+b
    def n_neg(a):
        return -a
    def n_mult(a,b):
        return a*b
    def n_multScalar(a,b):
        return a*b
    def n_iszero(a):
        return a==0
    def n_difference(a,b): #just for reporting
        return a-b
    def makenumber(a):
        return a

#END OF SECTION 1

#SECTION 2: FREE LIE ALGEBRAS
#in this section is defined a Polynomial which represents an element
#of an FLA over the field of numbers defined in section 1
#Note that this is all calculated on the fly as needed - so the size of the alphabet is unknown
#We are using 1 element strings with their own order as our alphabet
#Also a function makeListOfLyndonWords to enumerate a basis of the FLA

class LyndonWord:
    def __init__(self, a, l=None, r=None):
#        self.foliage = str(a)
        self.foliage=str(a) if (type(a)!=type(2) or a<10) else chr(97+a-10)
        self.l = l
        self.r = r
        if len(self.foliage)>1 and l is None:
            raise Exception( "bad initialisation of LyndonWord ",self.foliage)
    def __eq__(self,other):
        return self.foliage==other.foliage
    def __hash__(self):
        return hash(self.foliage)
    def __repr__(self):
        return "LW("+self.foliage+")"
    def __str__(self):
        return self.foliage

#convert a LyndonWord to a string of brackets
def stringBracketedWord(i):
    if(i.l is None):
        return i.foliage
    else:
        return "[" + stringBracketedWord(i.l)+","+stringBracketedWord(i.r)+"]"

class Polynomial:
    def __init__(self, **kwargs):
        if "word" in kwargs:
            a=kwargs["word"]
            if isinstance(a,LyndonWord):
                self.terms = {a:n_one}
            else:
                raise Exception( "bad initialisation of Polynomial ", kwargs)                
        elif "coords" in kwargs:
            self.terms={}
            x=kwargs["coords"]
            for i in range(len(x)):
                l=LyndonWord(i+1)
                self.terms[l]=makenumber(x[i])
        elif "fullhash" in kwargs:
            self.terms=kwargs["fullhash"]
        elif kwargs:
            raise Exception( "bad initialisation of Polynomial ", kwargs)                
        else:
            self.terms = {}
    def __repr__(self):
        return "Polynomial("+str(self.terms)+")"
    def isnonzero(self):
        return bool(self.terms)
    def negate(self):
        for key in self.terms:
           self.terms[key] = n_neg(self.terms[key])
    def scalarMultiplyByObject(self, d):
        for key in self.terms:
            self.terms[key]=n_mult(self.terms[key],d)
    def scalarMultiplyByConstant(self, d):
        if d==0:
            self.terms = {}
        else:
            for key in self.terms:
                self.terms[key]=n_multScalar(self.terms[key],d)
    def sumWith(self, b):
        a=b.terms
        for key in a:
            if key in self.terms:
                f=n_add(a[key],self.terms[key])
                if n_iszero(f):
                    del self.terms[key]
                else:
                    self.terms[key]=f
            else:
#                print "not found:" + str(key)
                self.terms[key]=a[key]

#calculates the lie product of polynomials x and y in the lyndon basis
#ignoring any words created longer than m
#returns an unreferenced polynomial
def productPolynomials(x,y,m):            
    def productLyndonWords(a,b): #returns polynomial which is otherwise unreferenced, or None
        alen=len(a.foliage)
        if len(b.foliage)+alen>m:
            return None#Polynomial()
        if b.foliage==a.foliage:
            return None#Polynomial()
        if b.foliage<a.foliage:
            #print "swapping " + str(a) + " and "+str(b)
            t=productLyndonWords(b,a)
            t.negate()
            return t
        conc = a.foliage+b.foliage
        if conc<b.foliage  and (alen==1 or a.r.foliage>=b.foliage):
            return Polynomial(word=LyndonWord(conc,a,b))
        #print "reordering"
        a1 = productPolynomials(Polynomial(word=a.r),productLyndonWords(b,a.l),m)
        a2 = productPolynomials(Polynomial(word=a.l),productLyndonWords(a.r,b),m)
        a1.sumWith(a2)
        return a1
    out=Polynomial()
    if x is None or y is None:
        return out
    #The following double loop wastes a lot of time
    ykeys=[i for i in y.terms if len(i.foliage)<m]
    for keyx in x.terms:
        if len(keyx.foliage)<m: #for large m, a lot of a general polynomial tends to be terms of degree m, so this check saves looping through y
            for keyy in ykeys: #y.a:
                t = productLyndonWords(keyx, keyy)
                if t is not None:#t.isnonzero():
                    t.scalarMultiplyByObject(n_mult(x.terms[keyx],y.terms[keyy]))
                    out.sumWith(t)
    return out

def makeListOfLyndonWords(d,m,justcount=False):
    a=[[LyndonWord(i+1) for i in range(d)]]
    for M in range(2,m+1):
        mm=[]
        for leftlength in range(1,M):
            rightlength=M-leftlength
            for left in a[leftlength-1]:
                for right in a[rightlength-1]:
                    if left.foliage<right.foliage and (leftlength==1 or left.r.foliage>=right.foliage):
                        mm.append(LyndonWord(left.foliage+right.foliage,left,right))
        a.append(mm)
    if justcount:
        return [len(s) for s in a]
    return [i for s in a for i in s]

#END OF SECTION 2

#SOME SANITY CHECKS
#l=LyndonWord(12,LyndonWord(1),LyndonWord(2))
#L=Polynomial(l)
#L.negate()
#L.sumWith(L)
#print L
#LL=productPolynomials(L,Polynomial(LyndonWord(3)))
#L2=productPolynomials(L,Polynomial(LyndonWord(4)))
#print productPolynomials(LL,L2)

#SECTION 2A: alternatives to the above which use the hall basis in Coropa instead of Lyndon words
#Note that this still represents basis elements using the class LyndonWord even though it wont be a Lyndon Word (TODO: rename LyndonWord)
#set match_coropa to True to actually use this basis

match_coropa=False
def Coropa_LyndonWord_Less(x,y):
    lx=len(x.foliage)
    ly=len(y.foliage)
    if lx==ly:
        if(lx==1):
            return x.foliage<y.foliage #single letter case
        #Comparing x.foliage and y.foliage lexicographically here roughly works - goes wrong e.g. at level 12 in dimension 3
        if(Coropa_LyndonWord_Less(x.l,y.l)):
               return True
        if(Coropa_LyndonWord_Less(y.l,x.l)):
               return False
        return Coropa_LyndonWord_Less(x.r, y.r)
    return lx<ly

def productPolynomialsCoropa(x,y,m):            
    def productLyndonWords(a,b): #returns polynomial which is otherwise unreferenced
        #print "prod("+stringBracketedWord(a)+","+stringBracketedWord(b)+")"
        alen=len(a.foliage)
        blen=len(b.foliage)
        if blen+alen>m:
            return Polynomial()
        if b.foliage==a.foliage:
            return Polynomial()
        if Coropa_LyndonWord_Less(b,a):
            t=productLyndonWords(b,a)
            t.negate()
            return t
        conc = a.foliage+b.foliage
        if blen==1 or not Coropa_LyndonWord_Less(a, b.l):
            return Polynomial(word=LyndonWord(conc,a,b))
        #print "reordering"
        a1 = productPolynomialsCoropa(productLyndonWords(a,b.l),Polynomial(word=b.r),m)
        a2 = productPolynomialsCoropa(productLyndonWords(b.r,a),Polynomial(word=b.l),m)
        a1.sumWith(a2)
        return a1
    out=Polynomial()
    ykeys=[i for i in y.terms if len(i.foliage)<m]
    for keyx in x.terms:
        if len(keyx.foliage)<m: #for large m, a lot of a general polynomial tends to be terms of degree m, so this check saves looping through y
            for keyy in ykeys: #y.a:
                t = productLyndonWords(keyx, keyy)
                if t.isnonzero():
                    t.scalarMultiplyByObject(n_mult(x.terms[keyx],y.terms[keyy]))
                    out.sumWith(t)
    return out

def makeListOfLyndonWordsCoropa(d,m,justcount=False):
    a=[[LyndonWord(i+1) for i in range(d)]]
    for M in range(2,m+1):
        mm=[]
        for leftlength in range(1,1+M//2):
            rightlength=M-leftlength
            for il, left in enumerate(a[leftlength-1]):
                for right in (a[rightlength-1] if leftlength<rightlength else a[rightlength-1][(il+1):]):
                   if rightlength==1 or not Coropa_LyndonWord_Less(left, right.l):
                        mm.append(LyndonWord(left.foliage+right.foliage,left,right))
        a.append(mm)
    if justcount:
        return [len(s) for s in a]
    return [i for s in a for i in s]

def printListOfLyndonWords(d,m,coropa=False):
    l=(makeListOfLyndonWords(d,m) if (not coropa) else makeListOfLyndonWordsCoropa(d,m))
    for i in l:
        print_ (i.foliage, stringBracketedWord(i))
        
def verifyMultiplication(m,product,listOfLyndonWords):
    #amateur check that listOfLyndonWords and productPolynomials are consistent
    #this function should return without printing anything
    lyndonWordHash=dict((i,1) for i in listOfLyndonWords) 
    enumeratedListOfLyndonWords = list(enumerate(listOfLyndonWords))
    poly1=Polynomial(fullhash=dict((l,N_polynomial(n_word=N_monomial(["a["+str(i)+"]"]))) for i,l in enumeratedListOfLyndonWords))
    poly2=Polynomial(fullhash=dict((l,N_polynomial(n_word=N_monomial(["n["+str(i)+"]"]))) for i,l in enumeratedListOfLyndonWords))
    poly3=product(poly1,poly2,m) #in case of inconsistency, expect this to recurse indefinitely
    for key in poly3.terms:
        if not (key in lyndonWordHash):
            print_ ("product doesn't work, it gives keys which are not lyndon words, e.g." + key)
def verifyMultiplicationLyndon(d,m):
    listOfLyndonWords = makeListOfLyndonWords(d,m)
    verifyMultiplication(m,productPolynomials,listOfLyndonWords)
def verifyMultiplicationCoropa(d,m):
    listOfLyndonWords = makeListOfLyndonWordsCoropa(d,m)
    verifyMultiplication(m,productPolynomialsCoropa,listOfLyndonWords)
def verifyUniqueFoliagesCoropa(d,m):
    listOfLyndonWords = makeListOfLyndonWordsCoropa(d,m)
    mapOfFoliages = dict((w.foliage,1) for w in listOfLyndonWords)
    print_ (len(mapOfFoliages), len(listOfLyndonWords))
#printListOfLyndonWords(2,8)
#printListOfLyndonWords(3,12,True)
#print makeListOfLyndonWordsCoropa(5,10,True)
#for i in range(1,6):
#    m=12
#    k = sum(i**M for M in range(1,m+1))
#    print k, sum(makeListOfLyndonWords(i,m,True))
#verifyMultiplicationCoropa(3,12)    
#go up to level 12, even in comparison with coropa    
#verifyUniqueFoliagesCoropa(5,10)
#sys.exit()
if match_coropa:
    def makeListOfLyndonWords(d,m):
        return makeListOfLyndonWordsCoropa(d,m)
    def productPolynomials(a,b,m):
        return productPolynomialsCoropa(a,b,m)
#END OF SECTION 2A

#SECTION 3 - BCH
#a function bch2 is defined which calculated the bch product in the Free Lie Algebra
#this depends crucially on loading the coefficients from a file which has been downloaded from http://www.ehu.eus/ccwmuura/bch.html

def readbchcoords(sourcefile):
    ta = np.loadtxt(sourcefile)#this will result in float type even though we have ints everywhere
    ta[:,3]=ta[:,3]/ta[:,4]#safe - this is division of floats
    ta[:,4]=0
    return ta

if len(sys.argv)==3:
    default_directory_for_bch_data=""
elif os.name=="nt":
    default_directory_for_bch_data="C:/play/MachineLearning/"
else:
    default_directory_for_bch_data="/storage/maths/phrnai/libs/"
directory_for_bch_data=default_directory_for_bch_data #USER CONTROL


ta=readbchcoords(directory_for_bch_data + "bchLyndon20.dat")
#ta1=readbchcoords(directory_for_bch_data + "bchHall20.dat")

#These two functions are not in use - just if you want to construct words from the file.
def makeListOfLyndonWordsFrom_bchHall_file(m):
    ta=readbchcoords(directory_for_bch_data + "bchHall20.dat")
    out=[LyndonWord(1),LyndonWord(2)]
    for bchterm in range(2,totallengths[m-1]):
        left = out[int(ta[bchterm,1])-1]
        right = out[int(ta[bchterm,2])-1]
        #out.append(LyndonWord(left.foliage+right.foliage,left,right))
        out.append(LyndonWord(right.foliage+left.foliage,right,left)) #do this the wrong way round to match coropa
    return out
def printListOfLyndonWordsFrom_bchHall_file(m):
    for i in makeListOfLyndonWordsFrom_bchHall_file(m):
        print_( i.foliage, stringBracketedWord(i))


def counts(x):
    o=np.zeros((20,),np.int_)
    for i in x:
        o[i-1]=1+o[i-1]
    return o

orders = np.zeros(ta.shape[0],np.int_)
orders[0]=orders[1]=1
for i in range(2,ta.shape[0]):
    orders[i]=orders[int(ta[i,1])-1]+orders[int(ta[i,2])-1] #ta is actually floats, but the first 3 columns are integers from the file
levellengths=counts(orders)
totallengths=np.cumsum(levellengths)

#up to level m
def bch(x,y,bracket,sum,multiply,ta,m):
    arr=[x,y]
    for bchterm in range(2,totallengths[m-1]):
        left = int(ta[bchterm,1])
        right = int(ta[bchterm,2])
        arr.append(bracket(arr[left-1],arr[right-1]))
    out=copy.deepcopy(x)#just a dictionary copy would do
    sum(out,y)
    for bchterm in range(2,totallengths[m-1]):
        scale = ta[bchterm,3]
        poly = arr[bchterm]
        multiply(poly,scale)
        sum(out,poly)
    return out

def bch2(x,y,ta,m):
    return bch(x,y,lambda x,y : productPolynomials(x,y,m),lambda x,y : x.sumWith(y),lambda x,y : x.scalarMultiplyByConstant(y),ta,m)
#END OF SECTION 3

#SECTION 4: WRITE CODE AND EXIT
#C++ code is written which calculates the bch of elements which each have components in every lyndon word. 
#NOTE the user control area

#parseCmdLine("38,2") => [38,2]
#parseCmdLine("38") => [38]
#parseCmdLine("2-4,23") => [2,3,4,23]
def parseCmdLine(a):
    out=set()
    for aa in string.split(a,","):
        aaa=string.split(aa,"-")
        if len(aaa)==1:
            out.add(int(aaa[0]))
        else:
            for i in range(int(aaa[0]),int(aaa[1])+1):
                out.add(i)
    return sorted([s for s in out])

tempIndex=0 # this should be inside outputCode, but there is no nonlocal declaration until python 3
def outputCode():
    file=None
    file=open("bch.cpp","w")
    header=None
    header=open("bch.h","w")
    def printContents(s,leftSeg, rightSeg,d,doreversed,inplace):
        global tempIndex
        it = enumeratedListOfLyndonWords
        split=len(listOfLyndonWords)>14
        statementsSinceLastSplit=0
        linesSinceLastSplit=0
        startingTempIndex=tempIndex
        if split:
            print_ ("static void NOINLINE internal"+str(tempIndex)+fn_args+"\n{",file=file)
            tempIndex+=1
        else:
            print_( fn_name+fn_args+"\n{",file=file)
        if doreversed:
            it = reversed(it)
        for i,l in it:
            ii=str(i)
            consts=[]
            if (not inplace) and ((not leftSeg) or i<d):
                consts.append(" a["+ii+"]")
            if(not rightSeg) or i<d:
                consts.append(" b["+ii+"]")
            consts=" +".join(consts)
            skipnewline = 0==len(consts)
            if inplace:
                print_( ("    a["+ii+"] +="+consts),end="",file=file)
            else:
                print_( ("    c["+ii+"] +="+consts),end="",file=file)
            for t in s.terms[l].terms_:
                skip=False
                coeff=s.terms[l].terms_[t]
                if coeff==1 and len(t.strings)==1:#this is one of the obvious terms!
                    skip=True
                elif skipnewline:
                    print_(str(coeff),end="",file=file)
                    skipnewline = False
                elif coeff<0:
                    print_( "\n            "+str(coeff),end="",file=file)
                else:
                    print_( "\n            +"+str(coeff),end="",file=file)
                if (not skip) and split:
                    linesSinceLastSplit += 1
                for tt in t.strings:
                    if skip:
                        pass
                    else:
                        print_("*"+tt, end="",file=file)
            print_(";",file=file)
            statementsSinceLastSplit+=1
            if split and (statementsSinceLastSplit>6 or linesSinceLastSplit>150):
                print_("}\nstatic void NOINLINE internal"+str(tempIndex)+fn_args+"\n{",file=file)
                tempIndex += 1
                statementsSinceLastSplit=0
                linesSinceLastSplit=0
        if split:
            print_("}\n"+fn_name+fn_args+"\n{",file=file)
            if inplace:
                for i in range(startingTempIndex,tempIndex):
                    print_("    internal"+str(i)+"(a,b);",file=file)
            else:
                for i in range(startingTempIndex,tempIndex):
                    print_( "    internal"+str(i)+"(a,b,c);",file=file)
        print_("}", file=file)
    
    sig="LogSignature<d,m>"
    seg="Segment<d>"
    #USER CONTROL AREA
    dimensionsToDo=[2]
    levelsToDo=[1,2,3,4,6]#range(1,10)
    if len(sys.argv)==3:
        dimensionsToDo=parseCmdLine(sys.argv[1])
        levelsToDo=parseCmdLine(sys.argv[2])
    levelsDimensionsToDo=[(m,d) for d in dimensionsToDo for m in levelsToDo]
    doFloat=False #separate function using floats instead of doubles. In fact, floats are usually good enough
    #functionsToDo=frozenset([1,2,3,4,5])
    functionsToDo=frozenset([2]) #i.e., only do the most important function
    #functionsToDo=frozenset([2,5]) #i.e., only do the in place functions
    #USER CONTROL AREA ENDS
    print_( "#ifndef JR_BCH_H\n#define JR_BCH_H",file=header)
    print_( "#include <array>",file=header)
    print_( "template<size_t d, size_t m> class SigLength;",file=header)
    print_( "template<size_t d> using Segment = std::array<double,d>;",file=header)
    if doFloat:
        print_( "template<size_t d> using FSegment = std::array<float,d>;",file=header)
    print_( "template<size_t d, size_t m> using LogSignature=std::array<double,SigLength<d,m>::value>;",file=header)
    if doFloat:
        print_( "template<size_t d, size_t m> using FLogSignature=std::array<float,SigLength<d,m>::value>;",file=header)
    if 1 in functionsToDo:
        print_( "template<size_t d, size_t m>\nvoid joinSegmentToSignature(const "+sig+"& a, const "+seg+"& b, "+sig+"& c);",file=header)
    if 2 in functionsToDo:
        print_( "template<size_t d, size_t m>\nvoid joinSegmentToSignatureInPlace("+sig+"& a, const "+seg+"& b);",file=header)
        if doFloat:
            print_( "template<size_t d, size_t m>\nvoid joinSegmentToSignatureInPlace(F"+sig+"& a, const F"+seg+"& b);",file=header)
    if 3 in functionsToDo:
        print_( "template<size_t d, size_t m>\nvoid joinTwoPaths(const "+seg+"& a, const "+seg+"& b, "+sig+"& c);",file=header)
    if 4 in functionsToDo:
        print_( "template<size_t d, size_t m>\nvoid joinTwoSignatures(const "+sig+"& a, const "+sig+"& b, "+sig+"& c);",file=header)
    if 5 in functionsToDo:
        print_( "template<size_t d, size_t m>\nvoid joinTwoSignaturesInPlace("+sig+"& a, const "+sig+"& b);",file=header)
    print_( '#include "bch.h"',file=file)
    #USER CONTROL - whether you want NOINLINE to do anything will depend on your compiler
    if(os.name=="nt"):
        print_( '#define NOINLINE ',file=file)
    else:
        print_( '#define NOINLINE  __attribute__ ((noinline))' ,file=file)

    for (m,d) in levelsDimensionsToDo:
        print_ ("//",d,m,time.asctime()) #to show progress
        listOfLyndonWords = makeListOfLyndonWords(d,m)
        enumeratedListOfLyndonWords = list(enumerate(listOfLyndonWords))
        arbitrarysignature_a = Polynomial(fullhash=dict((l,N_polynomial(n_word=N_monomial(["a["+str(i)+"]"]))) for i,l in enumeratedListOfLyndonWords))
        arbitrarysignature_b = Polynomial(fullhash=dict((l,N_polynomial(n_word=N_monomial(["b["+str(i)+"]"]))) for i,l in enumeratedListOfLyndonWords))
        linesignature_a = Polynomial(fullhash=dict((LyndonWord(i+1),N_polynomial(n_word=N_monomial(["a["+str(i)+"]"]))) for i in range(d)))
        linesignature_b = Polynomial(fullhash=dict((LyndonWord(i+1),N_polynomial(n_word=N_monomial(["b["+str(i)+"]"]))) for i in range(d)))
        sig="LogSignature<"+str(d)+","+str(m)+">"
        seg="Segment<"+str(d)+">"
        print_("template<> class SigLength<"+str(d)+","+str(m)+">{public: enum {value = "+str(len(listOfLyndonWords))+"}; };",file=header)
        if 1 in functionsToDo:
            fn_name = "template<>\nvoid joinSegmentToSignature<"+str(d)+","+str(m)+">"
            fn_args = "(const "+sig+"& a, const "+seg+"& b, "+sig+"& c)"
            print_( fn_name+fn_args+";",file=header)
            printContents(bch2(arbitrarysignature_a,linesignature_b,ta,m),False,True,d,False,False)
        if 2 in functionsToDo:
            vv = bch2(arbitrarysignature_a,linesignature_b,ta,m)
            fn_name = "template<>\nvoid joinSegmentToSignatureInPlace<"+str(d)+","+str(m)+">"
            fn_args = "("+sig+"& a, const "+seg+"& b)"
            print_( fn_name+fn_args+";",file=header)
            printContents(vv,False,True,d,True,True)
            if doFloat:
                fn_args="(F"+sig+"& a, const F"+seg+"& b)"
                print_( fn_name+fn_args+";",file=header)
                printContents(vv,False,True,d,True,True)
        if 3 in functionsToDo:
            fn_name="template<>\nvoid joinTwoPaths<"+str(d)+","+str(m)+">"
            fn_args="(const "+seg+"& a, const "+seg+"& b, "+sig+"& c)"
            print_( fn_name+fn_args+";",file=header)
            printContents(bch2(linesignature_a,linesignature_b,ta,m),True,True,d,False,False)
        if 4 in functionsToDo:
            fn_name="template<>\nvoid joinTwoSignatures<"+str(d)+","+str(m)+">"
            fn_args="(const "+sig+"& a, const "+sig+"& b, "+sig+"& c)"
            print_( fn_name+fn_args+";",file=header)
            printContents(bch2(arbitrarysignature_a,arbitrarysignature_b,ta,m),False,False,d,False,False)
        if 5 in functionsToDo:
            fn_name="template<>\nvoid joinTwoSignaturesInPlace<"+str(d)+","+str(m)+">"
            fn_args="("+sig+"& a, const "+sig+"& b)"
            print_( fn_name+fn_args+";",file=header)
            printContents(bch2(arbitrarysignature_a,arbitrarysignature_b,ta,m),False,False,d,True,True)
    print_("#endif", file=header)
    print_ ("//done",time.asctime())
outputCode()
sys.exit(0)

#END OF SECTION 4

#SECTION 5 - WRITE MATHEMATICA CODE AND EXIT
#This code makes a mathematica function in bch.m which calculates the 
#level 4 signature of a 2d path made of 4 straight segments
#-this is convenient because the log signature has 8 elements
#matching the 8 degrees of freedom of the input.
#To use this, comment out the outputCode() and the sys.exit(0) just above
def mathematica():
    listOfLyndonWords = makeListOfLyndonWords(2,4)
    #enumeratedListOfLyndonWords = list(enumerate(listOfLyndonWords))
    linesignature_a = Polynomial(fullhash=dict((LyndonWord(i),N_polynomial(n_word=N_monomial(["a[["+str(i)+"]]"]))) for i in range(1,3)))
    linesignature_b = Polynomial(fullhash=dict((LyndonWord(i),N_polynomial(n_word=N_monomial(["b[["+str(i)+"]]"]))) for i in range(1,3)))
    linesignature_c = Polynomial(fullhash=dict((LyndonWord(i),N_polynomial(n_word=N_monomial(["c[["+str(i)+"]]"]))) for i in range(1,3)))
    linesignature_d = Polynomial(fullhash=dict((LyndonWord(i),N_polynomial(n_word=N_monomial(["d[["+str(i)+"]]"]))) for i in range(1,3)))
    s=bch2(linesignature_a,linesignature_b,ta,4)
    s=bch2(s,linesignature_c,ta,4)
    s=bch2(s,linesignature_d,ta,4)
    file=open("bch.m","w")
    print_( "bch[a_,b_,c_,d_] := (",file=file)
    print_("n=ConstantArray[0,8];",file=file)
    for i,l in enumerate(listOfLyndonWords):
        print_( "n[["+str(i+1)+"]]=",file=file,end="")
        first=True
        for t in s.terms[l].terms_:
            if not first:
                print_( "\n       ",file=file,end="")
            coeff=s.terms[l].terms_[t]
            if abs(coeff)<0.00000000001:
                #CHEAT - the tiny results would not appear 
                #if we stuck with exact arithmetic
                #and the 1.3e-17 is not understood by mathematica
                #so just skip these ones.
                continue
            if(coeff<0):
                c=str(coeff)
            else:
                c="+"+str(coeff)
            print_( c,end="",file=file)
            for tt in t.strings:
                 print_( "* "+tt,file=file,end="")
            first=False
        print_( ";",file=file)
    print_( "n)",file=file)

#mathematica()
#sys.exit(0)

#END OF SECTION 5

#SECTION 6 - WRITE PYTHON CODE AND EXIT

def writePythonCode():
    file=open("bch.py","w")
    print_( "#this is bch.py; it is autogenerated by Jeremy's generateBCH code.",file=file)
    print_( "#If you have a 2D path p with n points, which is a numpy array of dimension [n,2],",file=file)
    print_( "#then you can calculate its signature up to level m with bch.getLogSigOfPath(p,m)",file=file)
    print_( "import numpy",file=file)
    levelsToDo=range(1,10) #range(1,6)
    dimensionsToDo=[2];  
    for m in levelsToDo:
        print_ ("//",m,time.asctime()) #to show progress
        for d in dimensionsToDo:
            print_( ("class LogSignature_"+str(d)+"_"+str(m)+":"),file=file)
            print_( "\tdef __init__(self):",file=file)
            listOfLyndonWords = makeListOfLyndonWords(d,m)
            print_( "\t\tself.sigarray=numpy.zeros(",str(len(listOfLyndonWords)),",dtype='float32')",file=file)
            print_( "\tdef __repr__(self):\n\t\treturn str(self.sigarray)",file=file)
            enumeratedListOfLyndonWords = list(enumerate(listOfLyndonWords))
            arbitrarysignature_a = Polynomial(fullhash=dict((l,N_polynomial(n_word=N_monomial(["a["+str(i)+"]"]))) for i,l in enumeratedListOfLyndonWords))
            linesignature_b = Polynomial(fullhash=dict((LyndonWord(i+1),N_polynomial(n_word=N_monomial(["b["+str(i)+"]"]))) for i in range(d)))
            s = bch2(arbitrarysignature_a,linesignature_b,ta,m)
            print_( "\tdef joinSegment(self,b):\n\t\ta=self.sigarray",file=file)
            for i,l in reversed(enumeratedListOfLyndonWords):
                ii=str(i)
                print_( "\t\ta["+ii+"]+=",end="",file=file)
                if i<d:
                    skipnewline = False
                    print_( "b["+ii+"]",end="",file=file)
                else:
                    skipnewline=True
                for t in s.terms[l].terms_:
                    skip=False
                    coeff=s.terms[l].terms_[t]
                    if coeff==1 and len(t.strings)==1:
                        skip=True
                    elif skipnewline:
                        print_( str(coeff),end="",file=file)
                        skipnewline=False
                    elif coeff<0:
                        print_( "\\\n\t\t\t"+str(coeff),end="",file=file)
                    else:
                        print_( "\\\n\t\t\t+"+str(coeff),end="",file=file)
                    for tt in t.strings:
                        if skip:
                            pass
                        else:
                            print_( "*"+tt,end="",file=file)
                print_( "",file=file)
    print_( "def getLogSigObject(d, m):",file=file)
    for m in levelsToDo:
        for d in dimensionsToDo:
            print_( "\tif d=="+str(d)+" and m=="+str(m)+":",file=file)
            print_( "\t\treturn LogSignature_"+str(d)+"_"+str(m),file=file)
    print_( "\traise ValueError('you asked for LogSignature at level '+str(m)+' at dimension '+str(d)+' which is not supported')",file=file)
    print_( "\n#returns Log signature of path p, which should be a numpy array of dimension [n,2] where n>=2, up to level m",file=file)
    print_( "def getLogSigOfPath(p,m):",file=file)
    print_( "\tlogsig=getLogSigObject(2,m)()",file=file)
    print_( "\tlogsig.sigarray[0:2]=p[1,:]-p[0,:]",file=file)
    print_( "\tfor i in range(2,p.shape[0]):",file=file)
    print_( "\t\tlogsig.joinSegment(p[i,:]-p[i-1,:])",file=file)
    print_( "\treturn logsig.sigarray",file=file)
writePythonCode()
sys.exit()
#END OF SECTION 6

#SECTION 7 - MISCELLANY WHICH WON'T HAPPEN - WE'VE EXITED

pathlength=2
np.random.seed(33)
d=2
m=2
def randompath(length = 10):
    a=np.random.uniform(0,1,(length,d))
    a[0,:]=0
    a=np.cumsum(a,0)
    return a;
def randomsteps(length = 10):
    a=np.random.uniform(0,1,(length,d))
    return a;

p=randomsteps(pathlength)
#print p
def stepsig(x):
    return Polynomial(coords=x)

steps = [stepsig(p[i,:]) for i in range(p.shape[0])]

t=timeit.timeit()
sig1= bch2(steps[0],steps[1],ta,m)
sig2= bch2(steps[0],steps[1],ta1,m)
t2= timeit.timeit()
#print "time: ", t2-t
t=t2
for i in range(2,len(steps)):
#    print i
    sig1=bch2(sig1,steps[i],ta,m)
    sig2=bch2(sig2,steps[i],ta1,m)
    t2= timeit.timeit()
    #print "time: ", t2-t
    t=t2

