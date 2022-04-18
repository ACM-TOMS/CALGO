#import signature #BG using swig
#import embed # my purpose-built thing, needs python 2
import numpy, sqlite3, timeit, roughPathSignatures
#, tosig
import esig.tosig as tosig
import iisignature

do10=True

con = sqlite3.connect("perf.sqlite")
con.cursor().execute("create table if not exists A (method text, D INT, M int, reps int, time real, performed TEXT, ignore, version, pathlength)")
con.commit()

if do10:
    default_length = 10
else:
    default_length = 100

def generatePath(d, length = default_length):
    return numpy.random.uniform(size=(length,d)).astype("float32")
def generatePath64(d, length = default_length):
    return numpy.random.uniform(size=(length,d))

#d always means dimension, m always means level

def s1(d,m):#this should perhaps be called sigtools ??
    if d==10:
        return
    z=numpy.zeros((signature.logsigdim(d,m),),dtype="float32")
    signature.logSignature(generatePath(d),m,z)
    return z

def s2(d,m):
    if d==10:
        return
    z=numpy.zeros((signature.sigdim(d,m),),dtype="float32")
    signature.signature(generatePath(d),m,z)
    return z

def s3(d,m):
    #if the length is too long, then this will complain about python recursion
    #100 is fine, 1000 is not
#    if d==3: #takes about 21s per calc
#       return
    return roughPathSignatures.signature(generatePath(d),m)

def s4(d,m):
    return embed.sig_notemplates(generatePath(d),m)

def s5(d,m):
    if(d==2):
        return embed.sig(generatePath(d),m)
    else:
        return embed.sig_general(generatePath(d),m)
def s6(d,m):
    if(d==2):
        return embed.logsig(generatePath(d),m)
    elif(d==3):
        return #I can#t do this case - compilation runs out of memory
    else:
        return embed.logsig_general(generatePath(d),m)
def s7(d,m):
    return embed.logsigfromsig(generatePath(d),m)
def s8(d,m):
    return tosig.stream2sig(generatePath64(d),m)
def s9(d,m):
    return tosig.stream2logsig(generatePath64(d),m)
def s10(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"c")
def s11(d,m):
    iisignature.prepare(d,m,"c")
def s12(d,m):#as of version 0.15, this is really slow, not worth reporting
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"s")
def s13(d,m):
    iisignature.prepare(d,m,"s")
def s14(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"o")
def s15(d,m):
    return iisignature.sig(generatePath(d),m)
def s16(d,m):#as of version 0.15, this is really slow, not worth reporting
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"s")
def s17(d,m):
    iisignature.prepare(d,m,"sh")
def s18(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"c")
def s19(d,m):
    iisignature.prepare(d,m,"ch")
def s20(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"o")
def s21(d,m):
    iisignature.prepare(d,m,"oh")
def s22(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"x")
def s23(d,m):
    global g_s
    return iisignature.logsig(generatePath(d),g_s,"z")
def s24(d,m):
    iisignature.prepare(d,m,"o")
methods = []
#methods=(("s1","signature.logSignature"),("s2","signature.signature"),
#         ("s3","BG_allpython_signature"),
#         ("s4","embed.sig_notemplates"),("s5","embed.sig"),("s6","embed.logsig"),("s7",",embed.logsigfromsig"),
#         ("s8","stream2sig"),("s9","stream2logsig"),("s10",",iisignature.logsig"))
#methods =   (("s3","BG_allpython_signature"),)
#methods = (("s8","stream2sig"),)
#methods.append (("s9","stream2logsig"),)
#methods.extend([("s8","stream2sig"),("s9","stream2logsig")])
#methods = (("s7",",embed.logsigfromsig"),)
#methods.append(("s10",",iisignature.logsig"))
#methods.append(("s11","(iisignature.prepare)"))
#methods.append(("s12","iisignature.fromsig"))
#methods = (("s17","iisignature.prepareSHall"),)
#methods.extend([("s16","iisignature.fromsigHall"),("s17","iisignature.prepareSHall")])
#methods.extend([("s16","iisignature.fromsigHall")])
#methods.extend([("s18","iisignature.logsigHall"),("s19","iisignature.prepareCHall")])
#methods.extend([("s18","iisignature.logsigHall")])
#methods.extend([("s20","iisignature.logsigOHall"),("s21","iisignature.prepareOHall")])
#methods.append(("s21","iisignature.prepareOHall"))
#methods.append(("s19","iisignature.prepareCHall"))
#methods.append(("s17","iisignature.prepareSHall"))
methods.append(("s13","iisignature.prepareS"))
#methods.append(("s24","iisignature.prepareO"))
#methods.extend([("s14","iisignature.logsigO")])
#methods.append(("s22","iisignature.x"))
#methods.extend([("s20","iisignature.logsigOHall")])
#methods.append(("s15","iisignature.sig"))
#methods.append(("s23","logsig_z"))
##methods.extend([("s10",",iisignature.logsig"),("s12","iisignature.fromsig")])
#methods = (("s4","embed.sig_notemplates"),)
#methods = (("s4","embed.sig_notemplates"),("s10",",iisignature.logsig"),("s12","iisignature.fromsig"),("s11","(iisignature.prepare)"),("s13","iisignature.prepareS"))
#methods = (("s10",",iisignature.logsig"),("s12","iisignature.fromsig"),("s14","iisignature.logsigO"),("s15","iisignature.sig"),)

DM_to_test = ((2,2),(2,6), (2,10), (3,6), (3,10),(5,5), (10,4))

#DM_to_test = ((2,2),(2,6), (2,10), (5,5), (10,4))
#DM_to_test = ((10,4),)
#DM_to_test = ((2,2 ),)
DM_to_test = [(2,i) for i in range(2,12)]
DM_to_test = [(3,i) for i in range(2,11)]
#DM_to_test = [(i,5) for i in range(2,11)]
DM_to_test = DM_to_test + [(i,5) for i in range(2,11) if i!=3]
DM_to_test = ((3,10),)
#DM_to_test = [(i,5) for i in [15,20,30,40]]
defaultReps=10 #not enough!

def getPreparation(method, d, m):
    if method in ["s10","s14"]:
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'c');" % (d,m)
    if method=="s12":
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'s');" % (d,m)
    if method=="s16":
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'sh');" % (d,m)
    if method=="s18":
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'ch');" % (d,m)
    if method=="s20":
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'oh');" % (d,m)
    if method in ["s22", "s23"]:
        return "import iisignature; go.g_s = iisignature.prepare(%d,%d,'x');" % (d,m)
    return ""
def isAPrepare(method):
    return method in ["s11", "s13", "s17", "s19", "s21", "s24"]

if __name__=="__main__": # so that when timeit calls import go nothing is run!
    for method,methodname in methods:
        for d,m in DM_to_test:
            reps = (1 if isAPrepare(method) else defaultReps)
            call = "go."+method+"("+str(d)+","+str(m)+")"
            prep = getPreparation(method,d,m)
            t=timeit.timeit(call,"import go; " + prep + call, number=reps)
            print (t, methodname, d, m, reps)
            con.cursor().execute("insert into A (method,d,m,reps,time,pathlength,performed) values (?,?,?,?,?,?,CURRENT_TIMESTAMP)",
                                 (methodname,d,m,reps,t, default_length))
            con.commit()
