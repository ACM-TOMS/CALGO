import tabulate, sqlite3, numpy

con = sqlite3.connect("perf.sqlite")
#con.cursor().execute("create table if not exists A (method text, D INT, M int, reps int, time real)")

length_reps=(10,1000)

a=con.execute("select method,d,m,reps,time,performed,pathlength from A where pathlength=? and reps=? and ignore is null",length_reps).fetchall()
b=[i for i in a]
#print (tabulate.tabulate(b,headers=["method","d","m", "reps","time"]))

methods=(("AddinLogSig","signature.logSignature"),("AddinSig","signature.signature"),
           ("PythonSig","BG_allpython_signature"),
          ("iisignatureSig","embed.sig_notemplates"),
          ("TemplateSig","embed.sig"),
#          ("TemplateSig1","embed.sig1"),
          ("GeneratedLogSig","embed.logsig"),("LibalgebraLogSig",",embed.logsigfromsig"),("SigToolsSig","stream2sig"),("SigToolsLogSig","stream2logsig"),
         ("iisignatureLogSigC",",iisignature.logsig"),#the comma was accidental
         ("iisignatureLogSigCH","iisignature.logsigHall"),
         ("iisignatureLogSigO","iisignature.logsigO"),
         ("iisignatureLogSigOH","iisignature.logsigOHall"),
         ("iisignatureLogSigS","iisignature.fromsig"),
         ("iisignatureLogSigSH","iisignature.fromsigHall"),
         ("prepare('c')","(iisignature.prepare)"),
         ("prepare('s')","iisignature.prepareS"))
def rename(x):
    for m,n in methods:
        if x==n:
            return m #r"\t{"+m+"}"
    raise Exception (x+" has no mapping to a method in the DB")
def fullname(x):
    for m,n in methods:
        if x==m:
            return n
    print (x)
    raise Exception ("nothing is called "+x+" as output")

DM_to_test = ((2,2),(2,6), (2,10), (3,10), (5,5), (10,4))
DM_to_test = ((2,6), (2,10), (3,10), (5,5), (10,4))
DM_to_test = [(2,i) for i in range(2,12)]
DM_to_test = [(3,i) for i in range(2,11)]
#DM_to_test = [(i,5) for i in range(2,11)]
#headers=["","D=2,M=2","D=2,M=6","D=2,M=10","D=3,M=10","D=5,M=5","D=10,M=4"]
headers=["(d,m)","(2,2)","(2,6)","(2,10)","(3,10)","(5,5)","(10,4)"]
headers=["(d,m)","(2,6)","(2,10)","(3,10)","(5,5)","(10,4)"]
headers=["(d,m)"]+[str(i) for i in DM_to_test]
def get(mm):
    a = [mm[0],mm[1]]
    method=mm[2]
    for d,m in DM_to_test:
        #print "d=",d
        found = False
        if (method.startswith("sign") and d==10):# or (d==3 and "embed.logsig"==method):
            a.append(".")
            continue
        if (d==3 and m==2 and "embed.logsig"==method):
            a.append("616.106") # I did this one separately
            continue
        if (d==10 and m==5 and "stream2logsig"==method):
            a.append("") # This one doesn't calculate
            continue
#        for line in b:
#            #print line
#            #print str(b[0]==method),b[1]==d,b[2]==m
#            if line[0]==method and line[1]==d and line[2]==m:
#                a.append("%.3f" % round(line[4],3))
#                #print "ij"
#                if found:
#                    raise Exception("two copies of "+str(d)+" "+method+" "+str(m))
#                found = True
#                break;
#        if not found:
#            raise Exception ("not found anything for "+method + " "+str(d)+" "+str(m))
        bb=[line for line in b if line[0]==method and line[1]==d and line[2]==m]
        if len(bb)==0:
            raise Exception ("not found anything for "+method + " "+str(d)+" "+str(m))
        if len(bb)>1:
            bbb=[line for line in bb if line[5] is not None]
            if len(bbb)==1:
                bb=bbb
            else:
                raise Exception("two copies of "+str(d)+" "+method+" "+str(m))
        a.append("%.2f" % round(bb[0][4],2))

    return a
Order = ["PythonSig","AddinSig","SigToolsSig","TemplateSig",#"TemplateSig1",
         "iisignatureSig",
         "AddinLogSig","SigToolsLogSig","GeneratedLogSig","LibalgebraLogSig",
         "iisignatureLogSigC",
         "iisignatureLogSigCH",
         "iisignatureLogSigO",
         "iisignatureLogSigS",
         "iisignatureLogSigSH",
         "prepare('c')","prepare('s')"
         #""
         ]
Order = ["iisignatureLogSigC",
         "iisignatureLogSigCH",
         "iisignatureLogSigO",
         "iisignatureLogSigOH",
         "iisignatureLogSigS",
         "iisignatureLogSigSH",
         "SigToolsLogSig",
         ]

Ord = [(r"\verb|C|","Lyndon",",iisignature.logsig"),#the comma was accidental
       (r"\verb|C|","standard","iisignature.logsigHall"),
       (r"\verb|O|","Lyndon","iisignature.logsigO"),
       (r"\verb|O|","standard","iisignature.logsigOHall"),
       (r"\verb|S|","Lyndon","iisignature.fromsig"),
       (r"\verb|S|","standard","iisignature.fromsigHall"),
       (r"\verb|esig|","standard","stream2logsig")]


ta = [get(mm) for mm in Ord]
#print (ta)

if 1:#bold the minima
    for i in range(2,len(ta[-1])):
        ind = numpy.argmin([(9999999 if t[i]=="" else float(t[i])) for t in ta])
        ta[ind][i]=r'\bftab '+ta[ind][i]

print (tabulate.tabulate(ta,headers=headers,tablefmt="latex_raw"))

    
