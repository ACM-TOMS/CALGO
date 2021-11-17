import tabulate, sqlite3, numpy
import numpy as np, math
from sklearn import svm
#import grapher
#https://www.bastibl.net/publication-quality-plots/
import matplotlib as mpl

useLibertine = True
save=1
if save:
    mpl.use("pdf")
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'monospace' : ['Computer Modern']})
#plt.rc('text', usetex=not useLibertine)
if useLibertine:
    pre="""                                                                                                           
\usepackage[T1]{fontenc}                                                                                              
\usepackage[tt=false,type1=true]{libertine}                                                                           
%\setmonofont{inconsolata}                                                                                            
\usepackage[varqu]{zi4}                                                                                               
\usepackage[libertine]{newtxmath}                                                                                     
"""
    plt.rcParams['text.latex.preamble'] = pre #'\usepackage{libertine},\usepackage[libertine]{newtxmath},\usepackage{sfmath},\usepackage[T1]{fontenc}'
else:
    plt.rc('font', **{'family': 'Linux Libertine O' if useLibertine else 'serif',
                      'serif': ['Computer Modern'], 'monospace' : ['Computer Modern']})
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)
plt.rc('axes', labelsize=8)
mylinewidth=0.6
plt.rc('axes',linewidth=mylinewidth)#default 1

width = 3.487
height = width / 1.618
fig, ax = plt.subplots()
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

con = sqlite3.connect("../perf.sqlite")

prep=True #draw prep graph instead of performance

if 0:#sweep d for m=5 instead of sweeping m for d=2 or d=3
    #Two types: log-linear and log-log
    loglog=not prep
    include_esig= not prep
    m=5
    reps=1 if prep else 1000
    pathlength=10
    plotbestfitline=loglog
    if prep:
        a=con.execute("select method,d,time, performed from A where ignore is null and m=?",(m,))
    else:
        a=con.execute("select method,d,time, performed from A where ignore is null and pathlength=? and reps=? and m=?",(pathlength,reps,m))
    a=a.fetchall()
    def x_values(method):
        if "esig"==method:
            return range(2,10)#esig can't do d=10, m=5
        return range(2,11)
    filetag="_level5" + ("loglog" if loglog else "")
    xlab="dimension"
else:
    d=2# 2 or 3
    pathlength=10#100 or 10

    include_esig = (not prep) and ((pathlength==100 and d==2) or (d==3))
    reps = 1000 if pathlength == 100 or d==3 else 10000
    max_m=11 if d==2 else 10
    plotbestfitline=False
    loglog=False
    def x_values(method):
        return range(2,max_m+1)
    if prep:
        a=con.execute("select method,m,time, performed from A where ignore is null and d=?",(d,))
    else:
        a=con.execute("select method,m,time, performed from A where ignore is null and pathlength=? and reps=? and d=?",(pathlength,reps,d))
    a=a.fetchall()
    filetag=("" if pathlength==100 else ("Length"+str(pathlength)+"-"))+str(d)+"d"
    xlab="level"
b=[i for i in a]
#print (tabulate.tabulate(b,headers=["method","m","time","performed"]))

methods=(("AddinLogSig","signature.logSignature"),("AddinSig","signature.signature"),
           ("PythonSig","BG_allpython_signature"),
          ("iisignatureSig","embed.sig_notemplates"),
          ("TemplateSig","embed.sig"),
#          ("TemplateSig1","embed.sig1"),
          ("GeneratedLogSig","embed.logsig"),("LibalgebraLogSig",",embed.logsigfromsig"),("SigToolsSig","stream2sig"),("SigToolsLogSig","stream2logsig"),
         ("esig","stream2logsig"),
         ("iisignatureLogSigC",",iisignature.logsig"),#the comma was accidental
         ("C",",iisignature.logsig"),#the comma was accidental
         ("iisignatureLogSigCH","iisignature.logsigHall"),
         ("iisignatureLogSigO","iisignature.logsigO"),
         ("O","iisignature.logsigO"),
         ("iisignatureLogSigOH","iisignature.logsigOHall"),
         ("iisignatureLogSigS","iisignature.fromsig"),
         ("S","iisignature.fromsig"),
         ("iisignatureLogSigSH","iisignature.fromsigHall"),
         ("prepare('c')","(iisignature.prepare)"),
         ("prepare('s')","iisignature.prepareS"))

def fullname(x,hall=False):
    if prep:
        if hall and (x in ("C", "O","S")):
            return "iisignature.prepare"+x+"Hall"
        if x=="C":
            return "(iisignature.prepare)"
        if x=="O":
            return "iisignature.prepareO"
        if x=="S":
            return "iisignature.prepareS"
        raise Exception("what is "+x)
    for m,n in methods:
        if x==m:
            return n
    print (x)
    raise Exception ("nothing is called "+x+" as output")

def get(method, m):
    bb=[line for line in b if line[0]==method and line[1]==m]
    if len(bb)==0:
        raise Exception ("not found anything for "+method + " "+str(m))
    if len(bb)>1:
        bbb=[line for line in bb if line[3] is not None]
        if len(bbb)==1:
            bb=bbb
        else:
            raise Exception("two copies of "+method+" "+str(m))
    return bb[0][2]
Order = ["C","O","S"]
if include_esig:
    Order.append("esig")

#series = [i for j in [[x,k,"+"] for k in y] for i in j]
#plt.plot(*series)
#d,x,o,8 and all the triangles look bigger , 8 and o look identical
#.,+,* look smaller
for nm,symbol,col in zip(Order,["v","o","x","d"],["r","b","g","k"]):
    if prep and nm =="O":
        continue
    i=[get(fullname(nm),m) for m in x_values(nm)]
    #print(i)
    #plt.plot(x,i,symbol,label=nm)
    nm1=r"\verb|"+nm+"|"
    x=list(x_values(nm))
    plt.scatter(x,i,c=col,marker=symbol,label=nm1,edgecolors='none')
    if prep:
        iHall=[get(fullname(nm,True),m) for m in x]
        plt.scatter(x,iHall,s=8,c=col,marker=symbol,label=nm1+" Hall",edgecolors='none')
    if plotbestfitline:
        def log(x):
            return math.log(float(x))
        x4fit=[ii for ii in x if ii>=6]
        y4fit=i[(-len(x4fit)):]
        lx4fit=list(map(log,x4fit))
        ly4fit=list(map(log,y4fit))
        fit=np.polyfit(lx4fit,ly4fit,1)
        pred4fit=np.poly1d(fit)(lx4fit)
        pred4fit=list(map(math.exp,pred4fit))
        plt.plot(x4fit,pred4fit,c=col,linewidth=0.5,alpha=0.5)
        textheight=y4fit[-1]+(25 if nm=="O" else 0)
        #text="%.2f" % round(fit[0],2)#2dp is too detailed
        text="%.1f" % round(fit[0],1)
        plt.text(11,textheight,text,{'color':col,'fontsize':8})
#prop=mpl.font_manager.FontProperties("monospace")
legend=plt.legend(loc="upper left")
legend.get_frame().set_linewidth(mylinewidth)
plt.yscale('log')
if loglog:
    plt.xscale('log')
    plt.xlabel(xlab + " - logarithmic scale")
    plt.xlim(1,15)
    x_s=list(x_values("C"))
    ax.set_xticks(x_s)
    ax.set_xticklabels([str(i) for i in x_s])
else:
    plt.xlabel(xlab)
#plt.xlim(1.5,max_m+0.5)
plt.ylabel('time(s) - logarithmic scale')
#grapher.pushplot(plt)
#plt.draw()
#dpi doesn't change the scale in the picture
filename = "/home/jeremyr/Dropbox/phd/graphs/"+("prep" if prep else "perf")+"sweeps"+filetag+("Lib" if useLibertine else "")

if save:
    fig.set_size_inches(width,height)
    plt.savefig(filename+".pdf")
else:
    plt.show()
if 0:
    #plt.savefig(filename+".png",dpi=300)
    from PIL import Image
    img = Image.open(filename+".png").convert('LA')
    img.save(filename+"bw.png")

    
