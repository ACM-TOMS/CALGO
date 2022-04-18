from manipulateCharacters import loadDataset
import iisignature
import numpy as np
#from sklearn import svm
#import grapher
import sys
import matplotlib as mpl
mpl.use("pdf")
#consider matplotlib2tikz?
import matplotlib.pyplot as plt
useLibertine = True
plt.rc('text', usetex=not useLibertine)
plt.rc('font', **{'family': 'Linux Libertine O' if useLibertine else 'serif',
                  'serif': ['Computer Modern'], 'monospace' : ['Computer Modern']})
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)
plt.rc('axes', labelsize=8)
mylinewidth=0.6
plt.rc('axes',linewidth=mylinewidth)#default 1
#plt.rc('lines',linewidth=mylinewidth)#default 1

width = 3.487
height = width / 1.618
fig, ax = plt.subplots()
fig.subplots_adjust(left=.15, bottom=.20, right=.99, top=.97)

trainX,trainY,testX,testY=loadDataset("/home/jeremyr/data/Handwriting/pendigits.pickle")

s=iisignature.prepare(2,2)

def makeData(x_in,y_in):
    eights=[i[0] for i,j in zip(x_in, y_in) if j==4 and len(i)==1]
    zeros=[i[0] for i,j in zip(x_in, y_in) if j==9 and len(i)==1]
    #sigs8=np.array([iisignature.logsig(i,s) for i in eights])
    #sigs0=np.array([iisignature.logsig(i,s) for i in zeros])
    sigs8=np.array([iisignature.sig(i,5) for i in eights])
    sigs0=np.array([iisignature.sig(i,5) for i in zeros])
    #print (len(trainY),sigs8.shape,sigs0.shape)
    x=np.vstack([sigs8,sigs0])
    n1=len(eights)
    n2=len(zeros)
    y=np.tri(n1+n2,1,-n1,"int8")[:,0] #n1 zeros followed by n2 ones
    return x,y

numbers = [8,0]
#numbers=[6,9]
numberss=[str(i) for i in numbers]

def getAreas(x_in,y_in):
    #eights=[i[0] for i,j in zip(x_in, y_in) if j==8 and len(i)==1]
    #zeros=[i[0] for i,j in zip(x_in, y_in) if j==0 and len(i)==1]
    eights=[i[0] for i,j in zip(x_in, y_in) if j==numbers[0]]
    zeros=[i[0] for i,j in zip(x_in, y_in) if j==numbers[1]]
    sigs8=np.array([iisignature.logsig(i,s) for i in eights])
    sigs0=np.array([iisignature.logsig(i,s) for i in zeros])
    #print me
    return sigs8[:,2],sigs0[:,2]

#if 0:
#    train_x,train_y=makeData(trainX,trainY)
#    clf=svm.SVC()
#    clf.fit(train_x, train_y)
#    test=makeData(testX,testY)
#    print("accuracy:",np.mean(clf.predict(test[0])==test[1]))
#    sys.exit()

areas8, areas0=getAreas(trainX,trainY)

#fig=plt.figure()
plt.hist([areas8,areas0],40,normed=1,label=numberss,color=["k","lightblue"],edgecolor="none")
legend = plt.legend()
legend.get_frame().set_linewidth(mylinewidth)
plt.xlabel('signed area')
plt.ylabel('density')
#ax.tick_params(length=0)
ax.tick_params(direction="out") # looks better with all the bars
#n8,bins8,patches8=plt.hist(areas8,50,normed=1,facecolor="green")
#n0,bins0,patches0=plt.hist(areas0,50,normed=1,facecolor="red")
#grapher.pushplot(plt)
#plt.draw()
#dpi doesn't change the scale in the picture
filename = "/home/jeremyr/Dropbox/phd/graphs/histogram"+numberss[0]+numberss[1]+("Lib" if useLibertine else "")
if 0:
    plt.savefig(filename+".png",dpi=300)
    from PIL import Image
    img = Image.open(filename+".png").convert('LA')
    img.save(filename+"bw.png")

fig.set_size_inches(width,height)
plt.savefig(filename+".pdf")

