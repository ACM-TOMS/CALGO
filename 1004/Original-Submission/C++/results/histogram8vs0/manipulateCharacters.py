#Ben Graham and Jeremy Reizenstein, University of Warwick, GPL v3
import scipy, scipy.interpolate
from copy import deepcopy
from six.moves import cPickle
import six
from six.moves import reduce, range
from math import *
from numpy import *
import numpy

##  __  __             _             _       _         _             _       _                   _       _
## |  \/  |           (_)           | |     | |       | |           (_)     (_)                 | |     | |
## | \  / | __ _ _ __  _ _ __  _   _| | __ _| |_ ___  | |_ _ __ __ _ _ _ __  _ _ __   __ _    __| | __ _| |_ __ _
## | |\/| |/ _` | '_ \| | '_ \| | | | |/ _` | __/ _ \ | __| '__/ _` | | '_ \| | '_ \ / _` |  / _` |/ _` | __/ _` |
## | |  | | (_| | | | | | |_) | |_| | | (_| | ||  __/ | |_| | | (_| | | | | | | | | | (_| | | (_| | (_| | || (_| |
## |_|  |_|\__,_|_| |_|_| .__/ \__,_|_|\__,_|\__\___|  \__|_|  \__,_|_|_| |_|_|_| |_|\__, |  \__,_|\__,_|\__\__,_|
##                      | |                                                           __/ |
##                      |_|                                                          |___/

def boundingBox(c):
    return [reduce(maximum,[s.max(0) for s in c]),\
            reduce(minimum,[s.min(0) for s in c])]

def normalize(strokes,scale=True):
    [mx,mn]=boundingBox(strokes)
    shrinkF=max(max(mx-mn),0.1) if scale else 1.0
    for i,stroke in enumerate(strokes):
        strokes[i]=(2*stroke-mx-mn)/shrinkF
    return strokes

def constantSpeed(path,density):
    # print path.shape,
    lengths=((path[1:,:]-path[:-1,:])**2).sum(1)**0.5
    cs=append(-10**(-30),cumsum(lengths))
    # print cs[-1], 1+int(round(cs[-1]/density))
    if path.shape[0]==1:
        return tile(path,(path.shape[0],1)).astype('float32')
    else:
        f=scipy.interpolate.interp1d(cs, path, axis=0)
        return f(linspace(0,cs[-1],1+int(round(cs[-1]/density)))).astype('float32')
def constantSpeedLookingAtFirstTwoDims(path,density):
    # print path.shape,
    lengths=((path[1:,:2]-path[:-1,:2])**2).sum(1)**0.5
    cs=append(-10**(-30),cumsum(lengths))
    # print cs[-1], 1+int(round(cs[-1]/density))
    if path.shape[0]==1:#only one point
        return tile(path,(path.shape[0],1)).astype('float32')
    else:
        f=scipy.interpolate.interp1d(cs, path, axis=0)
        return f(linspace(0,cs[-1],1+int(round(cs[-1]/density)))).astype('float32')

#as previous, but also give the cos and sin of the direction
def constantSpeedLookingAtFirstTwoDimsWithDir(path,density):
    # print path.shape,
    lengths=((path[1:,:2]-path[:-1,:2])**2).sum(1)**0.5
    cs=append(-10**(-30),cumsum(lengths))
    # print cs[-1], 1+int(round(cs[-1]/density))
    if path.shape[0]==1:#only one point
        return tile(path,(path.shape[0]+2,1)).astype('float32')
    else:
        f=scipy.interpolate.interp1d(cs, path, axis=0)
        out_xs = linspace(0,cs[-1],1+int(round(cs[-1]/density)))
        offsets = numpy.full(out_xs.shape,cs[-1]/(10*out_xs.shape[0]))
        offsets[-1]=-1*offsets[-1]
        out_ys = f(out_xs)#.astype('float32')
        out_perturbed = f(out_xs+offsets)
        grad = out_perturbed-out_ys
        lengths = hypot(grad[:,0],grad[:,1])
        lengths[lengths<0.00000001]=1
        cosAndSin = grad[:,:2]/lengths[:,newaxis]
        return hstack([out_ys,cosAndSin]).astype('float32')

def constantSpeedCharacter(strokes, density):
    return [constantSpeed(s,density) for s in strokes]

def fixedLengthStrokes(path,length):
    # print path.shape,
    lengths=((path[1:,:]-path[:-1,:])**2).sum(1)**0.5
    cs=append(-10**(-30),cumsum(lengths))
    if path.shape[0]==1:
        return tile(path,(length,1)).astype('float32')
    else:
        f=scipy.interpolate.interp1d(cs, path, axis=0)
        return f(linspace(0,cs[-1],length)).astype('float32')

def fixedLengthStrokesCharacter(strokes, length):
    return [fixedLengthStrokes(s,length) for s in strokes]


def distortCharacter(char,stretch=0.3,rotate=0.3,jiggle=0.03,norm=False):
    char=deepcopy(char)
    stretchX=random.uniform(1-stretch,1+stretch)
    stretchY=random.uniform(1-stretch,1+stretch)
    alpha=random.uniform(-rotate,rotate)
    r=random.randint(3)
    for i in range(len(char)):
        char[i][:,0]*=stretchX
        char[i][:,1]*=stretchY
        if r==0:
            char[i]=dot(char[i],[[cos(alpha),sin(alpha)],[-sin(alpha),cos(alpha)]])
        if r==1:
            char[i]=dot(char[i],[[1,0],[alpha,1]])
        if r==2:
            char[i]=dot(char[i],[[1,alpha],[0,1]])
        char[i]+=random.uniform(-jiggle,jiggle,2)
    if norm:
        return normalize(char)
    else:
        return char

def mypickleLoad(f):
    if six.PY3:
        return cPickle.load(f,encoding="bytes")
    else:
        return cPickle.load(f)
    
def loadDataset(filename):
    f=open(filename,"b")
    trainX=[normalize([array(stroke,dtype='float32')\
                       for stroke in character])\
            for character in mypickleLoad(f)]
    trainY=array(cPickle.load(f)).astype('int32')
    testX=[normalize([array(stroke,dtype='float32')\
                       for stroke in character])\
            for character in mypickleLoad(f)]
    testY=array(cPickle.load(f)).astype('int32')
    return trainX,trainY,testX,testY

def pickleSave(name,obj):
    pkl_file = open(name, 'wb')
    cPickle.dump(obj,pkl_file,protocol=-1)
    pkl_file.close()

def pickleLoad(name):
    pkl_file = open(name, 'rb')
    obj=mypickleLoad(pkl_file)
    pkl_file.close()
    return obj
