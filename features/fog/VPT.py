
import numpy as np
import math
import sys
import random

class Vec:
    def __init__(self,x1=0,y1=0,z1=0):
        self.x = x1
        self.y = y1
        self.z = z1
    
    def copy(self,other):
        return Vec(other.x,other.y,other.z)
    
    def __add__(self,other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Vec(x,y,z)
    
    def __sub__(self,other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Vec(x,y,z)
    
    def mulV(self,other):
        x = self.x*other.x
        y = self.y*other.y
        z = self.x*other.z
        return Vec(x,y,z)
    
    def mulN(self,num):
        return Vec(self.x*num,self.y*num,self.z*num)
    
    def xcross(self,other):
        x = self.y*other.z-self.z*other.y
        y = self.z*other.x-self.x*other.z
        z = self.x*other.y-self.y*other.x
        return Vec(x,y,z)
   
    def dot(self,other):
        return self.x*other.x + self.y*other.y + self.z*other.z
    
    def length(self):
        return ((self.dot(self))**0.5)
    
    def norm(self):
        return self.mulN((1/self.length()))

class Ray:
    def __init__(self,o,d):
        self.o = o
        self.d = d

class Num:
    def __init__(self,num):
        self.num = num
    def copy(self,num):
        return Num(num)

class Sphere:
    def __init__(self,r1,p,e,c,refl):
        self.rad = r1
        self.p = p
        self.e = e
        self.c = c
        #DIFF0,SPEC1,REFR2
        self.refl = refl

    def intersect(self,r,tin,tout):
        op = self.p - r.o
        eps = 0.0001
        b = op.dot(r.d)
        det = b*b - op.dot(op) + self.rad*self.rad
        if det<0:
            return 0
        else:
            det = det**0.5
        if not tin==None and not tout==None:
            if b-det<0:
                tin.num = 0
            else:
                tin.num = b-det
            tout.num = b+det
        if (b-det)>eps:
            return b-det
        elif (b+det)>eps:
            return (b+det)
        else:
            return 0
class Spheres:
    def __init__(self):
        self.size = 0
        self.data = np.empty((128,11))
    
    def getItem(self,key):
        r = self.data[key,0]
        p = Vec(self.data[key,1],self.data[key,2],self.data[key,3])
        e = Vec(self.data[key,4],self.data[key,5],self.data[key,6])
        c = Vec(self.data[key,7],self.data[key,8],self.data[key,9])
        refl = self.data[key,10]
        return Sphere(r,p,e,c,refl)
    
    def add(self,sphere):
        if self.size==128:
            pass
        else:
            self.data[self.size,0] = sphere.rad
            self.data[self.size,1] = sphere.p.x
            self.data[self.size,2] = sphere.p.y
            self.data[self.size,3] = sphere.p.z
            self.data[self.size,4] = sphere.e.x
            self.data[self.size,5] = sphere.e.y
            self.data[self.size,6] = sphere.e.z
            self.data[self.size,7] = sphere.c.x
            self.data[self.size,8] = sphere.c.y
            self.data[self.size,9] = sphere.c.z
            self.data[self.size,10] = sphere.refl
            self.size += 1


spheres = Spheres()
spheres.add(Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(0.0,0.0,0.0),Vec(.75,.25,.25),0))
spheres.add(Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(0.0,0.0,0.0),Vec(.25,.25,.75),0))
spheres.add(Sphere(1e5, Vec(50,40.8, 1e5),     Vec(0.0,0.0,0.0),Vec(.75,.75,.75),0))
spheres.add(Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(0.0,0.0,0.0),Vec(.75,.75,.75),0))
spheres.add(Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(0.0,0.0,0.0),Vec(.75,.75,.75),0))
spheres.add(Sphere(16.5,Vec(27,16.5,47),       Vec(0.0,0.0,0.0),Vec(1,1,1).mulN(0.999), 1))
spheres.add(Sphere(16.5,Vec(73,16.5,78),       Vec(0.0,0.0,0.0),Vec(1,1,1).mulN(0.999), 2))
spheres.add(Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(0.0,0.0,0.0), 0))

homogeneousMedium = Sphere(300, Vec(50,50,80), Vec(0,0,0), Vec(0,0,0), 0);
sigma_s = 0.009
sigma_a = 0.006
sigma_t = sigma_s+sigma_a

def clamp(x):
    if x<0:
        return 0
    elif x>1:
        return 1
    else:
        return x

def toInt(x):
    return int((clamp(x)**(1/2.2))*255 + 0.5)

def sampleSegment( epsilon, sigma, smax):
    return -math.log(1.0 - epsilon * (1.0 - math.exp(-sigma * smax))) / sigma

def sampleSphere(e1,e2):
    z = 1.0 - 2.0 * e1
    sint = sqrt(1.0 - z * z)
    return Vec(math.cos(2.0 * math.pi * e2) * sint, math.sin(2.0 * math.pi * e2) * sint, z);

def sampleHG(g,e1,e2):
    s = 1.0-2.0*e1
    cost = (s + 2.0*g*g*g * (-1.0 + e1) * e1 + g*g*s + 2.0*g*(1.0 - e1+e1*e1))/((1.0+g*s)*(1.0+g*s))
    sint = (1.0-cost*cost)**0.5
    return Vec(math.cos(2.0 * math.pi * e2) * sint, math.sin(2.0 * math.pi * e2) * sint, cost)

def generateOrthoBasis(u,v,w):
    coVec = w;
    if math.fabs(w.x) <= math.fabs(w.y):
        if math.fabs(w.x) <= math.fabs(w.z):
            coVec = Vec(0,-w.z,w.y)
        else:
            coVec = Vec(-w.y,w.x,0)
    elif math.fabs(w.y) <= math.fabs(w.z):
        coVec = Vec(-w.z,0,w.x)
    else:
        coVec = Vec(-w.y,w.x,0)
    coVec = coVec.norm()
    u.x = w.xcross(coVec).x
    u.y = w.xcross(coVec).y
    u.z = w.xcross(coVec).z
    v.x = w.xcross(u).x
    v.y = w.xcross(u).y
    v.z = w.xcross(u).z

def scatter(r,sRay,tin,tout,s):
    s.num = sampleSegment(random.uniform(0,1), sigma_s, tout.num - tin.num)
    x = r.o + r.d.mulN(tin.num) + r.d.mulN(s.num)
    dir = sampleHG(-0.5,random.uniform(0,1), random.uniform(0,1))
    u = Vec(0.0,0.0,0.0)
    v = Vec(0.0,0.0,0.0)
    generateOrthoBasis(u,v,r.d)
    dir = u.mulN(dir.x)+v.mulN(dir.y)+r.d.mulN(dir.z)
    if not sRay==None:
        sRay.o = x
        sRay.d = dir
    return (1.0 - math.exp(-sigma_s * (tout.num - tin.num)))

def intersect(r,t,id,tmax=1e20):
    n=spheres.size
    d=0
    inf=tmax
    t.num=tmax
    i = 0
    while i<n:
        d=spheres.getItem(i).intersect(r,None,None)
        if d and d<t.num:
            t.num=d
            id.num=i
        i += 1
    return t.num<inf

def radiance(r,depth):
    t = Num(0)
    id = Num(0)
    tnear = Num(0)
    tfar = Num(0)
    scaleBy=1.0
    absorption=1.0
    intrsctmd = homogeneousMedium.intersect(r,tnear,tfar) >0
    if intrsctmd:
        sRay = Ray(Vec(0.0,0.0,0.0),Vec(0.0,0.0,0.0))
        s = Num(0)
        ms = scatter(r, sRay, tnear, tfar, s)
        prob_s = ms
        scaleBy = 1.0/(1.0-prob_s)
        if random.uniform(0,1) <= prob_s:
            if not intersect(r, t, id, tnear.num + s.num):
                depth += 1
                return radiance(sRay, depth).mulN(ms*(1.0/prob_s))
            scaleBy = 1.0
        else:
            if not intersect(r, t,id):
                return Vec(0,0,0)
        if t.num >= tnear.num:
            dist = 0
            if t.num>tfar.num:
                dist = tfar.num - tnear.num
            else:
                dist = t.num - tnear.num
            absorption=math.exp(-sigma_t * dist)
    else:
        if not intersect(r, t, id):
            return Vec(0,0,0)
    obj = spheres.getItem(id.num);
    x=r.o+r.d.mulN(t.num)
    n=(x-obj.p).norm()
    nl = Vec(0,0,0)
    if n.dot(r.d)<0:
        nl = n
    else:
        nl = n.mulN(-1)
    f=obj.c
    Le=obj.e
    p = 0
    if f.x>=f.y and f.x>=f.z:
        p = f.x
    elif f.y>=f.x and f.y>=f.z:
        p = f.y
    else:
        p = f.z
    depth += 1
    if depth>5:
        if random.uniform(0,1)<p:
            f=f.mulN(1/p)
        else:
            return Vec(0,0,0)
    if n.dot(nl)>0 or obj.refl != 2:
        f = f.mulN(absorption)
        Le = obj.e.mulN(absorption)
    else:
        scaleBy=1.0
    if obj.refl == 0:
        r1=2*math.pi*random.uniform(0,1)
        r2=random.uniform(0,1)
        r2s=r2**0.5
        w=nl
        u = Vec(0,0,0)
        if math.fabs(w.x)>0.1:
            u = Vec(0,1,0).xcross(w).norm()
        else:
            u = Vec(1,0,0).xcross(w).norm()
        v=w.xcross(u)
        d = (u.mulN(math.cos(r1)*r2s) + v.mulN(math.sin(r1)*r2s)+ w.mulN(((1-r2)**0.5))).norm()
        return (Le + f.mulV(radiance(Ray(x,d),depth))).mulN(scaleBy)

    elif obj.refl == 1:
        return (Le + f.mulV(radiance(Ray(x,r.d-n.mulN(2*n.dot(r.d))),depth))).mulN(scaleBy)
    reflRay = Ray(x, r.d-n.mulN(2*n.dot(r.d)))
    into = n.dot(nl)>0
    nc=1.0
    nt=1.5
    nnt = 0
    if into:
        nnt = nc/nt
    else:
        nnt = nt/nc
    ddn=r.d.dot(nl)
    cos2t = 1-nnt*nnt*(1-ddn*ddn)
    if cos2t<0:
        return (Le + f.mulV(radiance(reflRay,depth)))
    tdir = Vec(0,0,0)
    if into:
        tdir = (r.d.mulN(nnt) - n.mulN((ddn*nnt+(cos2t)**0.5))).norm()
    else:
        tdir = (r.d.mulN(nnt) + n.mulN((ddn*nnt+(cos2t)**0.5))).norm()
    a=nt-nc
    b=nt+nc
    R0=a*a/(b*b)
    c = 0
    if into:
        c = 1+ddn
    else:
        c = 1-tdir.dot(n)
    Re=R0+(1-R0)*c*c*c*c*c
    Tr=1-Re
    P=0.25 + 0.5*Re
    RP=Re/P
    TP=Tr/(1-P)
    if depth>2:
        if random.uniform(0,1)<P:
            return (Le + radiance(reflRay,depth).mulN(RP)).mulN(scaleBy)
        else:
            return (Le + f.mulV(radiance(Ray(x,tdir),depth).mulN(TP))).mulN(scaleBy)
    else:
        return (Le + radiance(reflRay,depth).mulN(Re) + f.mulV(radiance(Ray(x,tdir),depth).mulN(Tr))).mulN(scaleBy)


def draw():
    w = 1024
    h = 768
    samples = 1
    cam = Ray(Vec(50,52,285), Vec(0,-0.042612,-1).norm())
    cx = Vec(w*0.5135/h,0,0)
    cy = (cx.xcross(cam.d).norm()).mulN(0.5135)
    c = np.empty((w*h,3))
    y = 0
    while y<h:
        sys.stdout.write('{}%\r'.format(round(100.0*y/(h-1),2)))
        sys.stdout.flush()
        x = 0
        while x<w:
            sy = 0
            i = (h-y-1)*w+x
            while sy<2:
                sx = 0
                while sx<2:
                    r = Vec(0.0,0.0,0.0)
                    s = 0
                    while s<samples:
                        dx = random.uniform(-1,1)
                        dy = random.uniform(-1,1)
                        d = cx.mulN((((sx + 0.5 + dx)/2 + x)/w - 0.5)) + \
                            cy.mulN((((sy + 0.5 + dy)/2 + y)/h - 0.5)) + \
                            cam.d
                        r = r + radiance(Ray(cam.o+d.mulN(140),d.norm()),0).mulN(1.0/samples)
                        s += 1
                    c[i,0] = c[i,0] + clamp(r.x)*0.25
                    c[i,1] = c[i,1] + clamp(r.y)*0.25
                    c[i,2] = c[i,2] + clamp(r.z)*0.25
                    sx += 1
                sy += 1
            x += 1
        y += 1
    
    fp = open("image.ppm","w")
    fp.write("P3\n1024 768\n255\n")
    i = 0
    wh = w*h
    while i<wh:
        fp.write(str(toInt(c[i,0])))
        fp.write(" ")
        fp.write(str(toInt(c[i,1])))
        fp.write(" ")
        fp.write(str(toInt(c[i,2])))
        fp.write(" ")
        i += 1
    fp.close()

draw()

