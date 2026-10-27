BDI.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs
    }
    else if (r==2) {
        rat=thetas[2]*xs
    }
    else if (r==3) {
        rat=thetas[3]
    }
    return (rat)
}
BDI.S=matrix(nrow=1,ncol=3,data=c(+1,-1,+1),byrow=TRUE)


Dimer.Rates<-function(r,xs,thetas) { 
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*(xs[1]-1)/2
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    return (rat)
}
Dimer.S=matrix(nrow=2,ncol=2,data=c(-2,+2,
                                    +1,-1),byrow=TRUE)

LV.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    else if (r==3) {
        rat=thetas[3]*xs[1]*xs[2]
    }
    return (rat)
}
LV.S=matrix(nrow=2,ncol=3,data=c(-1,0,+1,
                                 0,+1,-1),byrow=TRUE)
LV.thetas.Golightly=c(1,.6,.005)
LV.x0.Golightly=c(79,71)
LV.thetas.FGS=c(.6,.3,.001)
LV.x0.FGS=c(40,140)

LVim.Rates<-function(r,xs,thetas) {
    if (r==1) {
        rat=thetas[1]*xs[1]
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    else if (r==3) {
        rat=thetas[3]*xs[1]*xs[2]
    }
    else if (r==4) {
        rat=thetas[4]
    }
    else if (r==5) {
        rat=thetas[5]
    }
    return (rat)
}

LVim.S=matrix(nrow=2,ncol=5,data=c(-1,+0,+1,+1,+0,
                                   +0,+1,-1,+0,+1),byrow=TRUE)
LVim.thetas.bGolightly=c(1,.6,.005,10,10)

SIS.Rates<-function(r,xs,thetas) { ## state = # infected
    rat=-1
    if (r==1) {
        rat=thetas[1]/thetas[3]*xs[1]*(thetas[3]-xs[1])
    }
    else if (r==2) {
        rat=thetas[2]*xs[1]
    }
    return (rat)
}
SIS.S=matrix(nrow=1,ncol=2,data=c(1,-1),byrow=TRUE)
SIS.CSthetas=c(1,.9,1000) ## (beta,gamma,n) in population of size 1000

SIR.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*xs[2]
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    return (rat)
}
SIR.S=matrix(nrow=2,ncol=2,data=c(-1,0,
                                  1,-1),byrow=TRUE)
SIR.CSthetas=c(.01,.25)


SEIR.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*xs[3]
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    else if (r==3) {
        rat=thetas[3]*xs[3]
    }
    return (rat)
}
SEIR.S=matrix(nrow=3,ncol=3,data=c(-1,0,0,
                                   +1,-1,0,
                                   0,+1,-1),byrow=TRUE)
SEIR.thetas.Abak=c(.018,1/12,1/3)

SEIRS.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*xs[3]
    }
    else if (r==2) {
        rat=thetas[2]*xs[2]
    }
    else if (r==3) {
        rat=thetas[3]*xs[3]
    }
    else if (r==4) {
        rat=thetas[4]*(thetas[5]-xs[3]-xs[2]-xs[1])
    }
    return (rat)
}
SEIRS.S=matrix(nrow=3,ncol=4,data=c(-1,0,0,+1,
                                   +1,-1,0,0,
                                   0,+1,-1,0),byrow=TRUE)
SEIRS.thetas=c(1.5/40,1.0,0.3,.03,40)


Schlogel.S=matrix(nrow=1,ncol=4,data=c(+1,-1,+1,-1),byrow=TRUE)
Schlogel.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*(xs[1]-1)/2
    }
    else if (r==2) {
        rat=thetas[2]*xs[1]*(xs[1]-1)*(xs[1]-2)/6
    }
    else if (r==3) {
        rat=thetas[3]
    }
    else if (r==4) {
        rat=thetas[4]*xs[1]
    }
    return (rat)
}
Schlogel.thetas=c(3,.5,.25,2.95)

SchlogelF.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]*xs[1]*(xs[1]-1)/2*xs[2]
    }
    else if (r==2) {
        rat=thetas[2]*xs[1]*(xs[1]-1)*(xs[1]-2)/6
    }
    else if (r==3) {
        rat=thetas[3]*xs[3]
    }
    else if (r==4) {
        rat=thetas[4]*xs[1]
    }
    return (rat)
}
SchlogelF.S=matrix(nrow=3,ncol=4,data=c(+1,-1,+1,-1,
                                       -1,+1,0,0,
                                       0,0,-1,+1),byrow=TRUE)
Schlogel.CSthetas=c(.3,.6,.25,.295)

## x1= RNA, x2= P, x3=P2, x4=DNA.P2, x5=DNA.
Autoreg.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {  ## DNA+P2->DNA.P2
        rat=thetas[1]*xs[5]*xs[3]
    }
    else if (r==2) { ## DNA.P2->DNA+P2
        rat=thetas[2]*xs[4]
    }
    else if (r==3) { ## DNA->DNA+RNA
        rat=thetas[3]*xs[5]
    }
    else if (r==4) { ## RNA->RNA+P
        rat=thetas[4]*xs[1]
    }
    else if (r==5) { ## 2P->P2
        rat=thetas[5]*xs[2]*(xs[2]-1)/2
    }
    else if (r==6) {
        rat=thetas[6]*xs[3] ## P2->2P
    }
    else if (r==7) { ## RNA->0
        rat=thetas[7]*xs[1]
    }
    else if (r==8) { ## P->0
        rat=thetas[8]*xs[2]
    }
    return (rat)
}
Autoreg.S=matrix(nrow=5,ncol=8,data=c(0,0,1,0,0,0,-1,0,
                                      0,0,0,1,-2,2,0,-1,
                                      -1,1,0,0,1,-1,0,0,
                                      1,-1,0,0,0,0,0,0,
                                      -1,1,0,0,0,0,0,0),byrow=TRUE)
thetas.Autoreg.GW2005=c(.1,.7,.35,.2,.1,.9,.3,.1)


StathGir.Rates<-function(r,xs,thetas) {
    rat=-1
    if (r==1) {
        rat=thetas[1]
    }
    else if (r==2) {
        rat=thetas[2]*xs[1]*(xs[1]-1)/2
    }
    else if (r==3) {
        rat=thetas[3]*xs[2]
    }
    else if (r==4) {
        rat=thetas[4]*xs[2]
    }
    return (rat)
}
StathGir.S=matrix(nrow=3,ncol=4,data=c(-1,-2,+2,0,
                                       +0,+1,-1,-1,
                                       +0,+0,+0,1),byrow=TRUE)

Moran.Rates<-function(r,xs,thetas) {
    rat=-1
    alpha=thetas[1]; beta=thetas[2]; u=thetas[3]; v=thetas[4]; npop=thetas[5]
    f=1-xs[1]/npop
    if (r==1) {
        if ((f>0) && (f<=1)) { 
            rat=(1-f)*( alpha*f*(1-u) +  beta*(1-f)*v )
        }
        else {
            rat=0
        }
    }
    if (r==2) {
        if ((f>=0) && (f<1)) { 
            rat=f*( alpha*f*u +  beta*(1-f)*(1-v) )
        }
        else {
            rat=0
        }
    }
    return(rat)
}
Moran.S=matrix(nrow=1,ncol=2,data=c(-1,+1),byrow=TRUE)

