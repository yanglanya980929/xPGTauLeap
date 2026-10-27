
## Input a set of times and the simulated object.
## Output the state at each of the inputted times.
get.state.at.times<-function(times,path) {
    nt=length(times)
    nx=ncol(path$xs)
    npatht=nrow(path$xs) ## number of path times
    Xobs=matrix(nrow=nt,ncol=nx)
    ind=1
    for (i in 1:nt) {
        while ((ind<npatht-1) && (path$ts[ind+1]<times[i])) {
            ind<-ind+1
        }
        Xobs[i,]=path$xs[ind,]
    }
    return(Xobs)
}

## Internal function.
## Get a vector of times and points so that
## plotting gives only horizontal and vertical lines.
get.one<-function(i,a) {
    ts=a$ts; nts=length(ts)
    xs=a$xs[,i]
    tsnew<-rep(0,2*nts-2)
    xsnew=tsnew
    tsnew[1]=ts[1]; tsnew[2*nts-2]=ts[nts]
    xsnew[1]=xs[1]; xsnew[2*nts-2]=xs[nts]
    for (i in 1:(nts-2)) {
        tsnew[2*i]=ts[i+1]
        tsnew[2*i+1]=ts[i+1]
        xsnew[2*i]=xs[i]
        xsnew[2*i+1]=xs[i+1]
    }

    return(list(ts=tsnew,xs=xsnew))
}

plot.path<-function(i,a,lines=FALSE,...) {
    path=get.one(i,a)

    time=path$ts
    x=path$xs
    if (lines) {
        lines(time,x,...)
    }
    else {
        plot(time,x,type="l",...)
    }
}

## Paths of multiple species all on the same window.
plot.paths<-function(is,a,cols,...) {
    lo=min(a$xs[,is]); hi=max(a$xs[,is])
    ni=length(is)
    path=get.one(is[1],a)
    time=path$ts; x=path$xs

    plot(time,x,type="l",ylim=c(lo,hi),col=cols[1],...)
    if (ni>=2) {
        for (i in 2:ni) {
            path=get.one(is[i],a)
            lines(path$ts,path$xs,col=cols[i],...)
        }
    }
}

## Internal function.
## What is the best shape of plot window given the number of plots required?
get.shape<-function(n) {
    if (n==1) {
        return(c(1,1))
    }
    else if (n==2) {
        return(c(2,1))
    }
    else if (n==3) {
        return(c(3,1))
    }
    else if (n==4) {
        return(c(2,2))
    }
    else if (n<=6) {
        return(c(3,2))
    }
    else if (n<=8) {
        return(c(4,2))
    }
    else if (n==9) {
        return(c(3,3))
    }
    else if (n<=12) {
        return(c(4,3))
    }
    else {
        return(c(4,4))
    }
}
## Paths of multiple species, with oe plot window for each.
plot.paths.separate<-function(is,a,...) {
    ni=length(is)
    sh=get.shape(ni)
    par(mfrow=sh)

    for (i in 1:ni) {
        plot.path(is[i],a,...)
    }
}

get.max<-function(i,a) {
    return(max(a$xs[,i]))
}

get.min<-function(i,a) {
    return(min(a$xs[,i]))
}
