sim.path.Gillespie<-function(x0,tend,S,Rates,thetas) {
    nx=nrow(S); nr=ncol(S)

    tcurr=0
    xcurr=x0
    rates=rep(0,nr)
    path=x0
    tpath=0
    done=FALSE
    nreacts<-rep(0,nr) ## number of each reaction
    
    while (!done) {
        for (i in 1:nr) {
            rates[i]=Rates(i,xcurr,thetas)
        }
        totrate=sum(rates)
        deltat=rexp(1)/totrate
        if (tcurr+deltat<tend) {
            r=sample(1:nr,1,prob=rates/totrate)
            nreacts[r]=nreacts[r]+1
            xcurr=xcurr+S[,r]
            tcurr=tcurr+deltat
        }
        else {
            tcurr=tend
            done=TRUE
        }
        path=rbind(path,xcurr,deparse.level=0)
        tpath=c(tpath,tcurr)
    }
    
    return(list(ts=tpath,xs=path,nreacts=nreacts))
}

sim.path.NextReact<-function(x0,tend,S,Rates,thetas) {
    nx=nrow(S); nr=ncol(S)

    tcurr=0
    xcurr=x0
    rates=rep(0,nr)
    path=x0
    tpath=0
    Lams=rep(0,nr) # cumulative intensity
    Lamscurr=rep(0,nr) 
    Ecurr=rexp(nr)
    
    while (tcurr<tend) {
        for (i in 1:nr) {
            rates[i]=Rates(i,xcurr,thetas)
        }
        ttonext=Ecurr/rates
        rlo=which.min(ttonext)
        tlo=ttonext[rlo]

        if (tcurr+tlo<tend) {
            ## We have an occurrence of reaction rlo, taking time tlo
            xcurr=xcurr+S[,rlo]
            tcurr=tcurr+tlo
            Lamscurr=Lamscurr+rates*tlo
            Lams=rbind(Lams,Lamscurr,deparse.level=0)
            Ecurr=Ecurr-tlo*rates
            Ecurr[rlo]=rexp(1)
            path=rbind(path,xcurr,deparse.level=0)
            tpath=c(tpath,tcurr)
        }
        else {
            path=rbind(path,xcurr,deparse.level=0)
            tpath=c(tpath,tend)
            Lamscurr=Lamscurr+rates*(tend-tcurr)
            Lams=rbind(Lams,Lamscurr,deparse.level=0)
            tcurr=tend+1
        }
    }
    
    return(list(ts=tpath,xs=path,Lams=Lams))
}

## Slower - keeps track of # reactions and the exponentials used
sim.path.NextReact.Detail<-function(x0,tend,S,Rates,thetas,kblock=1000) {
    nx=nrow(S); nr=ncol(S)

    tcurr=0; xcurr=x0
    rates=rep(0,nr)
    path=x0; tpath=0
    totnreact=0
    nreacts=rep(0,nr)
    Nreacts=rep(0,nr) ## number of each reaction
    Lams=rep(0,nr) # cumulative intensity
    Lamscurr=rep(0,nr) 
    kcurr=kblock
    Exps=matrix(nrow=nr,ncol=kblock,data=rexp(nr*kblock))
    Eindices=rep(1,nr)
    Ecurr=Exps[,1]
    
    while (tcurr<tend) {
        for (i in 1:nr) {
            rates[i]=Rates(i,xcurr,thetas)
        }
        ttonext=Ecurr/rates
        rlo=which.min(ttonext)
        tlo=ttonext[rlo]

        if (tcurr+tlo<tend) {
            ## We have an occurrence of reaction rlo, taking time tlo
            xcurr=xcurr+S[,rlo]
            tcurr=tcurr+tlo
            Lamscurr=Lamscurr+rates*tlo
            Lams=rbind(Lams,Lamscurr,deparse.level=0)
            nreacts[rlo]=nreacts[rlo]+1
            Nreacts=rbind(Nreacts,nreacts,deparse.level=0)
            if (Eindices[rlo]+1>kcurr) {
                Expsnew=matrix(nrow=nr,ncol=kblock,data=rexp(nr*kblock))
                Exps=cbind(Exps,Expsnew,deparse.level=0)
                kcurr=kcurr+kblock
            }
            Eindices[rlo]=Eindices[rlo]+1
            Ecurr=Ecurr-tlo*rates
            Ecurr[rlo]=Exps[rlo,Eindices[rlo]]
            path=rbind(path,xcurr,deparse.level=0)
            tpath=c(tpath,tcurr)
        }
        else {
            path=rbind(path,xcurr,deparse.level=0)
            tpath=c(tpath,tend)
            Lamscurr=Lamscurr+rates*(tend-tcurr)
            Lams=rbind(Lams,Lamscurr,deparse.level=0)
            Nreacts=rbind(Nreacts,nreacts,deparse.level=0)
            tcurr=tend+1
        }
    }
    
    return(list(ts=tpath,xs=path,Lams=Lams,Nreacts=Nreacts,Exps=Exps))
}

