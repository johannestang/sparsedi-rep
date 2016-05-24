library(forecastexp)
library(macrods)
library(pryr)
library(xtable)

# Load the dataset and define data function
data(LN2010)
datafunc <- partial(getmacrodata, ds=LN2010)

# Dates:
dates <- cbind(rep(1975:2007, each=12), rep(1:12,(2007-1975+1)));
startdate <- fsd <- c(1965,1)

# Define the variables we wish to forecast and the experiment function
vars <- c("ips10", "a0m051", "mtq", "lhnag", "punew", "gmdc", "puxf", "pwfsa")
expfunc <- partial(forecastexp, datafunc=datafunc, startdate=startdate, 
				   dates=dates, vars=vars, verbose=TRUE)

# Define the experiment
runexperiments <- function(h) {
    filename <- paste("SDIResultsH", h, ".RData", sep="");
    if (file.exists(filename)) load(filename) else allres <- list()

    if (!exists("AR", allres)) {
        allres$AR <- expfunc(partial(arforecast, p=0:6), h=h)
        save(allres, file=filename)
    }

	for (ck in c("fixed", "IC1", "IC2", "IC3", "BIC", "CH")) {
		nm <- paste("PC", ck, sep="")
		if (!exists(nm, allres)) {
            allres[[nm]] <- expfunc(partial(pcfactorforecast, p=0:6, k=8, chooseK=ck, factorstartdate=fsd), h=h)
            save(allres, file=filename)
        }
    }

	for (ck in c("fixed", "IC1", "IC2", "IC3", "BIC", "CH")) {
		nm <- paste("SPC", ck, sep="")
		if (!exists(nm, allres)) {
            allres[[nm]] <- expfunc(partial(spcfactorforecast, p=0:6, k=8, chooseK=ck, factorstartdate=fsd), h=h)
            save(allres, file=filename)
        }
    }

	for (ck in c("fixed", "IC1", "IC2", "IC3", "BIC", "CH")) {
		nm <- paste("PostSPC", ck, sep="")
		if (!exists(nm, allres)) {
            allres[[nm]] <- expfunc(partial(postspcfactorforecast, p=0:6, k=8, chooseK=ck, factorstartdate=fsd), h=h)
            save(allres, file=filename)
        }
    }

	for (ck in c("fixed", "BIC")) {
		nm <- paste("CFPC", ck, sep="")
		if (!exists(nm, allres)) {
            allres[[nm]] <- expfunc(partial(cfpcforecast, p=0:6, k=8, chooseK=ck), h=h)
            save(allres, file=filename)
        }
    }

	for (ck in c("fixed", "BIC")) {
		nm <- paste("LARSEN", ck, sep="")
		if (!exists(nm, allres)) {
            allres[[nm]] <- expfunc(partial(larsenforecast, p=0:6, k=8, chooseK=ck), h=h)
            save(allres, file=filename)
        }
    }

    if (!exists("LASSO", allres)) {
        allres$LASSO <- expfunc(partial(lassoforecast, p=0:6), h=h)
        save(allres, file=filename)
    }

    save(allres, file=filename)
}

# Helper function for making the tables
mktable <- function(h, outfilename="")
{
    filename <- paste("SDIResultsH", h, ".RData", sep="");
	load(filename);
	res <- rbind( 
        relmse(allres$PCfixed, allres$AR),
        relmse(allres$PCBIC, allres$AR),
        relmse(allres$PCIC1, allres$AR),
        relmse(allres$PCIC2, allres$AR),
        relmse(allres$PCIC3, allres$AR),
        relmse(allres$PCCH, allres$AR),

        relmse(allres$SPCfixed, allres$AR),
        relmse(allres$SPCBIC, allres$AR),
        relmse(allres$SPCIC1, allres$AR),
        relmse(allres$SPCIC2, allres$AR),
        relmse(allres$SPCIC3, allres$AR),
        relmse(allres$SPCCH, allres$AR),

        relmse(allres$PostSPCfixed, allres$AR),
        relmse(allres$PostSPCBIC, allres$AR),
        relmse(allres$PostSPCIC1, allres$AR),
        relmse(allres$PostSPCIC2, allres$AR),
        relmse(allres$PostSPCIC3, allres$AR),
        relmse(allres$PostSPCCH, allres$AR),

        relmse(allres$CFPCfixed, allres$AR),
        relmse(allres$CFPCBIC, allres$AR),
        relmse(allres$LARSENfixed, allres$AR),
        relmse(allres$LARSENBIC, allres$AR),
        relmse(allres$LASSO, allres$AR),
        
        sqrt(predmse(allres$AR))
    )

	rownames(res) <- paste(c(
            "PC&8", "PC&BIC", "PC&IC$_1$", "PC&IC$_2$", "PC&IC$_3$", "PC&CH",
			"SPC&8", "SPC&BIC", "SPC&IC$_1$", "SPC&IC$_2$", "SPC&IC$_3$", "SPC&CH",
			"Post-SPC&8", "Post-SPC&BIC", "Post-SPC&IC$_1$", "Post-SPC&IC$_2$", "Post-SPC&IC$_3$", "Post-SPC&CH",
            "CFPC&8", "CFPC&BIC", "LARS-EN&8", "LARS-EN&BIC", "LASSO&",
            "\\multicolumn{2}{l}{RMSFE(AR)}"), "&", h, sep="");
	res2 <- formatC(res, digits=4, format="f");
	for (i in 1:ncol(res)) for (j in 1:(nrow(res)-1))
	{
		if (res[j,i] == min(res[-nrow(res),i])) res2[j,i] <- paste("\\textbf{", res2[j,i], "}", sep=""); 
	}
	for (i in 1:ncol(res)) for (j in 1:6)
	{
		if (res[j,i] == min(res[1:6,i])) res2[j,i] <- paste("\\uline{", res2[j,i], "}", sep=""); 
	}
	for (i in 1:ncol(res)) for (j in 7:12)
	{
		if (res[j,i] == min(res[7:12,i])) res2[j,i] <- paste("\\uline{", res2[j,i], "}", sep=""); 
	}
	for (i in 1:ncol(res)) for (j in 13:18)
	{
		if (res[j,i] == min(res[13:18,i])) res2[j,i] <- paste("\\uline{", res2[j,i], "}", sep=""); 
	}

	print(xtable(res2), hline.after=c(6,12,18,23), include.rownames=TRUE, include.colnames=FALSE, 
		  only.contents=TRUE, sanitize.text.function=function(x) {return(x);}, file=outfilename);
}

getPandK <- function(res)
{
    n1 <- length(res)
    n2 <- length(res[[1]]$forecasts)
    ps <- numeric(n2)
    ks <- numeric(n2)

    for (i in 1:n1)
    {
        ps <- ps + res[[i]]$p
        ks <- ks + res[[i]]$k
    }
    ps <- ps/n1
    ks <- ks/n1
    if (length(ks)==0) ks <- ps*NA
    return(c(rbind(ps, ks)))
}

mktablePK <- function(h, outfilename="")
{
    filename <- paste("SDIResultsH", h, ".RData", sep="");
	load(filename);
	res <- rbind( 
        getPandK(allres$PCfixed),
        getPandK(allres$PCBIC),
        getPandK(allres$PCIC1),
        getPandK(allres$PCIC2),
        getPandK(allres$PCIC3),
        getPandK(allres$PCCH),

        getPandK(allres$SPCfixed),
        getPandK(allres$SPCBIC),
        getPandK(allres$SPCIC1),
        getPandK(allres$SPCIC2),
        getPandK(allres$SPCIC3),
        getPandK(allres$SPCCH),

        getPandK(allres$PostSPCfixed),
        getPandK(allres$PostSPCBIC),
        getPandK(allres$PostSPCIC1),
        getPandK(allres$PostSPCIC2),
        getPandK(allres$PostSPCIC3),
        getPandK(allres$PostSPCCH),

        getPandK(allres$CFPCfixed),
        getPandK(allres$CFPCBIC),
        getPandK(allres$LARSENfixed),
        getPandK(allres$LARSENBIC),
        getPandK(allres$LASSO),

        getPandK(allres$AR))

	rownames(res) <- paste(c(
            "PC&8", "PC&BIC", "PC&IC$_1$", "PC&IC$_2$", "PC&IC$_3$", "PC&CH",
			"SPC&8", "SPC&BIC", "SPC&IC$_1$", "SPC&IC$_2$", "SPC&IC$_3$", "SPC&CH",
			"Post-SPC&8", "Post-SPC&BIC", "Post-SPC&IC$_1$", "Post-SPC&IC$_2$", "Post-SPC&IC$_3$", "Post-SPC&CH",
            "CFPC&8", "CFPC&BIC", "LARS-EN&8", "LARS-EN&BIC", "LASSO&",
            "\\multicolumn{2}{l}{AR}"), "&", h, sep="");
	res2 <- formatC(res, digits=2, format="f");
    res2[res2==" NA"] <- ""
	print(xtable(res2), hline.after=c(6,12,18,23), include.rownames=TRUE, include.colnames=FALSE, 
		  only.contents=TRUE, sanitize.text.function=function(x) {return(x);}, file=outfilename);
}

getPsi <- function(h)
{
    psi <- c()
	fdates <- head(dates, nrow(dates)-h)
    for (i in 1:nrow(fdates))
    { 
	    cname <- forecastexp:::.cachefilename("spcfactors", 8, fdates[i,], startdate, "mean", "sd")
	    factorres <- forecastexp:::.getfromcache(cname, "./rescache", NULL)
        psi <- c(psi, factorres$Psi)
    }
    psi
}

runexperiments(12)
runexperiments(6)
runexperiments(24)

if (!file.exists("tex")) dir.create("tex")
mktable(12, outfilename="tex/table12.tex")
mktable(6, outfilename="tex/table6.tex")
mktable(24, outfilename="tex/table24.tex")

mktablePK(12, outfilename="tex/pktable12.tex")
mktablePK(6, outfilename="tex/pktable6.tex")
mktablePK(24, outfilename="tex/pktable24.tex")

cat(round(mean(getPsi(12)), 2), file="tex/psi12.tex")
cat(round(mean(getPsi(6)), 2), file="tex/psi6.tex")
cat(round(mean(getPsi(24)), 2), file="tex/psi24.tex")

