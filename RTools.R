##This is file contains regularly used R functions that can be reused.
##


##Testing function
testing = function(){
    print("You can access RTools functions!")
}

#Function to run the final steps of Differential Gene Expression analysis
runDEG = function (expr, design, contrasts, pval=0.05, fc=1.5, padjst= "fdr") {
    
    vfit = lmFit(expr, design)
    vfit = contrasts.fit(vfit, contrasts = contrasts)
    efit = eBayes(vfit)
    
    show(summary(decideTests(efit, adjust.method = padjst, p.value = pval, lfc = log2(fc))))
    plotSA(efit, main="Final model: Mean-variance trend")
    
    return(efit)
}


#Function to extract summary stats for a specific comparison, and merge it with annotation file
extractStats=function(efit, anno, mergeKey="ID", coef, pval=0.05, fc=1.5, padjst= "fdr", outpath='./',outname=''){
    
    stats=topTable(efit, p.value = pval, lfc = log2(fc), coef=coef, n=Inf)
    
    stats[mergeKey]=rownames(stats)
    stats=merge(stats,anno,by=mergeKey)
    
    write.table(stats, file=paste0(outpath,outname), sep="\t", quote=F, row.names=F)
    head(stats)
    return(stats)
}
#Quickly perform a heatmap
heatmap=function(dt, x, y, path, outname){
    library(ggplot2)
    tb=data.frame(table(dt[,x],dt[,y]))
    #names(tb)=c("Severity","Cluster","Freq")
    ggplot(tb, aes(Var1, Var2, fill = Freq)) + geom_tile(aes(fill = Freq)) + xlab(x) +ylab(y)+
        geom_text(aes(label = round(Freq, 1)), size=10) +
        scale_fill_gradient(low = "grey", high = "red") +theme(text = element_text(size = 18, face='bold'))+ theme(legend.position = "none")
    ggsave(paste0(path,outname,".png"), width = 9, height = 9)
}

##Quick bar or box plot
quickBarBox = function(dt, pltType, xvar, yvar, splitBy="", fillBy="", plotNewBy="", addDots=FALSE, addError=FALSE, se="", outpath="./", prefix=""){
    
    pltList=unique(dt[,plotNewBy])
    #print(pltList)
    for (plt in pltList){
        #print(paste0("Processing ",plt))
        dtsub=dt[dt[,plotNewBy]==plt,]
        
        #print("Number of row is: ", nrow(dtsub))
        p=ggplot(dtsub, aes_string(x=xvar, y=yvar, fill=fillBy)) + 
        facet_wrap(~dtsub[,splitBy]) +
        {if(addDots) geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=.1)} +
        #{if(addError) geom_errorbar(aes_string(ymin=yvar-se, ymax=yvar+se), width=.2, position=position_dodge(.9))} +
        ggtitle(plt)
        
        if(pltType=="bar"){
            p=p + geom_bar(stat="identity", color=FALSE, position=position_dodge())
        }else if(pltType=="box"){
            p=p + geom_boxplot(position=position_dodge(1))
        }
        
        # + geom_errorbar(aes(ymin=RQ-sd, ymax=RQ+sd), width=.2, position=position_dodge(.9)) 
        ggsave(plot = p, paste0(outpath,prefix,"_",plt,".png"))
    }
    
}

##Quick scatter plot with labels
quickScatter=function(data, xvar, yvar, label="", labelSelect={}, title="", xcut=1, ycut=1, outpath="./", outname="scatter"){
    options(ggrepel.max.overlaps = Inf)
    p=ggscatter(data, x=xvar,y=yvar, label=label, title = title,
          color = "steelblue4", size=1, font.label = c(6, "bold","darkred"), 
          label.select = labelSelect, repel=TRUE, max.overlaps=TRUE
          )+geom_vline(xintercept = c(-xcut,xcut),linetype='dotted',col = 'purple1')+geom_hline(yintercept = ycut, linetype='dotted', col = 'purple1')+geom_point(alpha = 0.5)+theme(text = element_text(size=14))
    ggsave(plot=p, paste0(outpath, outname, ".png"))
    show(p)
}

##Calculate q values for eQTL mapping p values
calcQvals=function(dt, pval){
    library(qvalue)
    
    dt$qvalue<-NA
    pval=dt[,pval]
    qval=qvalue(p=pval)
    dt$qvalue=qval$qvalues
    
    return(dt)
    
}


##Write tsv file
write.tsv=function(dt, outpath, outname){
    
    write.table(dt, paste0(outpath,outname,".tsv"), quote=F, row.names=F, col.names=T, sep="\t")
    
}



