
library(bedr)

multiIntersectCNVs <- function(list){
  
  samples = unique(list[[1]]$sample)
  
  df_int=NULL
  
  df_int = lapply(samples, function(sample){
    
    print(sample)
    
    # duplications
    mysample.DUP.sorted = lapply(list, function(program){
      cnvs=program[program$sample==sample & program$cnv_type=="DUP",c(2:4)]
      if(nrow(cnvs)>0){
        bedr.sort.region(cnvs, check.chr = FALSE)
      }
    })
    
    mysample.DUP.sorted <- mysample.DUP.sorted[!(sapply(mysample.DUP.sorted, is.null))]
   
    
     # deletions
    
    mysample.DEL.sorted = lapply(list, function(program){
      cnvs=program[program$sample==sample & program$cnv_type=="DEL",c(2:4)]
      if(nrow(cnvs)>0){
       bedr.sort.region(cnvs, check.chr = FALSE)
      }
    })
    mysample.DEL.sorted <- mysample.DEL.sorted[!(sapply(mysample.DEL.sorted, is.null))]
    
  
  
    
    # multi-intersection
    
    multDUP <- bedr.join.multiple.region(
      x = mysample.DUP.sorted
    );
    multDEL <- bedr.join.multiple.region(
      x = mysample.DEL.sorted
    );
    
    
    # joined df
    
    df_int = rbind(df_int,rbind(data.frame(multDEL[,1:5], type="DEL", sample=sample), data.frame(multDUP[,1:5], type="DUP", sample=sample)))
    
    
    return(df_int)
    
  }
  
  )
  
  
  df_int_rbind = do.call("rbind", df_int)
  
  
}





