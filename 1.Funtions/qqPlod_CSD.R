#' Base de datos de entrada tiene que tener en este orden las
#' columnas SNP, Chromosome, Position. Las demas columnas son
#' los modelos que se hayan evaluado
#' 
#' @Parametros_del_Modelo
#' 
#' @data = Base de datos con las 3 primeras columnas llamadas 
#' SNP, Chromosome, Position en ese orden, y las demas de los modelos a evaluar
#' @modelos = nombre de todos los modelos que se quieren evaluar. Por defecto
#' esta en Nulo y por tanto grafica todos los modelos que aparecen en la
#' base de datos @data
#' @Q1 y @Q3 = Los cuantiles a analizar. Por defecto son el primer (0.025) y (0.95)
#' @escala = Valor en el cual se quiere que los ejes de la grafica incrementen. Por defecto 1.
#' @yim y @xim = Valor para incrementar el limite maximo del eje Y y X. Por defecto es 0.5 y 0.02 unidades
#' @colores = Listado de colores que se quieren usar en el grafico. Esto es para personalizar el grafico 


library("tidyverse")

qqCDS <- function(data = data, 
                  modelos = NULL, ef_fix = NULL, 
                  Q1 = 0.025, Q3 = 0.975, escala=1, 
                  yim = 0.5,xim = 0.02, 
                  colores = NULL) {
  
  lkjf = grep("SNP", colnames(data))
  if (is_empty(lkjf)) {
    stop("\nLa columna del marcador debe llamarse SNP")
  }
  
  lkjf2 = grep("Chromosome", colnames(data))
  if (is_empty(lkjf2)) {
    stop("\nLa columna del marcador debe llamarse Position")
  }
  
  lkjf1 = grep("Position", colnames(data))
  if (is_empty(lkjf1)) {
    stop("\nLa columna del marcador debe llamarse Chromosome")
  }
  
if (!is.null(modelos)) {
  ao = grep(paste(modelos, collapse = "|"), 
            colnames(data))
  ef_fix1 = c(ef_fix, ao)
  rm(ao)
} else { #If is.null(modelos)
  ef_fix1 = ef_fix
}
  n1 <- length(colnames(data)) - length(ef_fix1)
  
  if (n1 == 1) {
    if (is.null(modelos)) {
      ci = colnames(data)[-c(ef_fix)]
    } else {
      ci = modelos
    }
    as = grep(ci, colnames(data))
    d1 <- data[,c(ef_fix,as)]
    d1$Model <- colnames(d1)[grep(ci, colnames(d1))]
    colnames(d1)[grep(ci, colnames(d1))] <- "Scores"
    
    d1 <- subset(d1, !is.na(d1$Scores))
    d1$Pvalue <- 10^(-d1$Scores)
    d1 <- d1[order(d1[,"Pvalue"]),]
    d1$adPvalue <- p.adjust(d1$Pvalue,"fdr")
    d1$p_value_quantiles <- (1:length(d1$Pvalue))/(length(d1$Pvalue)+1)
    d1$log.Quantiles <- -log10(d1$p_value_quantiles)
    d1$Expected <- -log10(ppoints(nrow(d1)))
    
    d2 <- data[,c(ef_fix,as)]
    d2$Modelo <- colnames(d2)[grep(ci, colnames(d2))]
    colnames(d2)[grep(ci, colnames(d2))] <- "Scores"
    d2 <- subset(d2, !is.na(d2$Scores))
    
    d2$Pvalue <- 10^(-d2$Scores)
    
    d2 <- d2[order(d2[,"Pvalue"]),]
    
    d2$p_value_quantiles <- (1:length(d2$Pvalue))/(length(d2$Pvalue)+1)
    d2$log.Quantiles <- -log10(d2$p_value_quantiles)
    d2$Expected <- -log10(ppoints(nrow(d2)))
    
    
    d2 <- dplyr::bind_rows(d2,data.frame(SNP = c("Min1","Max1"), Chromosome = c("Min1","Max1"),
                              Position = c(0,100e12), 
                              Scores = c(0, (max(d2$Expected) + 0.2)),
                              Pvalue = c(1,10^(-(max(d2$Scores) + 0.02))),
                              p_value_quantiles = c(0,0),
                              log.Quantiles = c(0, ((max(d2$log.Quantiles) + 0.2) )),
                              Expected = c(0, ((max(d2$Expected))+0.02)),
                              Modelo = unique(d2$Modelo)))
    d2 <- d2[order(d2[,"Expected"], decreasing = T),]
    N <- length(d2$Scores)
    for (j in 1:N) {
      i <- ceiling((10^-d2$log.Quantiles[j])*N)
      if(i==0)i=1
      d2$c95[j] <- -log10(qbeta(Q3,i,N-i+1))
      d2$c05[j] <- -log10(qbeta(Q1,i,N-i+1) )
    }
    
    ym = ceiling(max(max(d2$c05), max(d1$Scores)))
    
    if (is.null(colores)) {
      clo = c("green3","red","blue","orange2","purple","brown", "cornflowerblue", "darkblue", "brown3",  "brown3", "darkmagenta", "deeppink1 ")
    } else {
        clo = colores
      }
  
    PHG <- ggplot(d1,aes(x = log.Quantiles, y = Scores)) + 
      geom_ribbon(data=d2, aes(x = log.Quantiles,ymin=c95,ymax=c05), fill="grey86") +
      geom_point(aes(col = Model)) + geom_abline(slope = 1) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(0,ym), 
                         breaks = seq(0,ym,by = escala) ) + 
      scale_x_continuous(expand = c(0,0), 
                         limits = c(0,ceiling(max(d1$log.Quantiles))), 
                         breaks = seq(1,ceiling(max(d1$log.Quantiles)), by = escala )) +
      theme_bw() + theme(panel.grid = element_blank(), 
                         legend.text = element_text(size = 15),
                         legend.title = element_blank(),
                         legend.position = "bottom",
                         axis.title = element_text(size = 25),
                         axis.text = element_text(size = 18),
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = 'black')) +
      #guides(fill = guide_legend(nrow = 1)) + 
      labs(x= expression(Expected~~-log[10](italic("P value"))),
           y=expression(Observed~~-log[10](italic("P value"))) ) +
      scale_color_manual(values = c(clo)) +
      scale_linetype_manual(values = c("dashed","dashed")) + 
      coord_cartesian(xlim = c(0,max(d1$log.Quantiles) + xim ), 
                      ylim = c(0, ym + yim))
      
  } else { # n1 > 1

    ####################################################
    ci <- apply(data[,-ef_fix],2,
                function(x) sum(is.na(x)))
    ci <- names(ci[which.min(ci)])
    
    as = grep(ci, colnames(data))
    d1 <- data[,c(ef_fix, as)]
    d1$Model <- colnames(d1)[grep(ci, colnames(d1))]
    colnames(d1)[grep(ci, colnames(d1))] <- "Scores"
    
    d1 <- subset(d1, !is.na(d1$Scores))
    d1$Pvalue <- 10^(-d1$Scores)
    d1 = dplyr::bind_rows(d1,data.frame(SNP = c("Min","Max"), 
                             Chromosome = c("Min","Max"),
                             Position = c(0,1e12), 
                             Scores = c(0, ceiling(max(d1$Scores))),
                             Model = unique(d1$Model),
                             Pvalue = c(1,10^(-ceiling(max(d1$Scores))))))
    
    
    
    d1 <- d1[order(d1[,"Pvalue"]),]
    d1$adPvalue <- p.adjust(d1$Pvalue,"fdr")
    d1$p_value_quantiles <- (1:length(d1$Pvalue))/(length(d1$Pvalue)+1)
    d1$log.Quantiles <- -log10(d1$p_value_quantiles)
    d1$Expected <- -log10(ppoints(nrow(d1)))
    
    
    d2 <- data[,c(ef_fix,as)]
    colnames(d2)[grep(ci, colnames(d2))] <- "Scores"
    d2 <- subset(d2, !is.na(d2$Scores))
    
    d2$Pvalue <- 10^(-d2$Scores)
  
    d2 <- d2[order(d2[,"Pvalue"]),]
    
    d2$p_value_quantiles <- (1:length(d2$Pvalue))/(length(d2$Pvalue)+1)
    d2$log.Quantiles <- -log10(d2$p_value_quantiles)
    d2$Expected <- -log10(ppoints(nrow(d2)))
    d2 <- dplyr::bind_rows(d2,data.frame(SNP = c("Min1","Max1"), 
                              Chromosome = c("Min1","Max1"),
                              Position = c(0,100e12), 
                              #Model = unique(d2$),
                              Scores = c(0, (max(d2$Expected) + 0.2)),
                              Pvalue = c(1,10^(-(max(d2$Scores) + 0.02))),
                              p_value_quantiles = c(0,0),
                              log.Quantiles = c(0, ((max(d2$log.Quantiles) + 0.2) )),
                              Expected = c(0, ((max(d2$Expected))+0.02)) ))
    d2 <- d2[order(d2[,"Expected"], decreasing = T),]
    
    N <- length(d2$Scores)
    for (j in 1:N) {
      i <- ceiling((10^-d2$log.Quantiles[j])*N)
      if(i==0)i=1
      d2$c95[j] <- -log10(qbeta(Q3,i,N-i+1))
      d2$c05[j] <- -log10(qbeta(Q1,i,N-i+1) )
    }
    
    ass = grep(ci, colnames(data))
    nm = colnames(data)[-c(ef_fix, ass)]
    for (i in 1:length(nm)) {
      as1 = grep(paste0("^",nm[i],"$"), colnames(data))
      asd <- data[,c( ef_fix,as1)]
      asd$Model <- colnames(asd)[grep(nm[i], colnames(asd))]
      colnames(asd)[grep(paste0("^",nm[i],"$"), 
                         colnames(asd))] <- "Scores"
      asd <- subset(asd, !is.na(asd$Scores))
      asd$Pvalue <- 10^(-(asd[,"Scores"]))
      asd = dplyr::bind_rows(asd,data.frame(SNP = c("Min","Max"), Chromosome = c("Min","Max"),
                               Position = c(0,1e12), 
                               Scores = c(0, (max(asd$Scores) + 0.02)),
                               Model = unique(asd$Model),
                               Pvalue = c(1,10^(-(max(asd$Scores) + 0.02)))))
      
      asd <- asd[order(asd[,"Pvalue"]),]
      asd$adPvalue <- p.adjust(asd$Pvalue,"fdr")
      
      asd$p_value_quantiles <- (1:length(asd$Pvalue))/(length(asd$Pvalue)+1)
      asd$log.Quantiles <- -log10(asd$p_value_quantiles)
      asd$Expected <- -log10(ppoints(nrow(asd)))
      d1 <- dplyr::bind_rows(d1,asd)
    }
    
    ym = ceiling(max(max(d2$c05), max(d1$Scores)))
    if (is.null(colores)) {
      clo = c("red","blue","orange2","purple","brown", "cornflowerblue", "darkblue", "brown3",  "brown3", "darkmagenta", "deeppink1 " ,"green3")
    } else {
      clo = colores
    }
    PHG <- ggplot(d1,aes(x = log.Quantiles, y = Scores)) + 
      geom_ribbon(data=d2, aes(x = log.Quantiles,ymin=c95,ymax=c05), fill="grey86") +
      geom_point(aes(col = Model)) +geom_abline(slope = 1) +
      scale_y_continuous(expand = c(0,0), 
                         limits = c(0,ym), 
                         breaks = seq(0,ym,by = escala) ) + 
      scale_x_continuous(expand = c(0,0), 
                         limits = c(0,(max(d1$log.Quantiles)+0.2)), 
                         seq(1,ceiling(max(d1$log.Quantiles)), by = escala )) +
      theme_bw() + theme(panel.grid = element_blank(), 
                         legend.text = element_text(size = 15),
                         legend.title = element_blank(),
                         legend.position = "bottom",
                         axis.title = element_text(size = 25),
                         axis.text = element_text(size = 18),
                         legend.background = element_blank(),
                         legend.box.background = element_rect(colour = 'black')) +
      #guides(fill = guide_legend(nrow = 1)) + 
      labs(x= expression(Expected~~-log[10](italic("P value"))),
           y=expression(Observed~~-log[10](italic("P value"))) ) +
      scale_color_manual(values = c(clo)) +
      scale_linetype_manual(values = c("dashed","dashed")) + 
      coord_cartesian(xlim = c(0,(max(d1$log.Quantiles) + 0.02)), 
                      ylim = c(0, ym + yim))
    
    
    
}
  
return(PHG)
}



#qqCDS(scores, modelos = c("general"), Q1 = 0.05, Q3 = 0.95)



#gen <- scores[,c("SNP","5-dom-ref")]
#gen$fdr <- p.adjust(gen$`5-dom-ref`,"fdr")
#sum(gen$fdr <= 0.05, na.rm = T)

#a <- set.threshold(BD_uso_EF_Agua_3variables_dMauricio, method = "Bonferroni")
#q <- get.QTL(a)

#a1 <- set.threshold(BD_uso_EF_Agua_1variables_dMauricio, method = "FDR")
#q1 <- get.QTL(a1, models = c(colnames(scores)[-c(1,2,3,5)]))

#faserf <- subset(q, q$Chrom %in% c(1,7))
#faserf$Genes <- ifelse(faserf$Chrom == 1, (faserf$Position - 67625260)/1000000,(faserf$Position - 12686414)/1000000)

#sum(q$Model == "general")
