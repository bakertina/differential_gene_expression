t_test_wilcox_lfc <- function(gathered, mydir, name, ttest = TRUE, 
                              wilcox = TRUE, denom, numer, scrdir=FALSE, 
                              up_down_save=FALSE, rm_below_zeros=FALSE){
  require(lazyeval)
  suppressMessages(require(dplyr))
  require(broom)
  if (!missing(scrdir)){  source(paste0(scrdir, "add_symb_col_id_df_obj_biomaRt.R"))
  }
  
  print(noquote(strrep("*", 81)))
  if (missing(ttest))
    stop("-----> Need to specify if you want T-test calculations done")
  if (missing(wilcox))
    stop("-----> Need to specify if you want Wilcox calculations done")
  if (missing(denom))
    stop("-----> Need to specify what the denominator is")
  if (missing(numer))
    stop("-----> Need to specify what the numerator is")
  as.factor(gathered$test_variable) -> gathered$test_variable
  if (!length(levels(gathered$test_variable)) ==2)
    stop("Too many levels in test_variable, subset data and/or droplevels./n e.g. df")
  print(paste0("Will run; wilcox:", wilcox, "......T-test:", ttest))
  print(paste0("Using group -> ", denom, " as the denominator (control)"))
  print(paste0("Using group -> ", numer, " as the numerator(treatment) e.g: ", numer,"/", denom))
  if (!missing(scrdir))
    print("Will also add gene symbols to the output files")
  print(noquote(strrep("*", 81)))
  
  ###################################
  ##      T-TEST WUTH BH P.ADJUST
  if(ttest){
    print("Running T-Test......")
    gathered_t_test <- gathered %>% 
      group_by_("X") %>%
      do(tidy(t.test(.$norm_value~ .$test_variable))) %>% # can use my.t.test.p.value
      ungroup() %>%
      mutate(padjust=p.adjust(p.value, method='BH'))
    
    gathered_t_test <- gathered_t_test[order(gathered_t_test$X, decreasing=FALSE), ]  
    if (!missing(scrdir)){
      add_symb_col_id_df_obj_biomaRt(gathered_t_test); gathered_t_test$symb <- symb
    }
    ## SUBSET THE T-TEST RESULTS <0.05
    gathered_t_test  %>%
      base::subset(padjust < 0.05) -> gathered_t_test_0.05
    
    write.csv(gathered_t_test, paste0(mydir, name, "_t_test",".csv"), row.names = F)
    write.csv(gathered_t_test_0.05, paste0(mydir, name, "_t_test_0.05",".csv"), row.names = F)
  }
  ###################################
  ##      WILCOX TEST WITH BH (FDR) ADJUST 
  if (wilcox){
    print("Running Wilcox-Test......")
    suppressWarnings(gathered_wilcox <- gathered %>% ## WARNINGS
      group_by_("X") %>%
      do(tidy(wilcox.test(.$norm_value ~ .$test_variable, data=., paired = FALSE))) %>%
      ungroup() %>%
      mutate(padjust=p.adjust(p.value, method='fdr')))
    
    gathered_wilcox <- gathered_wilcox[order(gathered_wilcox$X, decreasing=FALSE), ] # order the data by X
    if (!missing(scrdir)){
      suppressMessages(add_symb_col_id_df_obj_biomaRt(gathered_wilcox)); gathered_wilcox$symb <- symb
    }
    ##      SUBSET THE WILCOX RESULTS FOR <0.05 P.VALUE
    gathered_wilcox  %>%
      subset(padjust < 0.05) -> gathered_wilcox_0.05
    
    ##      COLLECT THE MEAN VALUES ON THE GROUPS
    dots <- list("test_variable", "X")# dots is to be used in the NSE below 
    gathered  %>%
      group_by_(.dots = dots) %>%
      dplyr::summarise(mean_value=mean(norm_value)) %>%
      spread(test_variable, "mean_value") -> mean_by_metadata
    
    filter_A <- interp(~ f1 > 0, f1 = as.name(denom))# create a filter for NSE  varaible in the metadata
    filter_B <- interp(~ f2 > 0, f2 = as.name(numer))# create a filter for NSE  varaible in the metadata
    #div_vA_vB <- interp(~f1/f2, f1 = as.name(denom), f2 = as.name(numer))# set the NSE for variable # FOR NON LOG DATA
    sub_vA_vB <- interp(~f2-f1, f1 = as.name(denom), f2 = as.name(numer))# subtract the treatment from control
    
    ##      MEAN FOLD CHANGE ALL DATA
    print("Calculating Fold Change all......")
    means_fc <- mean_by_metadata %>% 
      mutate_(log2fc = (sub_vA_vB))
    
    ##      FILTER THE MEAN FC 2 LOGS UP OR DOWN
    means_fc_2fold <- means_fc  %>%
      filter(log2fc>1 | log2fc<(-1))
    
    ##      ADD GENE SYM BEFORE OUTPUT  
    means_fc <- means_fc[order(means_fc$X, decreasing=FALSE), ]
    if (!missing(scrdir)){
      suppressMessages(add_symb_col_id_df_obj_biomaRt(means_fc)); means_fc$symb <- symb
    }
    
    means_fc_2fold <- means_fc_2fold[order(means_fc_2fold$X, decreasing=FALSE), ]
    if (!missing(scrdir)){
      suppressMessages(add_symb_col_id_df_obj_biomaRt(means_fc_2fold)); means_fc_2fold$symb <- symb
    }
    #write.csv(means_fc, paste0(mydir, name, "_mean_log2fc_all", ".csv"), row.names = F)
    #write.csv(means_fc_2fold, paste0(mydir, name, "_mean_log2fc_2fold", ".csv"), row.names = F)
    
    ##      BIND THE WILCOX TEST AND LOGFC TO FILE
    cbind(means_fc, gathered_wilcox[,2:ncol(gathered_wilcox)]) -> wilcox_and_logfc# from col 2 so remove duplicate X
    if (rm_below_zeros) {
      which(wilcox_and_logfc[denom] <0 & wilcox_and_logfc[numer]) -> both_below_zeros
      wilcox_and_logfc[-both_below_zeros,] -> wilcox_and_logfc
    }
    
    ##    ADD THE SIGNIFICANT COLS
    wilcox_and_logfc$Significant <- ifelse(wilcox_and_logfc$padjust < 0.05 & abs(wilcox_and_logfc$log2fc) > 1 , "FDR < 0.05 & 2fc", "Not Sig")
    wilcox_and_logfc$Significant_2 <- ifelse(wilcox_and_logfc$p.value < 0.05 & abs(wilcox_and_logfc$log2fc) > 0.5849625 , "P < 0.05 & 1.5fc", "Not Sig")
    
    write.csv(wilcox_and_logfc, paste0(mydir, name, "_complete", ".csv"), row.names = F)
    
    ##  SUBSET THE SIGNIFICANCE
    wilcox_and_logfc_sub <- wilcox_and_logfc[wilcox_and_logfc$Significant=="FDR < 0.05 & 2fc",]
    wilcox_and_logfc_sub_2 <<- wilcox_and_logfc[wilcox_and_logfc$Significant_2=="P < 0.05 & 1.5fc",]
    
    if (up_down_save){
      flag_up <- which(wilcox_and_logfc_sub$log2fc > 0); flag_down <- which(wilcox_and_logfc_sub$log2fc < 0)
      flag_up2 <- which(wilcox_and_logfc_sub_2$log2fc > 0); flag_down2 <- which(wilcox_and_logfc_sub_2$log2fc < 0)
      res_up <- wilcox_and_logfc_sub[flag_up,]; res_down <- wilcox_and_logfc_sub[flag_down,] # make smaller
      res_up2 <- wilcox_and_logfc_sub_2[flag_up2,]; res_down2 <- wilcox_and_logfc_sub_2[flag_down2,] # make smaller
      write.csv(res_up, paste0(mydir, name, "_FDR_0.05_2fc_UP", ".csv"), row.names = F)
      write.csv(res_down, paste0(mydir, name, "_FDR_0.05_2fc_DOWN", ".csv"), row.names = F)
      write.csv(res_up2, paste0(mydir, name, "_P_0.05_1.5fc_UP", ".csv"), row.names = F)
      write.csv(res_down2, paste0(mydir, name, "_P_0.05_1.5fc_DOWN", ".csv"), row.names = F)
    }
    ###################################
    # ### WRITE RESULTS
    write.csv(wilcox_and_logfc_sub, paste0(mydir, name, "_FDR_0.05_2fc", ".csv"), row.names = F)
    write.csv(wilcox_and_logfc_sub_2, paste0(mydir, name, "_P_0.05_1.5fc", ".csv"), row.names = F)
    
    print(noquote(strrep("*", 81)))
    print(paste0("Number of genes with FDR < 0.05 & 2fc: ", nrow(wilcox_and_logfc_sub)))
    print(paste0("Number of genes with P < 0.05 & 1.5fc: ", nrow(wilcox_and_logfc_sub_2)))
    print("--------> Lots of files saved, see mydir")
    print(noquote(strrep("*", 81)))
  }
}



##      RUN FUNCTION
# t_test_wilcox_lfc(gathered = gathered, 
#                   name =  "HIV_SHORT_VS_NONE", # set a name for the comparison
#                   mydir = mydir_DEG, # set an output directory for files
#                   ttest = FALSE, #  TRUE/FALSE 
#                   wilcox = TRUE,  # TRUE/FALSE
#                   denom ="HIV_negative", # control(NTxC)
#                   numer ="HIV_pos_short_ART", # treatment(condition)
#                   up_down_save = TRUE, # TRUE/FALSE if want the up&down files too
#                   rm_below_zeros = TRUE, # TRUE/FALSE, if the mean of both conditions are below zero, remove
#                   # This can remove #3-5k genes
#                   scrdir = "/media/sf_Shared_Folder_Ubuntu/scripts_r/") 
