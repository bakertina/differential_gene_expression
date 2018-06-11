################################################################################
# TINA BAKER
# add_symb_col_id_df_obj_biomaRt (add_symb_col_id_df_obj_biomaRt.R)
# 01/06/17
################################################################################
# This function will take a df object, with ensembl_gene_id in col "X", get the
# gene symbols and then push a vector of symbols into environment (symb) to be appended
# into df in enviroment.
# This script will take care of order and duplicates

add_symb_col_id_df_obj_biomaRt <- function(obj_in, col_id) {
  require(biomaRt)
  suppressMessages(require(plyr))
  suppressMessages(require(dplyr))
  suppressMessages(require(bannerCommenter))
  message(banner("TINA BAKER:  add_symb_col_id_df_obj_biomaRt"))
  message("NB: The id column must be named X")
  message("Requires internet; \nIf build mart fails then change path to use the saved mart")
  if (missing(col_id)) {
    col_id = "X"
  }
  message("building MART .........................")
  #file <- read.csv(myfile_in, header = FALSE, as.is = TRUE) ## CHECK HEADER
  as.data.frame(obj_in) -> file
  file_ord <- file[order(file[col_id], decreasing = FALSE),] # order on col X
  ids <- as.character(file[, col_id]) #get ids
  unique(ids) -> unique # only search Biomart for the unique entries
  
  ###################### CHANGE PATH ###########################################
  # below will make a new mart, if error use a saved mart
  #load(file = "/media/sf_Shared_Folder_Ubuntu/index_files_mart/mart2.mart") # saved mart location
  #load(file = "example_files/mart.mart") # github mart
  #load(file = "example_files/mart2.mart") # github mart #  # mart2 is 91
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # make new mart
  #---> randomaly gives error sometimes; known bug in early bioMArt release
  ##############################################################################
  # ### SET THE THE QUREY
  gb_filter <-
    getBM(
      attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = unique,
      mart = mart
    )# this needs the internet # taken off hngv_symbl
  gb_filter_sort <- gb_filter[order(gb_filter$ensembl_gene_id), ]
  rownames(gb_filter_sort) <- NULL
  detach("package:biomaRt", unload = TRUE) # unload biomart; problems with dplyr
  
  # ### FIND MISSING
  unique_gb <- subset(gb_filter_sort, !duplicated(gb_filter_sort$ensembl_gene_id))
  unique_bio <- unique(gb_filter_sort$ensembl_gene_id)
  #identical(unique, unique_bio) # check if true
  which(!unique %in% unique_bio) -> missing_ids
  which(!unique_bio %in% unique) -> missing_bio #?
  
  # ### LEFT JOIN FOR MULTIPLE VALUES
  # this section was added with the new biomaRt release as quary the datbase take
  # longer than the previous release. Left join the results back to list to 
  # reduce multiple quares for same gene
  as.data.frame(ids) -> ids_df
  colnames(ids_df)[1] <- "ensembl_gene_id"
  as.character(ids_df$ensembl_gene_id) -> ids_df$ensembl_gene_id
  dplyr::left_join(ids_df, gb_filter_sort, by = "ensembl_gene_id") -> gb_filter_sort_left
  
  # ### LEFT JOIN TO INPUT
  # Make a copy and left join for missing values; reorder and save
  as.data.frame(file[col_id]) -> full_list # copy
  colnames(full_list) <- "ensembl_gene_id" # rename header
  as.character(full_list$ensembl_gene_id) -> full_list$ensembl_gene_id
  plyr::join(full_list, gb_filter_sort_left, by = "ensembl_gene_id", type = "left", match = "first" ) -> gb_filter_full
  #gb_filter_full_sort <- gb_filter_full[order(gb_filter_full$ensembl_gene_id),] # already done
  symb  <<- gb_filter_full$external_gene_name # push
  
  # ### PRINT SOME INFO
  message(paste0("-> Values in: ", length(file[, col_id]), " , " ,
                 "Values unique: ", length(unique), " , " ,
                 "Values converted: ", length(unique_gb$ensembl_gene_id)," , " ,
                 "Values out: ", length(symb)
  ))
  message(paste0("-> Misssing values from File not in BIOMART: ", ids[missing_ids]))
  message(paste0("-> Misssing values from BIOMART but in file: ", ids[missing_bio]))
  message(paste0("-> Misssing values return NA"))
  message(boxup("")) # close
}

################################################################################
######################## EXAMPLE ###############################################
################################################################################
# read.csv("/media/sf_Shared_Folder_Ubuntu/protect/day_2_DEG_rm/ACT_vs_LTB_complete.csv", sep=",") -> obj_in
#
# add_symb_col_id_df_obj(obj_in = obj_in, col_id = "X") # run function
# obj_in$symb <- symb # append to object
# write.csv(obj_in, "/media/sf_Shared_Folder_Ubuntu/protect/day_2_DEG_rm/ACT_vs_LTB_complete_symb.csv")
# # save

################################################################################
######################## EXAMPLE GITHUB ########################################
################################################################################
# read.csv("example_files/ACT_vs_LTB_complete.csv") -> obj_in
#
# add_symb_col_id_df_obj(obj_in = obj_in, col_id = "X") # run function
# obj_in$symb <- symb # append to object
# write.csv(obj_in, "example_files/ACT_vs_LTB_complete_symb.csv")
###################### END #####################################################
