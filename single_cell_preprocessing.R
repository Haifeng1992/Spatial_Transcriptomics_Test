library(data.table)
library(readr)
# Define the directory containing the tar files
directory = "/Users/xuhaifeng/Downloads/spatial_analysis_test/single cell sample/"

# List all tar.gz files in the directory
tar_files = list.files(directory, pattern = "\\.tar(\\.gz)?$", full.names = TRUE)

combined_data_list <- list()
combined_cell_type_list <- list()
file_num = length(tar_files)

setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test/single cell sample/")
for (i in 1:file_num) {
  temp_dir = "/Users/xuhaifeng/Downloads/spatial_analysis_test/single cell sample/temp"
  dir.create(temp_dir)
  untar(tar_files[i], exdir = temp_dir)
  all_items = list.files(temp_dir, full.names = TRUE)
  if(length(all_items)!=1){
    print(paste("Something worng with folder", i))
    return()
  }else{
    file_path = all_items[1]
  }
  barcodes_file = file.path(file_path, "count_matrix_barcodes.tsv")
  genes_file = file.path(file_path, "count_matrix_genes.tsv")
  expression_file = file.path(file_path, "count_matrix_sparse.mtx")
  
  cell_type_file = file.path(file_path, "metadata.csv")
  this_cell_type = read.csv(cell_type_file)
  
  this_cell_type = cbind(this_cell_type$X, this_cell_type$celltype_major)
  this_cell_type = as.data.frame(this_cell_type)
  combined_cell_type_list[[i]] = this_cell_type
  
  barcodes = read_tsv(barcodes_file, col_names = FALSE)
  barcodes = as.matrix(barcodes)
  barcodes = as.character(barcodes)
  genes = read_tsv(genes_file, col_names = FALSE)
  genes = as.matrix(genes)
  genes = as.character(genes)
  expression = Matrix::readMM(expression_file)
  expression = as.matrix(expression)
  expression = as.data.frame(expression)
  rownames(expression) = genes
  colnames(expression) = barcodes
  
  combined_data_list[[i]] = expression
  unlink(temp_dir, recursive = TRUE)
  
  cell_type = file.path(file_path, "count_matrix_sparse.mtx")
  
}

final_data = combined_data_list[[1]]
for (i in 2:length(combined_data_list)) {
  final_data = cbind(final_data, combined_data_list[[i]])
}
setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test/")
temp = final_data[1:5,1:5]
final_data = cbind(rownames(final_data), final_data)
colnames(final_data)[1] = "GENES"
write.table(final_data, file = "test_scRNA_GEP.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)
write.table(temp, file = "temp_scRNA_GEP.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

final_celltype = combined_cell_type_list[[1]]
for (i in 2:length(combined_cell_type_list)) {
  final_celltype = rbind(final_celltype, combined_cell_type_list[[i]])
}
colnames(final_celltype) = c("Cell IDs", "CellType")
setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test/")
write.table(final_celltype, file = "test_scRNA_celllabels.txt", sep = "\t", row.names = FALSE, quote = FALSE)
