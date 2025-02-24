library(readr)
#BiocManager::install("hdf5r")
library("rhdf5")
library("Seurat")
library("hdf5r")
#h5ls("/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/")

setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/")
seurat_object = Read10X_h5("/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/GSM6433585_092A_filtered_feature_bc_matrix.h5")
seurat_object <- CreateSeuratObject(counts = seurat_object, project = "MyProject")

# get expression values
counts_matrix <- GetAssayData(seurat_object, slot = "counts")
counts_matrix = as.matrix(counts_matrix)
counts_matrix = as.data.frame(counts_matrix)
counts_matrix = cbind(rownames(counts_matrix), counts_matrix)
colnames(counts_matrix)[1] = "V1"

write.table(counts_matrix, file = "test_spatial_GEP.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# get coordinates
positions = read_csv("/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/GSM6433585_092A_tissue_positions_list.csv.gz",
                     col_names = FALSE)
# I assume col 5 and col 6 are the coordinates 
positions = cbind(positions$X1, positions$X5, positions$X6)
positions = as.data.frame(positions)
colnames(positions) = c("SpotID", "row", "col")
rownames(positions) = positions$SpotID
spots = colnames(counts_matrix)[2:ncol(counts_matrix)]
positions = positions[spots,1:3]
write.table(positions, file = "test_spatial_coordinates.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


### another way?
# Paths to your files
#csv_gz_path <- "path/to/your/file.csv.gz"



# Uncompress the CSV
#positions =  read_csv("/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/GSM6433585_092A_tissue_positions_list.csv.gz", col_names = FALSE)
#setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test")
#write_csv(positions, "positions.csv")

#csv_path = "/Users/xuhaifeng/Downloads/spatial_analysis_test/positions.csv"
#h5_path = "/Users/xuhaifeng/Downloads/spatial_analysis_test/spatial_sample/GSM6433585_092A_filtered_feature_bc_matrix.h5"
#setwd("/Users/xuhaifeng/Downloads/spatial_analysis_test/tarball")
#tar("the_tarball.tar.gz", files = c(csv_path, h5_path), compression = "gzip")

# Verify creation (optional)
#list.files(path = dirname(tarball_path))
