#!/usr/bin/env Rscript
library(optparse)
library(graphics)
library(stringr)

tree_width=1000
tree_heat_gap=50
sample_width=100
heat_cluster_gap = 50
cluster_width=50
cluster_text_gap = 20
text_width = 200
white_space = 50
gene_height = 50

paper_size_height = 11.7
paper_size_width = 8.3
# 主函数
main_function <- function(heat_file, tree_file, cluster_file, gff_file,
                          sample_color_file, region, output_prefix) {
    # 在这里编写您的代码逻辑
    # cat("heat_file:", heat_file, "\n")
    # cat("tree_file:", tree_file, "\n")
    # cat("cluster_file:", cluster_file, "\n")
    # cat("gff_file:", gff_file, "\n")
    # cat("sample_color_file:", sample_color_file, "\n")
    # cat("region:", region, "\n")
    # cat("output_prefix:", output_prefix, "\n")

    ##### 处理heat_file #####
    ## 定义heat file里面的value 与 color的对应关系
    ## 9：#FF0000，0：#FFFF00，-9：#0000FF，NA：#808080
    heat_color_mapping <- list(
        "9" = "#FF0000",
        "0" = "#FFFF00",
        "-9" = "#0000FF",
        "NA" = "#808080"
    )


    heat_data = process_heat_file(heat_file)
    # head(heat_data,10)

    ##### 处理tree_file #####
    tree_string <- readLines(tree_file, n = 1)
    tree_string <- sub(";+$", "", tree_string)
    tree_seq <- substr(tree_string, 2, nchar(tree_string)-1)
    # cat("tree_string:", tree_string, "\n")


    ##### 处理 region #####
    interval_size = get_interval_size(region)

    ## 根据interval_size计算出各个变量的值
    tree_width <<- interval_size * (1/7)
    tree_heat_gap <<- interval_size * (1/7) *(1/20)
    sample_width <<- interval_size * (1/7) *(1/10)
    heat_cluster_gap <<- interval_size * (1/7) *(1/20)
    cluster_width <<- interval_size * (1/7) *(1/20)
    cluster_text_gap <<- interval_size * (1/7) *(1/200)
    text_width <<- interval_size * (1/7) *(1/5)
    gene_height <<- nrow(heat_data) *(1/21)



    ## x轴最大值 等于树占用的1000 + 样本占用的100 + 基因区间大小 + cluster颜色50
    max_x = tree_width + sample_width + tree_heat_gap + interval_size + heat_cluster_gap + cluster_width+ cluster_text_gap + text_width 
    ## y轴最大值 等于 样本占用的1047 + 50基因结构
    white_space <<- nrow(heat_data) *(1/21)
    max_y = nrow(heat_data) + white_space
    # print(max_y)
    

    ##### 处理 cluster_file #####
    cluster_data = read.table(cluster_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    names(cluster_data) = c("sample", "cluster")


    ##### 处理 sample_color_file #####
    ## order   exp_id  vcf_id  group   color   pch
    sample_color_data = read.table(sample_color_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,comment.char = '')
    names(sample_color_data) = c("order", "exp_id", "vcf_id", "group", "color", "pch")
    sample_color_data2 = sample_color_data[,c("vcf_id", "color")]
    names(sample_color_data2) = c("sample", "color") 
    sample_color_dict <- setNames(sample_color_data2$color, sample_color_data2$sample)


    ##### 处理 gff_file #####
    combined_data = process_gff_file3(gff_file, region)
    # print(combined_data)
    # head(combined_data,40)

    

    ##### 准备画图 #####
    pdf(file=paste(output_prefix,"pdf",sep="."),height = paper_size_height,width = paper_size_width)
    par(mar=c(1, 1, 1, 1))
    plot(1,1,xlim=c(0,max_x),ylim=c(-max_y, 0),type="n",axes=F,xlab="",ylab="")
    # -----------> x = c(0,max_x)
    # |
    # |
    # ↓
    # y = c(0,-max_y)

    ##### 画树 and 样本 #####
    sample_list_max_depth <- post_traverse_tree(tree_seq, 0)
    sample_list <- sample_list_max_depth[[1]]
    max_depth <- sample_list_max_depth[[2]]
    cell_unit = (max_depth*1.0)/(tree_width)
    post_traverse_tree_and_plot(tree_seq, list(), 0, 0, cell_unit,sample_color_dict)

    ##### 画heat matrix(depth) #####
    draw_heat_matrix2(heat_data, interval_size, sample_list)


    ##### 画cluster #####
    draw_cluster_rectangles2(cluster_data, sample_list,interval_size)

    if (interval_size <= 20000){
      ##### 画基因结构 #####
      draw_gene_structure3(combined_data)
      ##### mark SNP/INDEL #####
      ## 用虚线画出SNP/INDEL的位置，与heat matrix的位置对应起来
      # draw_snp_indel_markers(heat_data, combined_data, interval_size)
    }else{
      draw_gene_structure4(combined_data)
      # draw_snp_indel_markers2(heat_data, combined_data, interval_size)
    }

    ##### 画legend #####
    draw_sample_group_legend(sample_color_data, tree_width * (1/8), -gene_height, max_x ,max_y)

    dev.off()
}

draw_sample_group_legend <- function(sample_color_data, start_x, start_y,max_x,max_y) {
  # Get unique groups and their colors
  group_colors <- unique(sample_color_data[, c("group", "color")])

  # Set initial x and y coordinates
  x <- start_x
  y <- start_y

  unit_x = paper_size_width/max_x
  unit_y = paper_size_height/max_y

  legend_rect_height = gene_height*(1/10)
  paper_legend_size_height = legend_rect_height * unit_y
  width_count = paper_legend_size_height / unit_x

  # Loop through each unique group
  for (i in 1:nrow(group_colors)) {
    # Draw a rectangle with the group's color
    rect(x, y, x + width_count, y - legend_rect_height, col = group_colors$color[i], border = "black")

    # Write the group's name next to the rectangle,左对齐
    text(x + width_count + cluster_text_gap, y - (legend_rect_height/2), group_colors$group[i], cex = 0.5, pos=4)

    # Update the y coordinate for the next group
    y <- y - 2*legend_rect_height
  }
}

draw_gene_structure3 <- function(combined_data) {
  # Prepare data
  gene_structure_df <- combined_data
  names(gene_structure_df)[4] <- "feature"
#   print(combined_data)

  # Set parameters
  height_UTR <- gene_height/10
  height_intron <- gene_height/25
  height_CDS <- gene_height/10
#   height_other <- 1
  base_y <- -gene_height / 2
  arrow_y <- base_y
  
  min_pos <- min(gene_structure_df$start) 
  base_x = tree_width + sample_width + tree_heat_gap
  # Draw gene structure
  for (i in seq(1, nrow(gene_structure_df), 1)) {
    arrow_y <- base_y
    plotx <- c()
    ploty <- c()
    cl <- NA
    if (gene_structure_df$feature[i] == "CDS") {
      ploty <- c(base_y + height_CDS, base_y + height_CDS, base_y - height_CDS, base_y - height_CDS)
      cl <- "#007BA7"
      arrow_y <- base_y - height_CDS/2
    } else if (gene_structure_df$feature[i] == "UTR") {
      ploty <- c(base_y + height_UTR, base_y + height_UTR, base_y - height_UTR, base_y - height_UTR)
      arrow_y <- base_y - height_UTR/2
    } else if (gene_structure_df$feature[i] == "intron") {
      ploty <- c(base_y + height_intron, base_y + height_intron, base_y - height_intron, base_y - height_intron)
      cl <- "black"
      arrow_y <- base_y - height_intron/2
    } else if (gene_structure_df$feature[i] == "intergenic") {
      ploty <- c(base_y, base_y)
    } else { ## gene name
      gene_name <- gene_structure_df$feature[i]
      gene_start <- gene_structure_df$start[gene_structure_df$feature == gene_name]
      gene_end <- gene_structure_df$end[gene_structure_df$feature == gene_name]
      gene_center <- (gene_start + gene_end) / 2 - min_pos
      gene_name = paste(gene_name,"(",gene_structure_df$strand[i],")",sep="")
      text(base_x+gene_center, base_y + 30, gene_name, cex = 0.7)
      next
    }

    
    plotx <- c(base_x+gene_structure_df$start[i] - min_pos, base_x+gene_structure_df$end[i] - min_pos, base_x+gene_structure_df$end[i] - min_pos, base_x+gene_structure_df$start[i] - min_pos)
    
    if (gene_structure_df$feature[i] == "intergenic") {
      # lines(plotx[1:2], ploty[1:2], col = "black", lwd = 1)
    } else if (gene_structure_df$feature[i] == "intron"){
        ## 如果是intron 画一个曲线
        intro_len = gene_structure_df$end[i] - gene_structure_df$start[i]
        lines(c(base_x+gene_structure_df$start[i] - min_pos, base_x+gene_structure_df$start[i] - min_pos + intro_len/2), c(base_y,base_y-height_intron*(1)), col = "black", lwd = 1)
        lines(c(base_x+gene_structure_df$start[i] - min_pos + intro_len/2, base_x+gene_structure_df$end[i] - min_pos), c(base_y-height_intron*(1),base_y), col = "black", lwd = 1)

    } else {
      polygon(plotx, ploty, col = cl, lwd = 0.5)
      arrow_x <- (gene_structure_df$start[i] + gene_structure_df$end[i]) / 2 - min_pos
      if (gene_structure_df$strand[i] == "+") {
        arrow_direction <- ">"
      } else {
        arrow_direction <- "<"
      }
      # 用红色标识方向
    #   text(base_x + arrow_x, arrow_y, arrow_direction, cex = 0.7, col = "red")
    }
  }

  # Draw axis
  region_length <- max(gene_structure_df$end) - min_pos + 1
  lines(c(base_x+0, base_x+region_length), c(base_y + 10, base_y + 10), col = "black", lwd = 1)
  for (i in seq(0, region_length, 1000)) {
    lines(c(base_x+i, base_x+i), c(base_y + 10, base_y + 11), lwd = 1)
    x = formatC(i + min_pos, format = "f", big.mark = ",", digits = 0)
    # print(x)
    text(base_x+i, base_y + 15, x, cex = 0.5)
  }

}

draw_snp_indel_markers <- function(heat_data, combined_data, interval_size) {
  base_y <- -gene_height / 2
  min_pos <- min(combined_data$start)
  max_pos <- max(combined_data$end)
  base_x <- tree_width + sample_width + tree_heat_gap
  height_CDS <- 5
  rect_width <- interval_size/ncol(heat_data)
  
  colnames_heat_data <- colnames(heat_data)
  
  for (colname in colnames_heat_data) {
    # 提取位置信息
    split_colname <- strsplit(colname, "_")[[1]]
    position <- as.integer(split_colname[2])
    
    # 计算 SNP/INDEL 在基因结构上的 x 坐标
    snp_indel_x <- base_x + position - min_pos
    
    # 在基因结构上用虚线标识 SNP/INDEL 位置
    lines(c(snp_indel_x, snp_indel_x), c(base_y - height_CDS-1, base_y + height_CDS+1), col = "gray", lwd = 1, lty = 2)
    
    # 计算 SNP/INDEL 在 heat matrix 的列位置
    col_index <- which(colnames_heat_data == colname)
    
    # 用虚线连接 SNP/INDEL 在 heat matrix 的列位置与基因结构上的虚线
    for (row_index in 1:nrow(heat_data)) {
      heat_matrix_x <- tree_width + sample_width + tree_heat_gap + (col_index - 1) * rect_width+rect_width/2
      heat_matrix_y <- -gene_height
      lines(c(heat_matrix_x, snp_indel_x), c(heat_matrix_y, base_y-height_CDS-1), col = "gray", lwd = 0.5, lty = 1)
    }
  }
}


draw_cluster_rectangles2 <- function(cluster_data, sample_list, interval_size) {
  cluster_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#FFFF00")
  
  cluster_data_sample_order <- cluster_data[match(sample_list, cluster_data$sample), ]
  n_samples <- nrow(cluster_data_sample_order)
  
  rect_width <- cluster_width
  rect_height <- 1

  current_cluster <- cluster_data_sample_order$cluster[1]
  n <- 1
  cluster_color <- cluster_colors[(n %% length(cluster_colors)) + 1]
  cluster_start_index <- 1

  for (i in 1:n_samples) {
    sample_cluster <- cluster_data_sample_order$cluster[i]
    sample_name <- cluster_data_sample_order$sample[i]
    
    if (current_cluster != sample_cluster || i == n_samples) {
      if (i == n_samples) {
        i <- i + 1
      }
      cluster_middle_index <- cluster_start_index + floor((i - cluster_start_index) / 2)
      text_x <- tree_width + sample_width + tree_heat_gap + interval_size + heat_cluster_gap + cluster_width+ cluster_text_gap
      text_y <- -cluster_middle_index + rect_height / 2 - gene_height
      text(text_x, text_y, paste0("hap ", current_cluster), cex = 0.6, pos = 4)

      if (i != n_samples) {
        n <- n + 1
        current_cluster <- sample_cluster
        cluster_color <- cluster_colors[(n %% length(cluster_colors)) + 1]
        cluster_start_index <- i
      }
    }

    x_left <- tree_width + sample_width + tree_heat_gap + interval_size + heat_cluster_gap
    x_right <- x_left + rect_width
    y_top <- -i + rect_height / 2 - gene_height
    y_bottom <- y_top - rect_height
    
    rect(x_left, y_bottom, x_right, y_top, col = cluster_color, border = NA)
  }
}



depth2color <- function(depth) {
  color_list <- c("#BEBEBE", "#C4C4AA", "#CBCB96", "#D2D282", "#D9D96E", "#E0E05A", "#E7E746", "#EDED32", "#F4F41E", "#FBFB0A", "#FFF100",
                  "#FFD600", "#FFBB00", "#FFA100", "#FF8600", "#FF6B00", "#FF5000", "#FF3500", "#FF1A00", "#FF0000")
  if (depth >= 2) {
    return(color_list[length(color_list)])
  } else {
    step <- 2.0 / length(color_list)
    return(color_list[floor(depth / step) + 1])
  }
}

draw_heat_matrix2 <- function(heat_data, interval_size, sample_list) {
  sorted_heat_data <- heat_data[match(sample_list, rownames(heat_data)), ]
  n_rows <- nrow(sorted_heat_data)
  n_cols <- ncol(sorted_heat_data)

  # 计算每个矩形的宽度和高度
  rect_width <- interval_size / n_cols
  rect_height <- 1

  # 遍历热图数据并绘制矩阵
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # 获取单元格的深度值并转换为颜色
      cell_depth <- sorted_heat_data[i, j]
      cell_color <- depth2color(cell_depth)

      # 计算矩形的边界
      x_left <- tree_width + sample_width + tree_heat_gap + (j - 1) * rect_width
      x_right <- x_left + rect_width
      y_top <- -1 * i + rect_height / 2 - gene_height
      y_bottom <- y_top - rect_height

      # 绘制无边框矩形
      rect(x_left, y_bottom, x_right, y_top, col = cell_color, border = NA)
    }
  }
}



draw_heat_matrix <- function(heat_data, interval_size, heat_color_mapping,sample_list) {
  sorted_heat_data <- heat_data[match(sample_list, rownames(heat_data)),]
  n_rows <- nrow(sorted_heat_data)
  n_cols <- ncol(sorted_heat_data)

  # 计算每个矩形的宽度和高度
  rect_width <- interval_size / n_cols
  rect_height <- 1

  # 遍历热图数据并绘制矩阵
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # 获取单元格的值和颜色
      cell_value <- as.character(sorted_heat_data[i, j])
      cell_color <- heat_color_mapping[[cell_value]]

      # 计算矩形的边界
      x_left <- tree_width + sample_width + tree_heat_gap +(j - 1) * rect_width
      x_right <- x_left + rect_width
      y_top <- -1 * i + rect_height/2- gene_height
      y_bottom <- y_top - rect_height
      

      # 绘制无边框矩形
      rect(x_left, y_bottom, x_right, y_top, col = cell_color, border = NA)
    }
  }
}


draw_line <- function(x1, y1, x2, y2, lty = 1) {
  if (x1 == x2 && y1 == y2) {
    return()
  }
  y1 = y1-gene_height
  y2 = y2-gene_height
  lines(x = c(x1, x2), y = c(y1, y2), lty = lty)
}

draw_dashed_line <- function(x1, y1, x2, y2) {
  lines(x = c(x1, x2), y = c(y1, y2), lty = 2)
}

draw_yellow_line <- function(x1, y1, x2, y2, lty = 1) {
  lines(x = c(x1, x2), y = c(y1, y2), col = "yellow", lty = lty)
}



substr_count <- function(string, pattern) {
  return(sum(str_count(string, fixed(pattern))))
}


find_comma_index <- function(tree_seq) {
  comma_index_list <- c()
  for (index in 1:nchar(tree_seq)) {
    if (substr(tree_seq, index, index) == ",") {
      comma_index_list <- c(comma_index_list, index)
    }
  }
  
  i_list <- c()
  for (i in comma_index_list) {

    if (substr_count(substr(tree_seq, 1, i-1), "(") == substr_count(substr(tree_seq, 1, i-1), ")") &&
        substr_count(substr(tree_seq, i + 1, nchar(tree_seq)), ")") == substr_count(substr(tree_seq, i + 1, nchar(tree_seq)), "(")) {
      i_list <- c(i_list, i)
    }
  }
  
  return(i_list)
}

post_traverse_tree_and_plot <- function(tree_seq, samples_up_to_you, father_depth, cur_depth, cell_unit,sample_color_dict) {
  comma_index <- find_comma_index(tree_seq)
  comma_index_len <- length(comma_index)
  
  if (comma_index_len == 0) {
    sample_name <- strsplit(tree_seq, ":")[[1]][1]
    samples_up_to_you <- append(samples_up_to_you, sample_name)
    row_num <- length(samples_up_to_you)
    father_col <- as.integer(father_depth / cell_unit)
    cur_col <- as.integer(cur_depth / cell_unit)
    
    draw_line(father_col , -row_num, cur_col , -row_num)
    # print(paste0("father_col: ", father_col, " cur_col: ", cur_col, " row_num: ", row_num))
    # print(samples_up_to_you)

    sample_color <- paste0("", sample_color_dict[[sample_name]])  # 获取样本的颜色并添加 "#" 号
    # print(paste0("sample_name: ", sample_name, " sample_color: ", sample_color))
  
    # Define the rectangle's position and dimensions
    x_left <- cur_col
    x_right <- x_left + sample_width
    y_bottom <- -row_num-0.5 - gene_height
    y_top <- -row_num+0.5 - gene_height

    # Draw the rectangle without a border
    rect(x_left, y_bottom, x_right, y_top, col = sample_color, border = NA)


    return(list(samples_up_to_you, row_num))
  }
  
  seq_list <- list()
  pre_start <- 1
  for (i in 1:comma_index_len) {
    indx <- comma_index[i]
    sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, indx - 1))
    seq_list[[i]] <- sub_tree_seq_branch_len
    pre_start <- indx + 1
    if (i == comma_index_len) {
      sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, nchar(tree_seq)))
      seq_list[[i + 1]] <- sub_tree_seq_branch_len
    }
  }
  
  
  row_num_list <- list()
  for (i in 1:length(seq_list)) {
    sub_tree_seq <- seq_list[[i]][[1]]
    sub_branch_len <- seq_list[[i]][[2]]
    sub_sample_list_row_num <- post_traverse_tree_and_plot(sub_tree_seq, samples_up_to_you, cur_depth, cur_depth + sub_branch_len, cell_unit,sample_color_dict)
    samples_up_to_you <- sub_sample_list_row_num[[1]]
    row_num <- sub_sample_list_row_num[[2]]
    row_num_list[[i]] <- row_num
  }
  
  father_col <- as.integer(father_depth / cell_unit)
  cur_col <- as.integer(cur_depth / cell_unit)
#   col_num_list <- (father_col + 1):(cur_col + 1)
  row_num <- (row_num_list[[1]] + row_num_list[[length(row_num_list)]]) / 2

  # plot node
  #    |
  #----|
  #    |

  draw_line(cur_col, -row_num_list[[1]], cur_col, -row_num_list[[length(row_num_list)]])
  draw_line(father_col, -row_num, cur_col, -row_num)
  

  #plot root
  # |
  # |
  # |
    #   if (father_depth == 0 && cur_depth == 0) {
    #     draw_line(1, -row_num_list[[1]], 1, -row_num_list[[length(row_num_list)]])
    #   }
  return(list(samples_up_to_you, row_num))
}

rm_bracket2 <- function(tree_seq) {
  if (substr(tree_seq, 1, 1) == "(") {
    pos <- str_locate_all(tree_seq, fixed(")"))[[1]]
    # 取出最后一个右括号的位置
    rbound <- tail(pos, 1)[1]
    if (rbound == nchar(tree_seq)) {
      return(list(substr(tree_seq, 2, nchar(tree_seq) - 1), 1))
    } else {
      return(list(substr(tree_seq, 2, rbound-1), as.numeric( strsplit(tree_seq, ":")[[1]][length(strsplit(tree_seq, ":")[[1]])] )))
    }
  } else {
    if (grepl(":", tree_seq)) {
      return(list(tree_seq, as.numeric(strsplit(tree_seq, ":")[[1]][length(strsplit(tree_seq, ":")[[1]])])))
    } else {
      return(list(tree_seq, 1))
    }
  }
}

post_traverse_tree <- function(tree_seq, cur_branch_len) {
  comma_index <- find_comma_index(tree_seq)
  comma_index_len <- length(comma_index)

  if (comma_index_len == 0) {
    sample_name <- strsplit(tree_seq, ":")[[1]][1]
    return(list(list(sample_name), cur_branch_len))
  }

  seq_list <- list()
  pre_start <- 1
  for (i in 1:comma_index_len) {
    indx <- comma_index[i]
    sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, indx - 1))
    seq_list[[i]] <- sub_tree_seq_branch_len
    pre_start <- indx + 1
    if (i == comma_index_len) {
      sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, nchar(tree_seq)))
      seq_list[[i + 1]] <- sub_tree_seq_branch_len
    }
  }

  sample_list <- list()
  cur_max_branch_len <- 0
  for (i in 1:length(seq_list)) {
    sub_tree_seq <- seq_list[[i]][[1]]
    sub_branch_len <- seq_list[[i]][[2]]
    sub_sample_list_max_depth <- post_traverse_tree(sub_tree_seq, cur_branch_len + sub_branch_len)
    sub_sample_list <- sub_sample_list_max_depth[[1]]
    sub_max_depth <- sub_sample_list_max_depth[[2]]
    sample_list <- c(sample_list, sub_sample_list)
    cur_max_branch_len <- max(cur_max_branch_len, sub_max_depth)
  }

  return(list(sample_list, cur_max_branch_len))
}



process_gff_file3 <- function(gff_file, region){
  gff_data <- read.table(gff_file, sep="\t", header = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "#",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

  # 筛选 region 对应的染色体数据
  # 潜在bug,提供的区间 和 gff 文件中的染色体名称不一致
  region_seqid <- unlist(strsplit(region, ":"))[1]
  gff_data_filtered <- gff_data[gff_data$seqid == region_seqid,]

  # 筛选 region 对应的区间数据
  interval_string <- unlist(strsplit(region, ":"))[2]
  interval_bounds <- as.numeric(unlist(strsplit(interval_string, "-")))
  gff_data_filtered <- gff_data_filtered[gff_data_filtered$start >= interval_bounds[1] & gff_data_filtered$end <= interval_bounds[2],]

  if (nrow(gff_data_filtered) > 0) {
    # 提取 gene、CDS、UTR 和 exon 数据
    gene_data <- gff_data_filtered[gff_data_filtered$type == "gene", c("start", "end", "attributes", "strand")]
    cds_data <- gff_data_filtered[gff_data_filtered$type == "CDS", c("start", "end", "strand", "attributes")]
    utr_data <- gff_data_filtered[gff_data_filtered$type %in% c("five_prime_UTR", "three_prime_UTR"), c("start", "end", "strand", "attributes")]
    exon_data <- gff_data_filtered[gff_data_filtered$type == "exon", c("start", "end", "strand", "attributes")]

    # 获取第一个转录本的ID
    # 潜在bug: 并不是所有物种都是用这种正则能匹配到 gene id
     first_transcript_ids <- sapply(strsplit(gene_data$attributes, ";"), function(x) {
      id_field <- x[grep("^ID=", x)]
      if (length(id_field) > 0) {
        sub("^ID=", "", id_field) %>% paste0(., ".1")
      } else {
        NA
      }
    })

    print(first_transcript_ids)
    # gene_data <- subset(gene_data, select = -attributes)
    cds_data_filtered <- data.frame()
    utr_data_filtered <- data.frame()
    exon_data_filtered <- data.frame()
    intron_data <- data.frame()

    for (transcript_id in first_transcript_ids) {
      cds_data_filtered <- rbind(cds_data_filtered,
                                gff_data_filtered[gff_data_filtered$type == "CDS" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      utr_data_filtered <- rbind(utr_data_filtered,
                                gff_data_filtered[gff_data_filtered$type %in% c("five_prime_UTR", "three_prime_UTR") & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      exon_data_filtered <- rbind(exon_data_filtered,
                                gff_data_filtered[gff_data_filtered$type == "exon" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      # 计算 intron 数据
      tmp_exon_data = gff_data_filtered[gff_data_filtered$type == "exon" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")]
      if (nrow(tmp_exon_data) > 1) {
        for (i in 1:(nrow(tmp_exon_data) - 1)) {
          if (tmp_exon_data[i, "end"] + 1 < tmp_exon_data[i + 1, "start"]) {
            intron_data <- rbind(intron_data, data.frame(start = tmp_exon_data[i, "end"] + 1, end = tmp_exon_data[i + 1, "start"] - 1))
          }
        }
      }
    }
    
    # 用筛选后的数据更新CDS、UTR和exon数据
    cds_data <- cds_data_filtered
    utr_data <- utr_data_filtered
    exon_data <- exon_data_filtered
    # # 筛选第一个转录本的CDS、UTR和exon数据
    # if (nrow(cds_data) > 0) {
    #   cds_data <- cds_data[grep(paste0("Parent=", first_transcript_id), cds_data$attributes), ]
    # }
    # if (nrow(utr_data) > 0) {
    #   utr_data <- utr_data[grep(paste0("Parent=", first_transcript_id), utr_data$attributes), ]
    # }
    # if (nrow(exon_data) > 0) {
    #   exon_data <- exon_data[grep(paste0("Parent=", first_transcript_id), exon_data$attributes), ]
    # }

    # 删除不再需要的attributes列
    # cds_data <- subset(cds_data, select = -attributes)
    # utr_data <- subset(utr_data, select = -attributes)
    # exon_data <- subset(exon_data, select = -attributes)

    # 从 attributes 列中提取 ID 并分配给 gene_data 的 type 列
    if (nrow(gene_data) > 0) {
        gene_data$type <- sapply(strsplit(gene_data$attributes, ";"), function(x) {
        id_field <- x[grep("^ID=", x)]
        ## gene id 加上.1 
        id_field <- paste0(id_field, ".1")
        if (length(id_field) > 0) {
          sub("^ID=", "", id_field)
        } else {
          NA
        }
      })
      gene_data <- subset(gene_data, select = -attributes)
    }

    

    # 为所有数据添加 type 列，前提是数据框不为空
    if (nrow(cds_data) > 0) cds_data$type <- "CDS"
    if (nrow(utr_data) > 0) utr_data$type <- "UTR"
    if (nrow(exon_data) > 0) exon_data$type <- "exon"
    if (nrow(intron_data) > 0) intron_data$type <- "intron"

    ## 给intron_data添加strand列，根据其他data的strand列来判断
    if (nrow(intron_data) > 0) {
      intron_data$strand <- exon_data$strand[1]
      # 并且将列排一下序，让其和其他data的列顺序一致
      intron_data <- intron_data[,c("start", "end", "strand", "type")]
    }

    # 合并数据
    combined_data <- rbind(gene_data, cds_data, utr_data, intron_data)

    # 添加 intergenic 区域
    # 潜在bug： gene区间和当前转录本占用的区间不一定一致，因为gene区间是所有转录本的区间的并集
    # 优化参考：https://en.define.sh/How-to-Calculate-Intersection-Union-And-Difference-of-Intervals
    intergenic_start <- interval_bounds[1]
    intergenic_end <- interval_bounds[2]
    intergenic_data <- data.frame()
    if (nrow(gene_data) > 0) {
      for (i in 1:nrow(gene_data)) {
        if (gene_data[i, "start"] > intergenic_start) {
          intergenic_data <- rbind(intergenic_data, data.frame(start = intergenic_start, end = gene_data[i, "start"] - 1, type = "intergenic"))
        }
        intergenic_start <- gene_data[i, "end"] + 1
      }
    }
    if (intergenic_start <= intergenic_end) {
      intergenic_data <- rbind(intergenic_data, data.frame(start = intergenic_start, end = intergenic_end, type = "intergenic"))
    }

    # 合并 intergenic_data 并按照 start 列排序
    if (nrow(intergenic_data) > 0) {
        intergenic_data$strand <- "x"
        combined_data <- rbind(combined_data, intergenic_data)
    }
    combined_data <- combined_data[order(combined_data$start), ]
  } else {
    combined_data <- data.frame(start = interval_bounds[1], end = interval_bounds[2], strand <- "x", type = "intergenic")
  }

  return(combined_data)
}


get_interval_size <- function(region){
    # 分割字符串并提取区间部分
    interval_string <- unlist(strsplit(region, ":"))[2]
    # 分割区间字符串并转换为数值
    interval_bounds <- as.numeric(unlist(strsplit(interval_string, "-")))
    # 计算区间大小（包含边界）
    interval_size <- interval_bounds[2] - interval_bounds[1] + 1
    # cat("interval_size:", interval_size, "\n")

    return(interval_size)
}


process_heat_file <- function(heat_file) {
  
  # 读取 heat_file 并将所有元素转换为字符串
  heat_data <- read.table(heat_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE,check.names = FALSE)
  # 删除前两列
  # heat_data <- heat_data[, -c(1:2)]
  # # 将所有数据框中的元素转换为字符串
  # heat_data[] <- lapply(heat_data, as.character)
  # # 转置数据框
  # heat_data <- t(heat_data)
  # # NA 替换为 "NA"
  # heat_data[is.na(heat_data)] <- "NA"
  
  return(heat_data)
}

draw_gene <- function(start, end, strand, base_y, height,len) {
  mid =  (start+end)/2
  start = mid - len
  end = mid + len
  if (strand == "+") {  # 如果链是"+"
    # 画一个向右的三角形
    polygon(c(start, end, start), c(base_y+height, base_y, base_y-height), col="#007BA7",border = NA)
  } else {  # 如果链是"-"
    # 画一个向左的三角形
    polygon(c(end, start, end), c(base_y-height, base_y, base_y+height), col="#cf5843", border = NA)
  }
}

draw_gene_structure4 <- function(combined_data) {
  # Prepare data
  gene_structure_df <- combined_data
  names(gene_structure_df)[4] <- "feature"
  print(head(gene_structure_df))

  # Set parameters
  height_UTR <- 5
  height_intron <- 2
  height_CDS <- 5
#   height_other <- 1
  base_y <- -gene_height / 2
  arrow_y <- base_y
  
  min_pos <- min(gene_structure_df$start) 
  base_x = tree_width + sample_width + tree_heat_gap
  
  region_start = min(gene_structure_df$start)
  region_end = max(gene_structure_df$end)
  ## 如果区间小于2M就画基因结构，如果大于2M 所有基因用▶或者◀代替
  
  #draw gene as ▶ or ◀
  lines(c(base_x+region_start - min_pos , base_x+region_end - min_pos), c(base_y,base_y), col = "black", lwd = 1)
  for (i in seq(1, nrow(gene_structure_df), 1)) {
      if (gene_structure_df$feature[i] == "CDS" || gene_structure_df$feature[i] == "UTR" || 
          gene_structure_df$feature[i] == "intron" || gene_structure_df$feature[i] == "intergenic")
      {
          next
      }else { ## gene name
          arrow_x <- (gene_structure_df$start[i] + gene_structure_df$end[i]) / 2 - min_pos
          gene_len <- (gene_structure_df$end[i] - gene_structure_df$start[i] ) / 2 
          # print(gene_structure_df$feature[i])
          ## 画个箭头
          draw_gene(base_x + gene_structure_df$start[i] - min_pos, base_x + gene_structure_df$end[i] - min_pos,
                    gene_structure_df$strand[i], base_y, height_CDS,gene_len)
          print("111111111")
          start_a = base_x + gene_structure_df$start[i] - min_pos
          end_a = base_x + gene_structure_df$end[i] - min_pos
          mid_a = (start_a + end_a)/2
          # print(mid_a, start_a,end_a)
          # arrows(x0 = mid_a, y0 = base_y-10, x1 = mid_a, y1 = base_y-5, col="black")
          
      }
  }

  

  # Draw axis
  region_length <- max(gene_structure_df$end) - min_pos + 1
  lines(c(base_x+0, base_x+region_length), c(base_y + 10, base_y + 10), col = "black", lwd = 1)
  
  ## 确定一个合适的间隔，stepsize，这个值是1000的倍数，保证这个区间最多只有10个刻度
  stepsize = 1000
  while (region_length/stepsize > 10){
      stepsize = stepsize + 1000
  }


  for (i in seq(0, region_length, stepsize)) {
    lines(c(base_x+i, base_x+i), c(base_y + 10, base_y + 11), lwd = 1)
    x = formatC(i + min_pos, format = "f", big.mark = ",", digits = 0)
    # print(x)
    text(base_x+i, base_y + 15, x, cex = 0.5)
  }

}


# 处理命令行参数
process_command_line <- function() {
    # 定义命令行选项
    option_list <- list(
        make_option(c("-H", "--heat_file"), type = "character", help = "Path to the heat_file."),
        make_option(c("-t", "--tree_file"), type = "character", help = "Path to the tree_file."),
        make_option(c("-c", "--cluster_file"), type = "character", help = "Path to the cluster_file."),
        make_option(c("-g", "--gff_file"), type = "character", help = "Path to the gff_file."),
        make_option(c("-s", "--sample_color_file"), type = "character", help = "Path to the sample_color_file."),
        make_option(c("-r", "--region"), type = "character", help = "Region."),
        make_option(c("-o", "--output_prefix"), type = "character", help = "Output file prefix.")
    )

    # 解析命令行参数
    args <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)

    if (!is.null(args$help) && args$help) {
        print_help(OptionParser(option_list = option_list))
        quit(save = "no", status = 0)
    }

    # 检查参数是否完整
    required_args <- c("heat_file", "tree_file", "cluster_file", "gff_file", "sample_color_file", "region", "output_prefix")
    for (arg_name in required_args) {
        if (is.null(args$options[[arg_name]])) {
            cat("Error: Argument", arg_name, "is missing.\n")
            print_help(OptionParser(option_list = option_list))
            quit(save = "no", status = 1)
        }
        if ( !(arg_name %in% c("region", "output_prefix")) && !file.exists(args$options[[arg_name]])) {
            cat("Error: File", args$options[[arg_name]], "for argument", arg_name, "does not exist.\n")
            print_help(OptionParser(option_list = option_list))
            quit(save = "no", status = 1)
        }
    }

    main_function(args$options$heat_file, args$options$tree_file, args$options$cluster_file, args$options$gff_file,
                  args$options$sample_color_file, args$options$region, args$options$output_prefix)
}

# 如果脚本是通过 Rscript 直接运行的，处理命令行参数并调用主函数
if (interactive() == FALSE) {
    process_command_line()
}
