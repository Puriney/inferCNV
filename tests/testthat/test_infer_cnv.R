# Global data

make_fake_infercnv_obj <- function(some_matrix) {

    num_cells = ncol(some_matrix)
    num_genes = nrow(some_matrix)

    if (num_cells < 2) {
        stop("Error, need at least 2 cells in the matrix")
    }
    
    gene_order <- data.frame(chr=rep("chr1", num_genes),
                             start=1:num_genes,
                             stop=1:num_genes)
    
    midpt_cells = floor(num_cells/2)
    
    normal_cells = 1:midpt_cells
    tumor_cells = (midpt_cells+1):num_cells

    infercnv_obj <- new(
        Class = "infercnv",
        expr.data = some_matrix, 
        count.data = some_matrix,
        gene_order = gene_order,
        reference_grouped_cell_indices = list(normal=normal_cells),
        observation_grouped_cell_indices = list(tumor=tumor_cells) )

    return(infercnv_obj)

}


matrix_zeros <- matrix(rep(0,5), ncol=1)

matrix_one <- matrix(1:5, ncol=1)

matrix_one_long <- matrix(1:20, ncol=1)

matrix_one_long_2 <- matrix(c(1,2,4,7,9,11,12,14,17,19,16,14,
                              13,11,10,7,6,4,3,1), ncol=1)
matrix_two <- matrix(1:10, ncol=2)

matrix_two_long <- matrix(1:40, ncol=2)

matrix_two_long_2 <- matrix(c(1,2,4,7,9,11,12,14,17,19,16,14,13,11,10,7,6,4,3,
                              1,1,2,4,7,9,11,12,14,17,19,16,14,13,11,10,7,6,4,
                              3,1), ncol=2)

matrix_three <- matrix(1:15, ncol=3)

matrix_five <- matrix(1:25, ncol=5)

context("Test subtract_ref")

matrix_averef_five <- matrix(c(c(-101, -100, -100, -100, -99),
                               c(-101, -100, -99, -98, -99),
                               c(1, 1, 2, 3, 0),
                               c(110, 103, 90, 80, 70),
                               c(0, 0, 0, 0, 0),
                               c(100, 102, 100, 102, 102),
                               c(0, -1, -4, -1, -1),
                               c(105, 95, 80, 97, 80),
                               c(100, 99, 100, 101, 100),
                               c(0, 0, 0, 0, 0)),
                               ncol=5, byrow=FALSE)

avref_answer_1 <- matrix(0:4, ncol=1)
avref_answer_2 <- matrix(c(0:4,0:4), ncol=2)
avref_answer_3 <- matrix(c(-1:3, -1:3, -1:3), ncol=3)
avref_answer_4 <- matrix(rep(-3:1 + .5, 5),ncol=5)
avref_answer_5 <- matrix_zeros
matrix_averef_five_answer <- matrix(c(c(-1,0,0,0,0,-1,0,0,1,0),
                                      c(0,0,0,0,-1,40,33,20,10,0),
                                      c(0,0,0,0,0,0,0,0,0,0),
                                      c(0,0,-3,0,0,25,15,0,17,0),
                                      c(1,0,1,2,1,0,0,0,0,0)),
                                    ncol=10,
                                    byrow=TRUE)

test_that("subtract_ref works with one observation, one reference",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_one),
                              ref_groups=list(c(1))),
                 t(avref_answer_1))
          })
test_that("subtract_ref works with two observations, one reference",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_two),
                                           ref_groups=list(c(1))),
                 t(avref_answer_2))
          })

test_that("subtract_ref updated works with 3 observaions, two reference",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_three),
                                           ref_groups=list(c(1,3))),
                 t(avref_answer_3))
    })
test_that("subtract_ref works with 5 observations, two reference",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_five),
                                   ref_groups=list(c(2,5))),
                 t(avref_answer_4))
          })
test_that("subtract_ref works with 1 observation, 1 reference",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_zeros),
                                           ref_groups=list(c(1))),
                 t(avref_answer_5))
})

test_that("subtract_ref works with 10 obs, 5 references, 3 groups",{
    expect_equal(infercnv:::.subtract_expr(t(matrix_averef_five),
                                           ref_groups=list(c(2),c(4,6,8),c(10))),
                 matrix_averef_five_answer)
          })




context("Test center_columns")

center_sm_1 <- matrix(1:10, ncol=1)
center_sm_1_answer <- matrix(c(-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5),
                             ncol=1)
center_sm_3 <- matrix(1:21, ncol=3)
center_sm_3_answer <- matrix(rep(c(-3,-2,-1,0,1,2,3),3), ncol=3)

test_that("center_columns works with 1 observations",{
    expect_equal(infercnv:::.center_columns(expr_data=center_sm_3, method="mean"),
                 center_sm_3_answer)
    })

test_that("center_smoothed works with 3 observations",{
    expect_equal(infercnv:::.center_columns(center_sm_3, method="mean"),
                 center_sm_3_answer)
    })


context("Test below_min_mean_expr_cutoff")

below_answer_1 <- 1:5
below_answer_2 <- 1:4
below_answer_3 <- c(1)
below_answer_4 <- 1:3
below_answer_5 <- integer(0)
below_answer_6 <- 1:5

test_that(paste("below_min_mean_expr_cutoff works with one observation,",
                "cutoff too large to affect"),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_one,
                                                                        min_mean_expr=10),
                                 below_answer_1)
                })
test_that(paste("below_min_mean_expr_cutoff works with three observations,",
                "threshold too large to affect"),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_three,
                                                                        min_mean_expr=10),
                                 below_answer_2)
                })
test_that(paste("below_min_mean_expr_cutoff works with one observation,",
                "threshold excluding two."),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_one,
                                                                        min_mean_expr=2),
                                 below_answer_3)
                })
test_that(paste("below_min_mean_expr_cutoff works with three observations,",
                "threshold excluding three."),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_three,
                                                                        min_mean_expr=8.4),
                                 below_answer_4)
                })
test_that(paste("above_min_mean_expr_cutoff works with one observation,",
                "threshold including all."),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_one,
                                                            min_mean_expr=0),
                                 below_answer_5)
                })
test_that(paste("above_min_mean_expr_cutoff works with three observations,",
                "threshold excluding all."),{
                    expect_equal(infercnv:::.below_min_mean_expr_cutoff(expr_data=matrix_three,
                                                            min_mean_expr=100),
                                 below_answer_6)
                })


context("Test clear_noise")

noise_answer_1 <- matrix_one
noise_answer_2 <- matrix(c(0,0,0,4,5), ncol=1)
noise_answer_3 <- matrix_zeros
noise_answer_4 <- matrix_three
noise_answer_5 <- matrix(c(rep(0,11),12:15), ncol=3)
noise_answer_6 <- matrix(rep(0,15), ncol=3)

test_that("remove_noise works with one observation, threshold 0",{
    expect_equal(infercnv:::.clear_noise(expr_data=matrix_one,
                                         threshold=0),
                 noise_answer_1)
})
test_that(paste("remove_noise works with one observation, one ref,",
                "threshold removing some"),{
                    expect_equal(infercnv:::.clear_noise(expr_data=matrix_one,
                                                         threshold=4),
                                 noise_answer_2)
                })
test_that(paste("remove_noise works with one observation,",
                "threshold removing all"),{
                    expect_equal(infercnv:::.clear_noise(expr_data=matrix_one,
                                                         threshold=6),
                                 noise_answer_3)
                })
test_that("remove_noise works with three observation, threshold 0",{
    expect_equal(infercnv:::.clear_noise(expr_data=matrix_three,
                                         threshold=0),
                 noise_answer_4)
})
test_that("remove_noise works with three observation, threshold some",{
    expect_equal(infercnv:::.clear_noise(expr_data=matrix_three,
                                         threshold=12),
                 noise_answer_5)
})
test_that("remove_noise works with three observation, threshold all",{
    expect_equal(infercnv:::.clear_noise(expr_data=matrix_three,
                                         threshold=100),
                 noise_answer_6)
})


context("Test remove_tails")

tail_answer_1 <- c()
tail_answer_2 <- c(1:5, 16:20)
tail_answer_3 <- c(2:6, 13:17)
tail_answer_4 <- c(5:9, 11:15)
tail_answer_5 <- c(1, 5)

test_that(paste("remove tails works with one contig,",
                "one observation, no tail length"),{
    expect_equal(infercnv:::.remove_tails(smooth_matrix=matrix_one,
                                          chr=1:length(matrix_one),
                                          tail_length=0),
                 tail_answer_1)
                })
test_that("remove tails works with one contig, one observation, tail length 5",{
    expect_equal(infercnv:::.remove_tails(smooth_matrix=matrix_one_long,
                                          chr=1:length(matrix_one_long),
                                          tail_length=5),
                 tail_answer_2)
})
test_that("remove tails works with 3 contig, one observation, tail length 5",{
    expect_equal(infercnv:::.remove_tails(smooth_matrix=matrix_one_long,
                                          chr=2:17,
                                          tail_length=5),
                 tail_answer_3)
})
test_that("remove tails works with 3 contig, two observations, tail length 5",{
    expect_equal(infercnv:::.remove_tails(smooth_matrix=matrix_two_long,
                                          chr=5:15,
                                          tail_length=5),
                 tail_answer_4)
})
test_that(paste("remove tails works with one contig, one observation,",
                "tail length longer than contig"),{
                    expect_equal(infercnv:::.remove_tails(smooth_matrix=matrix_one,
                                                          chr=1:length(matrix_one),
                                                          tail_length=100),
                                 tail_answer_5)
                })


context("smooth_window")

smooth_answer_1 <- matrix_one

smooth_answer_2 <- matrix_one

# smooth_answer_3 <- matrix(c(1.00,2.53,4.60,6.60,8.60,10.60,12.60,14.60,
#                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
#                            7.60,6.00,4.20,2.73,1.00), ncol=1)
smooth_answer_3 <- matrix(c(2.88, 4.44, 6.67, 8.78, 10.67, 12.44, 14.44, 16.11, 16.78, 16, 14.44, 12.78, 11.11, 9.44, 7.56, 5.89, 4.22, 3.13, 2.17), ncol=1)


#smooth_answer_4 <- matrix(c(1.00,2.53,4.60,6.60,8.60,10.60,12.60,14.60,
#                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
#                            7.60,6.00,4.20,2.73,1.00,
#                            1.00,2.53,4.60,
#                            6.60,8.60,10.60,12.60,14.60,
#                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
#                            7.60,6.00,4.20,2.73,1.00), ncol=2)

smooth_answer_4 <- matrix(c(2.88, 4.44, 6.67, 8.78, 10.67, 12.44, 14.44, 16.11, 16.78, 16, 14.44, 12.78, 11.11, 9.44, 7.56, 5.89, 4.22, 3.13, 2.17, 2.88, 4.44, 6.67, 8.78, 10.67, 12.44, 14.44, 16.11, 16.78, 16, 14.44, 12.78, 11.11, 9.44, 7.56, 5.89, 4.22, 3.13, 2.17), ncol=2)

smooth_answer_5 <- matrix(c(1.67, 2.25, 3, 3.75, 4.33), ncol=1)

test_that("smooth_window works with one observation, window_length 0",{
    expect_equal(infercnv:::.smooth_window(data=matrix_one,
                                           window_length=0),
                 smooth_answer_1)
         })

test_that("smooth_window works with one observation, window_length 1",{
    expect_equal(infercnv:::.smooth_window(data=matrix_one,
                                           window_length=1),
                 smooth_answer_2)
})

test_that("smooth_window works with one observation, window_length 5", {isTRUE(all.equal(
                                           infercnv:::.smooth_window(data=matrix_one_long_2,
                                           window_length=5),
                                           smooth_answer_3))
})

test_that("smooth_window works with two observations, window_length 5",{isTRUE(all.equal(
                                                 infercnv:::.smooth_window(data=matrix_two_long_2,
                                                 window_length=5),
                                                 smooth_answer_4))
})

test_that(paste("smooth_window works with one observation,",
                "window_length longer than measurements"),{
                    isTRUE(all.equal(infercnv:::.smooth_window(data=matrix_one,
                                                               window_length=100),
                                                               smooth_answer_5))
})


context("create_sep_list")
create_sep_list_answer_1 <- list()
create_sep_list_answer_1[[1]] <- list()
create_sep_list_answer_1[[1]][[1]] <- c(0,0,0,0)
create_sep_list_answer_1[[2]] <- list()
create_sep_list_answer_1[[2]][[1]] <- c(0,0,0,0)

create_sep_list_answer_2 <- list()
# Column
create_sep_list_answer_2[[1]] <- list()
create_sep_list_answer_2[[1]][[1]] <- c(2,0,2,5)
create_sep_list_answer_2[[1]][[2]] <- c(8,0,8,5)
# Row
create_sep_list_answer_2[[2]] <- list()
create_sep_list_answer_2[[2]][[1]] <- c(0,4,9,4)
create_sep_list_answer_2[[2]][[2]] <- c(0,2,9,2)
 
test_that("create_sep_list works with 0 row and column and no seps",{
    expect_equal(create_sep_list(0, 0),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 0 row and column and seps",{
    expect_equal(create_sep_list(0, 0, 1:5, 3:6),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 10 row, 0 column and no seps",{
    expect_equal(create_sep_list(10, 0),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 0 row, 10 column and seps",{
    expect_equal(create_sep_list(0, 10),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 5 row, 9 column and seps",{
    expect_equal(create_sep_list(row_count=5,
                                 col_count=9,
                                 row_seps=c(1,3),
                                 col_seps=c(2,8)),
                 create_sep_list_answer_2)
})
 
context("remove_outliers_norm")
remove_outlier_norm_in_1 <- matrix(1:20, ncol=4)
remove_outlier_norm_out_1 <- matrix(c(rep(5,5),6:14,rep(15,6)), ncol=4)
remove_outlier_norm_in_2 <- matrix(c(1:15,
                                     c(-5,-4,3:13,21,26),
                                     1:15,
                                     1:15), ncol=4)
remove_outlier_norm_out_2 <- matrix(c(1:15,
                                    c(-.5,-.5,3:13,17.75,17.75),
                                    1:15,
                                    1:15), ncol=4)

test_that("remove_outliers_norm for hard threshold outside of values",{
    expect_equal(infercnv:::.remove_outliers_norm(remove_outlier_norm_in_1,
                                                      lower_bound=-1,
                                                      upper_bound=30),
                 remove_outlier_norm_in_1)
})
test_that("remove_outliers_norm for hard threshold within values",{
    expect_equal(infercnv:::.remove_outliers_norm(remove_outlier_norm_in_1,
                                                      lower_bound=5,
                                                      upper_bound=15),
                 remove_outlier_norm_out_1)

})
test_that("remove_outliers_norm for average bound",{
    expect_equal(infercnv:::.remove_outliers_norm(remove_outlier_norm_in_2,
                                                      out_method="average_bound"),
                 remove_outlier_norm_out_2)
})


context("order_reduce")
order_reduce_data_1 <- matrix(rep(1:10,2), ncol=2)
colnames(order_reduce_data_1) <- c("Sample_1","Sample_2")
row.names(order_reduce_data_1) <- paste("gene",1:10,sep="_") 
order_reduce_exp_1 <- matrix(rep(c(10,5,8,3,4,9,1,7,6,2),2), ncol=2)
row.names(order_reduce_exp_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")
colnames(order_reduce_exp_1) <- c("Sample_1","Sample_2")
order_reduce_pos_1 <- data.frame(chr=factor(c(1,1,2,2,3,3,4,4,5,5),levels=1:5),
                                 start=c(1,5,1,5,1,5,1,5,1,5),
                                 stop=c(4,9,4,9,4,9,4,9,4,9))
row.names(order_reduce_pos_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")
order_reduce_chr_1 <- data.frame(chr=factor(c(1,1,2,2,3,3,4,4,5,5),levels=1:5))
row.names(order_reduce_chr_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")

order_reduce_pos_2 <- data.frame(chr=factor(c(1,1,2,3,4,4),levels=1:4),
                                 start=c(1,5,5,5,1,5),
                                 stop=c(4,9,9,9,4,9))
row.names(order_reduce_pos_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")
order_reduce_exp_2 <- matrix(rep(c(10,5,3,9,1,7),2), ncol=2)
row.names(order_reduce_exp_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")
colnames(order_reduce_exp_2) <- c("Sample_1","Sample_2")
order_reduce_chr_2 <- data.frame(chr=factor(c(1,1,2,3,4,4),levels=1:4))
row.names(order_reduce_chr_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")

order_reduce_pos_3 <- data.frame(chr=c(1,1,2,3,4,4),
                                 start=c(1,5,5,5,1,5),
                                 stop=c(4,9,9,9,4,9))
row.names(order_reduce_pos_3) <- paste("GENE",c(10,5,3,9,1,7),sep="_")

test_that("order_reduce for NULL input.",{
    expect_equal(infercnv:::.order_reduce(NULL,NULL),
                 list(expr=NULL,order=NULL,chr_order=NULL))
})
test_that("order_reduce for happy path",{
    expect_equal(infercnv:::.order_reduce(order_reduce_data_1, order_reduce_pos_1),
                 list(expr=order_reduce_exp_1,
                      order=order_reduce_pos_1,
                      chr_order=order_reduce_chr_1))
})
test_that("order_reduce for happy path but dropping genes",{
    expect_equal(infercnv:::.order_reduce(order_reduce_data_1, order_reduce_pos_2),
                 list(expr=order_reduce_exp_2,
                      order=order_reduce_pos_2,
                      chr_order=order_reduce_chr_2))
})
test_that("order_reduce for no matching gene names",{
    expect_equal(infercnv:::.order_reduce(order_reduce_data_1, order_reduce_pos_3),
                 list(expr=NULL,
                      order=NULL,
                      chr_order=NULL))
})
