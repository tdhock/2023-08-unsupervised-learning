N.data <- 64L
min.seg.len <- 1L
max.segments <- 10L
cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
  N.data, min.seg.len, max.segments)
set.seed(1)
mean.vec <- seq(1, N.data)
mean.vec[ cost.dt[.N, seq(s2+1, N.data)] ] <- 0
data.vec <- rnorm(N.data, mean.vec, sd=0.01)
data.vec <- rep(c(0,1),l=N.data)
fit <- binsegRcpp::binseg_normal(data.vec, max.segments)
tree.list <- list(
  best=binsegRcpp::get_complexity_best_optimal_tree(cost.dt),
  empirical=binsegRcpp::get_tree_empirical(fit))
library(data.table)
tree.dt <- data.table(type=names(tree.list))[, {
  binsegRcpp::tree_layout(tree.list[[type]])
}, by=type]
total.dt <- tree.dt[, .(
  candidate.splits=sum(binsegRcpp::size_to_splits(size, min.seg.len))
), by=type]
join.dt <- total.dt[tree.dt, on="type"]
if(require(ggplot2)){
  ggplot()+
    facet_grid(. ~ type + candidate.splits, labeller=label_both)+
    geom_segment(aes(
      x, depth, 
      xend=parent.x, yend=parent.depth),
      data=join.dt)+
    geom_label(aes(
      x, depth, label=size),
      data=join.dt)+
    scale_y_reverse()
}


N.data <- 8L
min.seg.len <- 1L
max.segments <- 2L
cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
  N.data, min.seg.len, max.segments)
set.seed(1)
best.data <- if(max.segments<=4){
  c(2:max.segments, rep(0,N.data-max.segments+1))
}else{
  c(1:4, (4:1)*10)
}
tree.list <- list(
  worst=binsegRcpp::get_tree_empirical(binsegRcpp::binseg_normal(rep(0:1,l=N.data), max.segments)),
  best_guess=binsegRcpp::get_tree_empirical(binsegRcpp::binseg_normal(best.data, max.segments)),
  best_theory=binsegRcpp::get_complexity_best_optimal_tree(cost.dt),
  equal=binsegRcpp::get_tree_empirical(binsegRcpp::binseg_normal(1:N.data, max.segments)))
library(data.table)
tree.dt <- data.table(splits=names(tree.list))[, {
  binsegRcpp::tree_layout(tree.list[[splits]])
}, by=splits]
total.dt <- tree.dt[, .(
  candidates=sum(binsegRcpp::size_to_splits(size, min.seg.len))
), by=splits]
join.dt <- total.dt[tree.dt, on="splits"]
if(require(ggplot2)){
  ggplot()+
    facet_grid(. ~ splits + candidates, labeller=label_both)+
    geom_segment(aes(
      x, depth, 
      xend=parent.x, yend=parent.depth),
      data=join.dt)+
    geom_label(aes(
      x, depth, label=size),
      data=join.dt)+
    scale_y_reverse(breaks=seq(0,N.data))+
    scale_x_continuous("", breaks=NULL)+
    theme(panel.grid.minor=element_blank())
}




