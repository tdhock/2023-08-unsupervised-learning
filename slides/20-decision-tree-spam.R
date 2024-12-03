library(animint2)
library(data.table)
if(!file.exists("spambase.zip")){
  download.file("https://archive.ics.uci.edu/static/public/94/spambase.zip","spambase.zip")
}
if(!file.exists("spambase.names")){
  unzip("spambase.zip")
}
(spam.names <- nc::capture_all_str(
  "spambase.names",
  "\n+",
  name="[^\\|]+?",
  ":"
)[
, abbrev := sub(".*_", "", name)
][])
spam.dt <- fread("spambase.data")
setnames(spam.dt, c(spam.names$abbrev, "spam"))
spam.dt[1]
N.folds <- 3
valid.fold <- 1
set.seed(1)
fold.vec <- spam.dt[, .(
  fold=sample(rep(1:N.folds,l=.N))
), by=spam]$fold
is.subtrain <- fold.vec!=valid.fold
y.col <- ncol(spam.dt)
X.mat <- as.matrix(spam.dt[,-y.col,with=FALSE])
y.vec <- spam.dt[[y.col]]
X.subtrain <- X.mat[is.subtrain,]
y.subtrain <- y.vec[is.subtrain]
table(y.vec, fold.vec)
N.features <- ncol(X.mat)
candidate_loss <- function(indices){
  split.loss.list <- list()
  y.sub <- y.subtrain[indices]
  for(feature in 1:N.features){
    x.sub <- X.subtrain[indices,feature]
    ord.indices <- order(x.sub)
    x.rle <- rle(x.sub[ord.indices])
    cum.lengths <- cumsum(x.rle$lengths)
    data.before <- x.rle$values[-length(x.rle$values)]
    data.after <- x.rle$values[-1]
    y.cum <- c(0,cumsum(y.sub[ord.indices])[cum.lengths])
    N.cum <- c(0,cum.lengths)
    Loss <- function(first,last){
      num.data <- N.cum[last+1]-N.cum[first]
      num.pos <- y.cum[last+1]-y.cum[first]
      num.neg <- num.data-num.pos
      prob.pos <- num.pos/num.data
      prob.neg <- 1-prob.pos
      ll <- function(num,prob)ifelse(num, num*log(prob), 0)
      -ll(num.pos,prob.pos)-ll(num.neg,prob.neg)
    }
    before.start <- 1
    after.end <- length(x.rle$values)
    after.start <- seq(2, after.end)
    before.end <- after.start-1
    loss.split <- Loss(before.start, before.end)+Loss(after.start, after.end)
    loss.constant <- Loss(before.start, after.end)
    split.loss.list[[feature]] <- data.table(
      feature,
      N.before=N.cum[before.end+1]-N.cum[before.start],
      N.after=N.cum[after.end+1]-N.cum[after.start],
      data.before,
      split.point=(data.before+data.after)/2,
      data.after,
      loss.split,
      loss.constant,
      loss.diff=loss.split-loss.constant)
  }
  rbindlist(split.loss.list)
}
min.data.per.node <- 20
new_node <- function(indices,parent=NULL,node.type="root"){
  node <- new.env()
  node$lo <- node$hi <- NULL
  node$id <- next.node.id
  node$parent <- parent
  node.list[[next.node.id]] <<- node
  next.node.id <<- next.node.id + 1L
  node$depth <- if(is.null(parent))0L else parent$depth+1L
  node$indices <- indices
  node$pred <- mean(y.subtrain[indices])
  node$feasible <- candidate_loss(indices)[
    N.before >= min.data.per.node & N.after >= min.data.per.node & loss.constant != 0
  ]
  node$best <- if(nrow(node$feasible)){
    node$feasible[
    , is.min := loss.split == min(loss.split)
    ][]
    i.vec.min <- which(node$feasible$is.min==TRUE)
    i.selected <- i.vec.min[1]
    node$feasible[
    , loss.status := ifelse(is.min, "optimal", "sub-optimal")
    ][i.selected, loss.status := "can_split"]
    node$feasible[i.selected]
  }
  node$terminal <- TRUE
  node$split <- function(){
    node$terminal <- FALSE
    is.lo <- node$best[, X.subtrain[indices,feature] < split.point]
    logical.list <- list(lo=is.lo, hi=!is.lo)
    out.list <- list()
    for(child.name in names(logical.list)){
      is.child <- logical.list[[child.name]]
      child.indices <- indices[is.child]
      node[[child.name]] <- new_node(child.indices,node,child.name)
      if(is.data.table(node[[child.name]]$best)){
        out.list[[child.name]] <- node[[child.name]]
      }
    }
    out.list
  }
  node$prune_children <- function(){
    node$terminal <- TRUE
    node$lo <- node$hi <- NULL
    if(isTRUE(node$parent$lo$terminal) && isTRUE(node$parent$hi$terminal))node$parent
  }
  node$label_dt <- function(){
    data.table(
      terminal=node$terminal,
      prob1=node$pred,
      N=length(node$indices),
      feature=if(is.data.table(node$best))node$best$feature else NA,
      optimal=if(nrow(node$feasible))sum(node$feasible[, loss.constant!=0 & loss.status=="optimal"]) else 0,
      split.point=if(is.data.table(node$best))node$best$split.point else NA)
  }
  node$expr <- function(){
    if(node$terminal)node$pred
    else substitute(
      ifelse(X[,FEAT]<VAL, LO, HI), with(node$best, list(
        FEAT=feature, VAL=split.point,
        LO=node$lo$expr(),
        HI=node$hi$expr())))
  }
  node$predict <- function(X){
    eval(node$expr())
  }
  node
}
tree_layout <- function(node.dt){
  id.tab <- table(node.dt$id)
  stopifnot(all(id.tab==1))
  tree.dt <- data.table(node.dt, key="id")[, x := NA_real_]
  for(d in unique(tree.dt$depth)){
    if(d==0)tree.dt[depth==0, x := 0] else{
      d.nodes <- tree.dt[depth==d]
      px <- tree.dt[J(d.nodes$parent), x, nomatch=0L]
      d.nodes[, parent.x := px]
      ord.nodes <- d.nodes[order(parent.x)]
      new.x <- ord.nodes[, qp.x(parent.x,parent.x+space,parent.x-space)]
      tree.dt[J(ord.nodes$id), x := new.x]
    }
  }
  px <- tree.dt[J(tree.dt$parent), x, nomatch=NA]
  tree.dt[, let(
    parent.x = px,
    parent.depth = ifelse(is.na(parent), NA, depth-1)
  )][]
}
qp.x <- function
### Solve quadratic program to find x positions.
(target, y.up, y.lo){
  k <- length(target)
  D <- diag(rep(1, k))
  Ik <- diag(rep(1, k - 1))
  A <- rbind(0, Ik) - rbind(Ik, 0)
  b0 <- (y.up - target)[-k] + (target - y.lo)[-1]
  sol <- quadprog::solve.QP(D, target, A, b0)
  sol$solution
}
current_tree_layout <- function(){
  (node.parent.dt <- rbindlist(lapply(node.list, function(x)if(!is.null(x))with(x, data.table(
    id, depth, parent=if(is.null(parent))NA_integer_ else parent$id,
    label_dt(), iteration
  ))))[, let(
    label = ifelse(
      terminal,
      sprintf("p=%.2f", prob1),
      sprintf("%s<%.2f", colnames(X.mat)[feature], split.point)),
    Nlab = paste0(
      ifelse(optimal==0, "", "*"),
      "N=",N)
  )][
  , space := pmax(nchar(label),nchar(Nlab))
  ][])
  tree_layout(node.parent.dt)
}
current_set_cost <- function(){
  (train.pred.it <- data.table(
    iteration,
    row_id=1:nrow(spam.dt),
    set.name=ifelse(is.subtrain, "subtrain", "validation"),
    y=y.vec,
    pred.prob1=dtree$predict(X.mat)
  )[, let(
    pred.class = ifelse(pred.prob1<0.5, 0, 1),
    pred.prob.label = ifelse(y==1, pred.prob1, 1-pred.prob1)
  )][, let(
    correct = pred.class == y,
    log.loss=-log(pred.prob.label)
  )][])
  loss.long <- rbind(
    data.table(set.name="diff", loss=current.loss, error.percent=NA),
    train.pred.it[, .(
      loss=sum(log.loss),
      error.percent=100*mean(!correct)
    ), by=set.name])
  data.table(
    iteration, leaves=current.leaves, loss.long)
}

node.list <- list()
next.node.id <- 1L
dtree <- new_node(seq_along(y.subtrain))
current.loss <- dtree$best$loss.constant
can.split.list <- list(dtree)
node.layout.list <- list()
tree.info.list <- list()
candidate.dt.list <- list()
iteration <- 0
current.leaves <- 1
done <- FALSE
while(!done){
  tree.info.list[[paste(iteration)]] <- current_set_cost()
  node.layout.list[[paste(iteration)]] <- data.table(
    pass="forward",
    current_tree_layout())
  if(
    length(can.split.list)==0
    ##|| iteration==20
    ##|| any(is.infinite(loss.long$loss))
  ){
    done <- TRUE
  }else{
    (can.split.best <- rbindlist(lapply(can.split.list, with, best)))
    split.i <- can.split.best[,which.min(loss.diff)]
    can.split.list[[split.i]]$feasible[
      loss.status == "can_split",
      loss.status := "chosen_split"]
    candidate.dt.list[[paste(iteration)]] <- rbindlist(lapply(
      can.split.list, with, data.table(iteration, id, feasible)))
    current.loss <- current.loss+can.split.best[split.i,loss.diff]
    can.split.list <- c(can.split.list[-split.i], can.split.list[[split.i]]$split())
    print(iteration <- iteration+1)
    current.leaves <- current.leaves+1
  }
}
max.it <- iteration
(tree.info <- rbindlist(tree.info.list))
tree.info.wide <- dcast(tree.info, iteration ~ set.name, value.var="loss")
tree.info.wide[, all.equal(diff, subtrain)]

initial.prune.dt <- rbindlist(lapply(node.list, with, data.table(
  id, can_prune=isTRUE(lo$terminal)&&isTRUE(hi$terminal))))
can.prune.list <- node.list[initial.prune.dt[can_prune==TRUE]$id]
prune.cost.list <- tree.info.list[length(tree.info.list)]
chosen.hilite.prune.list <- list()
while(length(can.prune.list)){
  can.prune.dt <- rbindlist(lapply(can.prune.list, with, data.table(id, best)))
  prune.i <- can.prune.dt[,which.max(loss.diff)]
  parent_of_pruned <- can.prune.list[[prune.i]]
  chosen.hilite.prune.list[[paste(iteration)]] <- data.table(
    iteration, Node=parent_of_pruned$id)
  iteration <- iteration+1
  current.leaves <- current.leaves-1
  node.list[parent_of_pruned$lo$id] <- list(NULL)
  node.list[parent_of_pruned$hi$id] <- list(NULL)
  can.prune.list <- c(
    can.prune.list[-prune.i],
    parent_of_pruned$prune_children())
  node.layout.list[[paste(iteration)]] <- data.table(
    pass="pruning",
    current_tree_layout())
  current.loss <- current.loss-can.prune.dt[prune.i, loss.diff]
  prune.cost.list[[paste(iteration)]] <- current_set_cost()
}
(prune.cost <- rbindlist(prune.cost.list))

best.valid <- tree.info[set.name=="validation"][which.min(error.percent)]
ggplot()+
  geom_vline(aes(
    xintercept=iteration),
    data=best.valid)+
  geom_hline(aes(
    yintercept=error.percent),
    data=best.valid)+
  geom_point(aes(
    iteration, error.percent, color=set.name),
    data=tree.info)+
  geom_line(aes(
    iteration, error.percent, color=set.name),
    data=tree.info)+
  scale_y_log10()

both.cost <- melt(rbind(
  data.table(pass="forward", tree.info),
  data.table(pass="pruning", prune.cost)
)[set.name!='diff'], measure.vars=c("loss","error.percent"))
both.selection <- both.cost[
  variable=="loss" & set.name=="subtrain",
  penaltyLearning::modelSelection(.SD, "value", "leaves"),
  by=pass]
both.points <- both.cost[both.selection[,.(pass,leaves)], on=.(pass,leaves)]
(both.join <- both.selection[
,.(pass,leaves,min.log.lambda,max.log.lambda)
][
  both.cost[pass=="pruning"], on=.(pass,leaves)
][
, selected := !is.na(min.log.lambda)
][])
both.best <- both.points[
  set.name=="validation" & pass=="pruning", .SD[min(value)==value][which.min(leaves)], by=variable
][, .(iteration, leaves, variable, value, set.name)]
ggplot()+
  geom_vline(aes(
    xintercept=leaves),
    data=both.best)+
  geom_hline(aes(
    yintercept=value),
    data=both.best)+
  geom_point(aes(
    leaves, value, color=pass),
    shape=21,
    data=both.points)+
  geom_line(aes(
    leaves, value, color=pass, group=paste(pass,set.name)),
    data=both.cost)+
  scale_y_log10()+
  facet_grid(variable ~ set.name, labeller=label_both, scales="free")


ggplot()+
  geom_vline(aes(
    xintercept=iteration),
    data=both.best)+
  geom_vline(aes(
    xintercept=iteration),
    data=data.table(iteration=max.it))+
  geom_hline(aes(
    yintercept=value),
    data=both.best)+
  geom_point(aes(
    iteration, value, color=set.name),
    shape=21,
    data=both.points)+
  geom_line(aes(
    iteration, value, color=set.name, group=paste(pass,set.name)),
    data=both.cost)+
  scale_y_log10()+
  facet_grid(variable ~ ., labeller=label_both, scales="free")+
  scale_x_continuous(
    breaks=seq(0,200,by=20))

(node.layout <- rbindlist(node.layout.list))
(candidate.dt <- rbindlist(candidate.dt.list)[
, Feature := paste0("X",feature)
])

last.layout <- node.layout[iteration==max.it]
rect.h <- 0.35
rect.w.fac <- 0.9
ggplot()+
  ggtitle("Tree at last iteration")+
  scale_x_continuous("<-yes(feature<threshold), no(feature>=threshold)->", breaks=NULL)+
  geom_segment(aes(
    x, depth,
    key=paste(id,parent),
    xend=parent.x, yend=parent.depth),
    data=last.layout)+
  scale_fill_gradient2(
    "Prob(y=1)",
    low="deepskyblue",high="red",midpoint=0.5,
    na.value="grey")+
  geom_rect(aes(
    xmin=x-space, xmax=x+space,
    ymin=depth+rect.h, ymax=depth-rect.h,
    key=id,
    fill=ifelse(terminal, prob1, NA)),
    color="black",
    data=last.layout)+
  geom_text(aes(
    x, depth+0.1, label=label,
    key=id),
    data=last.layout)+
  geom_text(aes(
    x, depth-0.1, label=Nlab,
    key=id),
    data=last.layout)+
  scale_y_reverse()

node.layout[, Node := id]
candidate.dt[, let(
  Node = id,
  Split = sprintf("X%d<%f", feature, split.point)
)][]
(chosen.hilite <- rbind(
  candidate.dt[loss.status=="chosen_split",  .(iteration,Node)],
  rbindlist(chosen.hilite.prune.list)))
chosen.layout <- node.layout[
  chosen.hilite,
  on=.(iteration,Node)]
## > RColorBrewer::brewer.pal(3,"Set1")
## [1] "#E41A1C" "#377EB8" "#4DAF4A"
chosen.color <- "#4DAF4A"#"green"
chosen.size <- 5
tree.text.size <- 10
it.vline <- data.table(iteration=max.it)
both.vline <- rbind(
  it.vline[,.(variable="error.percent", iteration, label="largest tree")],
  both.best[,.(variable, iteration, label = paste("best validation", variable))])
vline.color <- "grey50"
both.join.sel <- both.join[!is.na(min.log.lambda)]
pass.cost <- dcast(
  both.cost[set.name=="subtrain"], leaves ~ pass, first, value.var="value"
)[, better := ifelse(forward<pruning, "growing", ifelse(pruning<forward, "pruning", "same"))][]
u.it.leaves <- unique(both.cost[,.(iteration,leaves)])
(leaves.hilite <- u.it.leaves[u.it.leaves,on="leaves"][order(iteration)][
, rank := rank(i.iteration)
##, by=iteration # for smooth transitions.
][pass.cost, on="leaves"])
viz <- animint(
  title="Cross-validation for Breiman's CART algorithm on SPAM data",
  source="https://github.com/tdhock/2023-08-unsupervised-learning/blob/main/slides/20-decision-tree-spam.R",
  selection=ggplot()+
    ggtitle("Select iteration via pruning penalty")+
    theme_bw()+
    theme(legend.position="none")+
    theme_animint(height=300)+
    scale_y_log10("")+
    geom_segment(aes(
      min.log.lambda, value,
      xend=max.log.lambda, yend=value,
      color=set.name),
      data=both.join.sel,
      showSelected="set.name",
      size=4)+
    ## geom_segment(aes(
    ##   mid.log.lambda, -Inf,
    ##   xend=mid.log.lambda, yend=Inf),
    ##   size=4,
    ##   clickSelects="iteration",
    ##   data=both.join.sel[, mid.log.lambda := ifelse(
    ##     min.log.lambda == -Inf, max.log.lambda-1,
    ##     ifelse(max.log.lambda == Inf, min.log.lambda+1, (min.log.lambda+max.log.lambda)/2))],
    ##   alpha=0.5)+
    scale_x_continuous("log(alpha = pruning penalty parameter)")+
    geom_tallrect(aes(
      xmin=min.log.lambda, xmax=max.log.lambda,
      ymin=0, ymax=Inf),
      data=both.join.sel,
      alpha=0.2,
      clickSelects="iteration")+
    facet_grid(variable ~ ., scales="free"),
  loss=ggplot()+
    theme_bw()+
    ggtitle("Select iteration")+
    theme_animint(width=800, height=300)+
    geom_vline(aes(
      xintercept=iteration),
      color=vline.color,
      data=both.best)+
    geom_vline(aes(
      xintercept=iteration),
      color=vline.color,
      data=it.vline)+
    geom_text(aes(
      iteration, Inf, label=label),
      hjust=0,
      data=both.vline)+
    geom_hline(aes(
      yintercept=value),
      color=vline.color,
      data=both.best)+
    geom_vline(aes(
      xintercept=iteration,
      linetype=better),
      data=unique(leaves.hilite[order(iteration),.(iteration,better)]),
      color="blue",
      alpha=0.2)+
    scale_linetype_manual(values=c(
      same=0,
      growing=1,
      pruning=3))+
    geom_tallrect(aes(
      xmin=i.iteration-0.5,
      xmax=i.iteration+0.5,
      key=rank),
      data=leaves.hilite,
      fill="black",
      alpha=0.2,
      color="transparent",
      showSelected="iteration")+
    geom_line(aes(
      iteration, value, color=set.name, group=paste(pass,set.name)),
      data=both.cost)+
    geom_point(aes(
      iteration, value, color=set.name, fill=selected),
      shape=21,
      data=both.join)+
    scale_y_log10("")+
    scale_fill_manual(values=c(
      "TRUE"="black",
      "FALSE"="transparent"))+
    make_tallrect(both.cost, "iteration")+
    facet_grid(variable ~ ., scales="free")+
    scale_x_continuous(
      breaks=seq(0,200,by=20)),
  tree=ggplot()+
    ggtitle("Tree at selected iteration")+
    theme_bw()+
    scale_x_continuous("<-yes(feature<threshold), no(feature>=threshold)->", breaks=NULL)+
    ##theme(legend.position="none")+
    theme_animint(width=1200, height=500)+
    geom_segment(aes(
      x, depth-rect.h,
      key=paste(id,parent),
      xend=parent.x, yend=parent.depth+rect.h),
      showSelected="iteration",
      size=1,
      data=node.layout[is.finite(parent.x)])+
    scale_fill_gradient2(
      "Prob(y=1)",
      low="deepskyblue",high="red",midpoint=0.5,
      na.value="grey")+
    geom_rect(aes(
      xmin=x-rect.w.fac*space, xmax=x+rect.w.fac*space,
      ymin=depth+rect.h, ymax=depth-rect.h,
      key=id,
      fill=ifelse(terminal, prob1, NA)),
      showSelected="iteration",
      color="transparent",
      data=node.layout)+
    geom_point(aes(
      x-space*0.8, depth,
      key=1),
      data=chosen.layout,
      fill=chosen.color,
      size=chosen.size,
      color="black",
      showSelected="iteration")+
    geom_text(aes(
      x, depth+0.3, label=label,
      key=id),
      size=tree.text.size,
      showSelected="iteration",
      data=node.layout)+
    geom_text(aes(
      x, depth-0.05, label=Nlab,
      key=id),
      size=tree.text.size,
      showSelected="iteration",
      data=node.layout)+
    coord_cartesian(expand=FALSE)+
    ## geom_rect(aes(
    ##   xmin=x-space, xmax=x+space,
    ##   ymin=depth+rect.h, ymax=depth-rect.h,
    ##   key=id),
    ##   color="black",
    ##   fill="transparent",
    ##   showSelected="iteration",
    ##   color_off="transparent",
    ##   clickSelects="Node",
    ##   data=node.layout)+
    scale_y_reverse(),
  duration=list(
    iteration=1000),
  out.dir="20-decision-tree-spam"
)
viz

if(FALSE){
  animint2pages(viz, "2024-11-27-decision-tree-spam")
}

## in feature plot, add aes(color=correct) to geom_point, by computing
## predicted probabilities for each data point in the train set
## (threshold at 0.5 to obtain predicted class, then define corect if
## predicted class same as label/y class). Also use the predicted
## probabilities to compute the total log likelihood over all data in
## the train set, and compare that to the number that comes from the
## iterative learning algorithm loss update rule.

## in feature plot, use geom_rect instead of geom_tile to display
## predictions.

## in tree and feature plots, for geom_rect (nodes in tree, and
## predictions in feature), remove clickSelects="Node" and set
## color="transparent" but keep aes(fill). Add a new geom_rect with
## fill="transparent" and clickSelects="Node" as the last layer, so
## that clicking anywhere in the region will change the selection
## (even over text/points).

## in candidates plot, add green geom_point and geom_hline to
## emphasize best split for each iteration. add green geom_point to
## tree plot.

## in candidates plot, add details for each candidate split point: how
## many data/pos/neg before/after split?

## in candidates plot, add geom_point and geom_segment for each split
## point, with clickSelects="Split" (in addition to geom_tallrect).
