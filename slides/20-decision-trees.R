N.col <- 2
N.row <- 200
set.seed(2)
X.mat <- matrix(runif(N.col*N.row), N.row, N.col)
hidden <- (X.mat[,1]<0.4 & X.mat[,2]<0.5) |
  (X.mat[,1] > 0.7 & X.mat[,2] > 0.4)
y.vec <- ifelse(hidden-0.5+rnorm(N.row,sd=0.5) < 0, 0, 1)
table(y.vec, hidden)
library(animint2)
library(data.table)

xy.dt <- data.table(y=y.vec, X.mat)
ggplot()+
  geom_point(aes(
    V1, V2, fill=factor(y)),
    shape=21,
    data=xy.dt)+
  scale_fill_manual(values=c(
    "0"="white",
    "1"="black"))+
  coord_equal()

candidate_loss <- function(indices){
  split.loss.list <- list()
  y.sub <- y.vec[indices]
  for(feature in 1:N.col){
    x.sub <- X.mat[indices,feature]
    ord.indices <- order(x.sub)
    x.ord <- x.sub[ord.indices]
    y.cum <- c(0,cumsum(y.sub[ord.indices]))
    Loss <- function(first,last){
      num.data <- last-first+1
      num.pos <- y.cum[last+1]-y.cum[first]
      num.neg <- num.data-num.pos
      prob.pos <- num.pos/num.data
      prob.neg <- 1-prob.pos
      ll <- function(num,prob)ifelse(num, num*log(prob), 0)
      -ll(num.pos,prob.pos)-ll(num.neg,prob.neg)
    }
    before.start <- 1
    after.end <- length(ord.indices)
    after.start <- seq(min.data.per.node+1, after.end-min.data.per.node+1)
    before.end <- after.start-1
    data.before <- x.ord[before.end]
    data.after <- x.ord[after.start]
    loss.split <- Loss(before.start, before.end)+Loss(after.start, after.end)
    loss.constant <- Loss(before.start, after.end)
    split.loss.list[[feature]] <- data.table(
      feature,
      data.before,
      split.point=(data.before+data.after)/2,
      data.after,
      loss.split,
      loss.constant,
      loss.diff=loss.split-loss.constant)
  }
  rbindlist(split.loss.list)[
  , is.min := loss.split == min(loss.split)
  ][]
}
min.data.per.node <- 5
new_node <- function(indices,parent=NULL,node.type="root"){
  node <- new.env()
  if(node.type=="root"){
    node$bounds <- data.table(
      V1min=0,V1max=1,
      V2min=0,V2max=1)
  }else{
    node$bounds <- parent$bounds
    node$bounds[[sprintf(
      "V%d%s",
      parent$best$feature,
      if(node.type=="lo")"max" else "min"
    )]] <- parent$best$split.point
  }
  node$id <- next.node.id
  node$parent <- parent
  node.list[[next.node.id]] <<- node
  next.node.id <<- next.node.id + 1L
  node$depth <- if(is.null(parent))0L else parent$depth+1L
  node$indices <- indices
  node$pred <- mean(y.vec[indices])
  if((#only compute loss if possible to split later.
    ! node$pred %in% c(0,1)    
  )&&(
    length(indices) >= min.data.per.node*2
  )){
    node$loss <- candidate_loss(indices)
    i.vec.min <- which(node$loss$is.min==TRUE)
    i.selected <- i.vec.min[1]
    node$loss[
    , loss.status := ifelse(is.min, "optimal", "sub-optimal")
    ][i.selected, loss.status := "can_split"]
    node$best <- node$loss[i.selected]
  }else{
    node$loss <- node$best <- NULL
  }
  node$terminal <- TRUE
  node$split <- function(){
    node$terminal <- FALSE
    is.lo <- node$best[, X.mat[indices,feature] < split.point]
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
  node$label_dt <- function(){
    data.table(
      terminal=node$terminal,
      prob1=node$pred,
      N=length(node$indices),
      feature=if(is.data.table(node$best))node$best$feature else NA,
      optimal=if(is.data.table(node$loss))sum(node$loss[, loss.constant!=0 & loss.status=="optimal"]) else 0,
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
tree_layout <- function(node.dt, space=0.5){
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
  px <- tree.dt[J(tree.dt$parent), x]
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

node.list <- list()
next.node.id <- 1L
dtree <- new_node(1:N.row)
current.loss <- dtree$best$loss.constant
can.split.list <- list(dtree)
rect.pred.list <- list()
train.pred.list <- list()
node.layout.list <- list()
tree.info.list <- list()
candidate.dt.list <- list()
for(iteration in 0:9){
  (train.pred.it <- data.table(
    iteration,
    row_id=1:nrow(xy.dt),
    xy.dt,
    pred.prob1=dtree$predict(X.mat)
  )[, let(
    pred.class = ifelse(pred.prob1<0.5, 0, 1),
    pred.prob.label = ifelse(y==1, pred.prob1, 1-pred.prob1)
  )][, let(
    correct = pred.class == y,
    log.loss=-log(pred.prob.label)
  )][])
  tree.info.list[[paste(iteration)]] <- data.table(
    iteration,
    loss=current.loss,
    train.log.loss=sum(train.pred.it$log.loss),
    nodes=length(node.list))
  train.pred.list[[paste(iteration)]] <- train.pred.it
  rect.pred.list[[paste(iteration)]] <- rbindlist(lapply(
    node.list, with, data.table(iteration, bounds, id, pred, terminal)))
  (node.parent.dt <- rbindlist(lapply(node.list, with, data.table(
    id, depth, parent=if(is.null(parent))NA_integer_ else parent$id,
    label_dt(), iteration
  ))))
  node.layout.list[[paste(iteration)]] <- tree_layout(
    node.parent.dt
  )[
  , label := ifelse(
    terminal,
    sprintf("p=%.2f", prob1),
    sprintf("X%d<%.2f", feature, split.point)
  )][]
  (can.split.best <- rbindlist(lapply(can.split.list, with, best)))
  split.i <- can.split.best[,which.min(loss.diff)]
  can.split.list[[split.i]]$loss[
    loss.status == "can_split",
    loss.status := "chosen_split"]
  candidate.dt.list[[paste(iteration)]] <- rbindlist(lapply(
    can.split.list, with, data.table(iteration, id, loss)))
  current.loss <- current.loss+can.split.best[split.i,loss.diff]
  can.split.list <- c(can.split.list[-split.i], can.split.list[[split.i]]$split())
}
(tree.info <- rbindlist(tree.info.list))
tree.info[, all.equal(loss, train.log.loss)]
(rect.pred <- rbindlist(rect.pred.list))
(train.pred <- rbindlist(train.pred.list))
(node.layout <- rbindlist(node.layout.list))
(candidate.dt <- rbindlist(candidate.dt.list)[
, Feature := paste0("X",feature)
])
ggplot()+
  facet_wrap("iteration", labeller=label_both,nrow=2)+
  geom_rect(aes(
    xmin=V1min, xmax=V1max,
    ymin=V2min, ymax=V2max,
    fill=pred),
    data=rect.pred)+
  geom_point(aes(
    V1, V2, fill=y, color=correct),
    shape=21,
    data=train.pred)+
  scale_color_manual(values=c("TRUE"="white","FALSE"="black"))+
  scale_fill_gradient2(low="blue",high="red",midpoint=0.5)+
  scale_x_continuous(breaks=seq(0,1,by=0.2))+
  scale_y_continuous(breaks=seq(0,1,by=0.2))+
  coord_equal()

ggplot()+
  facet_wrap("iteration", labeller=label_both,nrow=2)+
  geom_segment(aes(
    x, depth, 
    xend=parent.x, yend=parent.depth),
    data=node.layout)+
  scale_fill_gradient2(
    "Prob(y=1)",
    low="deepskyblue",high="red",midpoint=0.5,
    na.value="grey")+
  geom_label(aes(
    x, depth, label=label, fill=ifelse(terminal, prob1, NA)),
    data=node.layout)+
  scale_y_reverse()

rect.w <- 0.5
rect.h <- 0.3
node.layout[, Node := id]
candidate.dt[, let(
  Node = id,
  Split = sprintf("X%d<%f", feature, split.point)
)][]
last.pred <- rect.pred[
, Node := id
][iteration==max(iteration), .(Node,id,V1min,V1max,V2min,V2max)]
cand.join <- candidate.dt[
  last.pred, on="id"
][, `:=`(
  x=ifelse(feature==1, split.point, V1min),
  xend=ifelse(feature==1, split.point, V1max),
  y=ifelse(feature==1, V2min, split.point),
  yend=ifelse(feature==1, V2max, split.point)
)][]
(chosen.hilite <- candidate.dt[loss.status=="chosen_split"])
chosen.layout <- node.layout[
  chosen.hilite[, .(iteration,Node)],
  on=.(iteration,Node)]
## > RColorBrewer::brewer.pal(3,"Set1")
## [1] "#E41A1C" "#377EB8" "#4DAF4A"
chosen.color <- "#4DAF4A"#"green"
chosen.size <- 5
viz <- animint(
  title="Greedy decision tree learning algorithm for binary classification (Breiman's CART)",
  source="https://github.com/tdhock/2023-08-unsupervised-learning/blob/main/slides/20-decision-trees.R",
  video="https://vimeo.com/1048533733",
  loss=ggplot()+
    ggtitle("Select iteration")+
    theme_animint(width=300)+
    geom_point(aes(
      iteration, loss),
      help="One dot for each tree size",
      data=tree.info)+
    scale_x_continuous("Iteration/split", breaks=tree.info$iteration)+
    scale_y_continuous("Total logistic loss")+
    make_tallrect(tree.info, "iteration"),
  tree=ggplot()+
    ggtitle("Tree at selected iteration")+
    scale_x_continuous("<-yes(feature<threshold), no(feature>=threshold)->", breaks=NULL)+
    theme(legend.position="none")+
    geom_segment(aes(
      x, depth,
      key=paste(id,parent),
      xend=parent.x, yend=parent.depth),
      help="Segments connecting each parent node to its two child nodes",
      showSelected="iteration",
      data=node.layout)+
    scale_fill_gradient2(
      "Prob(y=1)",
      low="deepskyblue",high="red",midpoint=0.5,
      na.value="grey")+
    geom_rect(aes(
      xmin=x-rect.w, xmax=x+rect.w,
      ymin=depth+rect.h, ymax=depth-rect.h,
      key=id,
      fill=ifelse(terminal, prob1, NA)),
      help="One rect for each node in decision tree, gray for decision node, various shades of blue/red for prediction node",
      showSelected="iteration",
      color="transparent",
      data=node.layout)+
    geom_point(aes(
      x-rect.w*0.8, depth,
      key=1),
      help="Green dot shows node which will be split",
      data=chosen.layout,
      fill=chosen.color,
      size=chosen.size,
      color="black",
      showSelected="iteration")+
    geom_text(aes(
      x, depth+0.2, label=paste0(
        ifelse(optimal==0, "", "*"),
        label),
      key=id),
      help="Split rule (decision nodes) or proportion of positive labels (prediction nodes)",
      showSelected="iteration",
      data=node.layout)+
    geom_text(aes(
      x, depth-0.05, label=paste0("N=",N),
      key=id),
      help="Number of data points in node",
      showSelected="iteration",
      data=node.layout)+
    geom_rect(aes(
      xmin=x-rect.w, xmax=x+rect.w,
      ymin=depth+rect.h, ymax=depth-rect.h,
      key=id),
      help="Selected node shown with black border",
      color="black",
      fill="transparent",
      showSelected="iteration",
      color_off="transparent",
      clickSelects="Node",
      data=node.layout)+
    scale_y_reverse(),
  features=ggplot()+
    ggtitle("Decisions at selected iteration")+
    geom_rect(aes(
      xmin=V1min, xmax=V1max,
      ymin=V2min, ymax=V2max,
      key=Node,
      fill=pred),
      help="One rect for each node in decision tree",
      showSelected="iteration",
      color="transparent",
      data=rect.pred[terminal==TRUE])+
    geom_point(aes(
      V1, V2, fill=y,
      key=row_id,
      color=correct),
      help="One dot for each data point",
      shape=21,
      showSelected="iteration",
      data=train.pred)+
    scale_color_manual(values=c("TRUE"="white","FALSE"="black"))+
    scale_fill_gradient2(
      "Prob(y=1)",
        low="deepskyblue",high="red",midpoint=0.5,
        na.value="transparent")+
    scale_x_continuous("Feature X1", breaks=seq(0,1,by=0.2))+
    scale_y_continuous("Feature X2", breaks=seq(0,1,by=0.2))+
    geom_rect(aes(
      xmin=V1min, xmax=V1max,
      ymin=V2min, ymax=V2max,
      key=Node),
      help="Black border for selected node",
      clickSelects="Node",
      showSelected="iteration",
      color="black",
      fill="transparent",
      color_off="transparent",
      data=rect.pred)+
    geom_segment(aes(
      x, y,
      key=Split,
      xend=xend, yend=yend),
      data=cand.join,
      help="One grey segment for each split point in selected node",
      alpha=0.5,
      showSelected=c("iteration","Node"),
      clickSelects="Split")+
    coord_equal(),
  candidates=ggplot()+
    ggtitle("Loss decrease for selected iteration; select Node and Split")+
    theme_bw()+
    theme_animint(
      width=1100)+
    facet_grid(. ~ Feature, labeller=label_both)+
    scale_y_continuous(
      "<- more significant, Loss difference, less significant ->")+
    scale_x_continuous(
      "Threshold (split point)",
      breaks=seq(0,1,by=0.2))+
    geom_hline(aes(
      yintercept=loss.diff,
      key=1),
      help="Green horizontal line shows min loss in this iteration",
      data=chosen.hilite[,.(iteration,loss.diff)],
      color=chosen.color,
      showSelected="iteration")+
    geom_text(aes(
      split.point, loss.diff-3,
      hjust=ifelse(split.point<0.5, 0, 1),
      label=sprintf("_Diff=%.4f_", loss.diff),
      key=1),
      help="Green text shows min loss value in this iteration",
      data=chosen.hilite,
      color=chosen.color,
      showSelected="iteration")+
    geom_text(aes(
      split.point, -24,
      hjust=ifelse(split.point<0.5, 0, 1),
      label=sprintf("_Diff=%.4f_", loss.diff),
      key=paste(id, split.point)),
      help="Black text shows loss value of selected split",
      data=candidate.dt,
      showSelected=c("iteration","Node","Split"))+
    geom_tallrect(aes(
      key=paste(id, split.point),
      xmin=data.before,
      xmax=data.after),
      help="One rect for each split in selected node",
      data=candidate.dt,
      color=NA,
      alpha=0.5,
      clickSelects="Split",
      showSelected=c("iteration","Node"))+
    geom_segment(aes(
      split.point, -Inf,
      xend=split.point, yend=Inf,
      key=paste(id, split.point)),
      help="One segment for each split in selected node",
      data=candidate.dt,
      clickSelects="Split",
      alpha=0.5,
      showSelected=c("iteration","Node"))+
    geom_line(aes(
      split.point, loss.diff,
      group=id,
      key=id),
      help="One line for each feature/node that could be split in selected iteration",
      size=3,
      alpha=1,
      alpha_off=0.2,
      data=candidate.dt,
      showSelected="iteration",
      clickSelects="Node")+
    scale_fill_manual(values=c(
      "sub-optimal"="transparent",
      chosen_split=chosen.color,
      optimal="green",
      can_split="orange"))+
    geom_point(aes(
      split.point, loss.diff,
      fill=loss.status,
      key=paste(id, split.point)),
      help="One dot for each split in selected node",
      size=3,
      data=candidate.dt,
      alpha=1,
      alpha_off=1,
      color="black",
      color_off="transparent",
      clickSelects="Split",
      showSelected=c("iteration","Node")),
  duration=list(
    iteration=2000),
  out.dir="20-decision-trees"
)
viz

if(FALSE){
  animint2pages(viz, "2024-11-23-greedy-decision-tree")
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
