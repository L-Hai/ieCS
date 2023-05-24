###############################################################################
### Define class

### Pair-wise similarity score object
setOldClass("pheatmap")
setOldClass("phylo")
scoreCSClass <- setClass(
  "scoreCSObject",
  slots = c(
    marker.table = "data.frame",
    pair.score.matrix = "matrix",
    p1 = "pheatmap",
    edges = "data.frame",
    nodes = "data.frame"
  )
)


### Tree aggregation method object
CS <- setClass(
  "ClusterSimilarity",
  slots = c(
    score = "vector",
    cluster = "list",
    gene = "list",
    text = "list"
  )
)
CSCoreClass <- setClass("CSCoreFuncObject",
                        slots = c(newickCS = "ClusterSimilarity",
                                  treeObject = "phylo"))

###############################################################################
### Calculate pair-wise similarity scores

# @scoreMode
# Integer 1,Query x Query; 2,(Query + Reference) x (Query + Reference); 3,Query x Reference

scoreCS <- function(scoreMode) {
  ModeString <-
    switch(
      scoreMode,
      "Query x Query",
      "(Query + Reference) x (Query + Reference)",
      "Query x Reference"
    )
  withProgress(
    message = paste0('Scoring cluster similarity on ', ModeString),
    value = 0.1,
    {
      g1 <- read.csv('ieCS.Que.Markers.csv', stringsAsFactor = F)
      if (scoreMode == 1) {
        g2 <- g1
        marker.table <- g1
      }
      if (scoreMode != 1) {
        g2 <- read.csv('ieCS.Ref.Markers.csv', stringsAsFactor = F)
        g2$tag <- 'Reference'
        g2 <- g2[, c('cluster', 'gene', 'order', 'tag')]
        g1$tag <- 'Query'
        g1 <- g1[, c('cluster', 'gene', 'order', 'tag')]
        marker.table <- rbind(g1, g2)
        if (scoreMode == 2) {
          g1 <- marker.table
          g2 <- marker.table
        }
      }
      pair.rank <- list()
      pair.score <- list()
      graph.list <-
        data.frame(
          from = character(),
          to = character(),
          weight = character(),
          stringsAsFactors = FALSE
        )
      g1.c <- unique(g1$cluster)
      g2.c <- unique(g2$cluster)
      pair.score.matrix <-
        matrix(0, nrow = length(g1.c), ncol = length(g2.c))
      i = 1
      nall <- 0.7 / length(g1.c)
      for (c1 in g1.c) {
        incProgress(nall)
        g1.g <- g1[g1$cluster == c1, 'gene']
        j = 1
        for (c2 in g2.c) {
          g2.g <- g2[g2$cluster == c2, 'gene']
          g.g <- intersect(g1.g, g2.g)
          if (length(g.g) > 0) {
            g1.r <- match(g.g, g1.g)
            g2.r <- match(g.g, g2.g)
            pair.rank[[paste0(c1, '_', c2)]] <-
              paste0(g1.r, '_', g2.r)
            ascore <- floor(sum(50 / (g1.r + g2.r)))
            pair.score[[paste0(c1, '_', c2)]] <- ascore
            pair.score.matrix[i, j] <- ascore
            graph.list <-
              rbind(graph.list,
                    data.frame(
                      from = c1,
                      to = c2,
                      weight = ascore
                    ))
          } else{
            pair.rank[[paste0(c1, '_', c2)]] <- 'NA'
            pair.score[[paste0(c1, '_', c2)]] <- 0
            graph.list <-
              rbind(graph.list, data.frame(
                from = c1,
                to = c2,
                weight = 0
              ))
          }
          j = j + 1
        }
        i = i + 1
      }
      pair.score <- unlist(pair.score)
      pair.score <- pair.score[order(pair.score, decreasing = T)]

      rownames(pair.score.matrix) <- g1.c
      colnames(pair.score.matrix) <- g2.c

      pair.rank.l <-
        lapply(pair.rank, function(x)
          x = paste(x, collapse = ";"))
      p1 <- pheatmap::pheatmap(pair.score.matrix)
      incProgress(0.1)
      if(scoreMode !=3){
      graph.df <- graph.list
      nodes <-
        data.frame(
          'id' = 1:length(g1.c),
          'label' = g1.c,
          stringsAsFactors = FALSE
        )
      graph.df[, 1] <-
        plyr::mapvalues(graph.df[, 1], nodes$label, nodes$id, warn_missing = FALSE)
      graph.df[, 2] <-
        plyr::mapvalues(graph.df[, 2], nodes$label, nodes$id, warn_missing = FALSE)

      edges <- graph.df
      reor <- data.frame(
        from = character(),
        to = character(),
        weight = character(),
        stringsAsFactors = FALSE
      )
      for (i in 1:nrow(edges)) {
        x <- edges[i, ]
        if (as.numeric(x[1]) > as.numeric(x[2])) {
          temp <- x[2]
          x[2] <- x[1]
          x[1] <- temp
        }
        if (x[1] == x[2]) {
          next
        }
        reor <- rbind(reor, x)
      }
      reor <- unique(reor)
      edges <- reor
    }}
  )
  if(scoreMode == 3){
    edges <- data.frame()
    nodes <- data.frame()
  }
  scoreCSObj <- scoreCSClass(
    marker.table = marker.table,
    pair.score.matrix = pair.score.matrix,
    p1 = p1,
    edges = edges,
    nodes = nodes
  )
  return(scoreCSObj)
}

###############################################################################
### Identify overlapping gene of cell clusters in each supercluster

# @marker.table
# Markers table
# @supercluster.composition
# Supercluster composition
# Return
# overlapping gene of cell clusters in each supercluster

findOverlapGene <-
  function(marker.table, supercluster.composition) {
    overlappingGeneList <- list()
    for (i in 1:length(supercluster.composition)) {
      oneSuperCluster <- supercluster.composition[i]
      if (length(oneSuperCluster[[1]]) > 1) {
        cluster1 <- oneSuperCluster[[1]][1]
        cluster2 <- oneSuperCluster[[1]][2]
        cluster1.g <-
          marker.table[marker.table$cluster == cluster1, 'gene']
        cluster2.g <-
          marker.table[marker.table$cluster == cluster2, 'gene']
        overlappingGene <- intersect(cluster1.g, cluster2.g)
        if (length(oneSuperCluster[[1]]) > 2) {
          for (j in 3:length(oneSuperCluster[[1]])) {
            cluster.more.g <-
              marker.table[marker.table$cluster == oneSuperCluster[[1]][j], 'gene']
            overlappingGene <-
              intersect(overlappingGene, cluster.more.g)
          }
        }
        if (length(overlappingGene) > 0) {
          overlappingGeneList[[names(oneSuperCluster)]] <- overlappingGene
        }
      }
    }
    return(overlappingGeneList)
  }

###############################################################################
### Identify superclusters with tree aggregation method

### Find out sets of similar clusters which past the minimum similarity score cutoff
# @df
# Pair-wise cluster similarity matrix
# @n
# minimum similarity score

findSimilarClusterSet <- function(df, n) {
  df.binary <- df > n
  df.index <- which(df.binary == TRUE, arr.ind = T)
  df.index <- as.data.frame(df.index)
  df.index$rowname <- rownames(df)[df.index$row]
  df.index$colname <- gsub('\\.', '-', colnames(df)[df.index$col])
  similarClusterSet <- list()
  for (i in unique(df.index$colname)) {
    similarClusterSet[[i]] <- df.index[df.index$colname == i, 'rowname']
  }
  similarClusterSet <- unique(similarClusterSet)
  return(similarClusterSet)
}

### Perform tree aggregation of cell clusters
# @marker.table
# Markers table
# @seqn
# Integer sequence from user-defined minimum similarity cutoff to maximum similarity score
# @df
# Pair-wise cluster similarity matrix
# Return
# Tree aggregation result object

CSCoreFunc <- function(marker.table, seqn, df) {
  #== Identify sets of similar cluster in different degree of similarity
  comclist <- list()
  for (n in seqn) {
    comclist[[as.character(n)]] <- findSimilarClusterSet(df, n)
  }
  unique.name <-
    names(comclist)[!duplicated(comclist)]
  un.comclist <- unique(comclist)
  names(un.comclist) <- unique.name

  condi <- list()
  for (i in 1:length(un.comclist)) {
    ll1 <- un.comclist[i]
    ll2n <- 1
    for (ll2 in ll1[[1]]) {
      if (length(ll2) > 1) {
        condin <- paste(names(ll1), ll2n, sep = '_')
        condi[[condin]] <- ll2
      }
      ll2n <- ll2n + 1
    }
  }

  #== Find overlapping genes
  condi <- rev(condi) # From high score to low score
  un.condi.n <-
    names(condi)[!duplicated(condi)]
  un.condi <- unique(condi)
  names(un.condi) <- un.condi.n
  un.condi.fi <-
    findOverlapGene(marker.table, un.condi)
  un.condi.fi.n <- names(un.condi.fi)
  un.condi.fi.score <-
    unlist(lapply(strsplit(un.condi.fi.n, '_'), function(x) {
      x[1]
    }))
  un.condi.fi.cluster <-
    un.condi[un.condi.fi.n]

  #== Result object of sets of similar clusters
  multi <-
    CS(score = un.condi.fi.score,
       cluster = un.condi.fi.cluster,
       gene = un.condi.fi)

  #== Construct tree structure
  cluster.all <- rownames(df)
  cluster.all.score <-
    floor(diag(as.matrix(df)))
  cluster.all.text <-
    paste(cluster.all, cluster.all.score, sep = ':')
  newick.text <- list()
  newick.element <- list()

  #== Start = identify supercluster
  for (i in 1:length(multi@cluster)) {
    acluster <- multi@cluster[[i]]
    score <- multi@score[i]
    if (length(newick.element) == 0) {
      newick.element[[1]] <-  acluster
      acluster.text <-
        plyr::mapvalues(acluster,
                        cluster.all,
                        cluster.all.text,
                        warn_missing = FALSE)
      text <-
        paste0('(', paste(acluster.text, collapse = ','), '):', score)
      newick.text[[1]] <- text
    } else{
      #== Check inclusion of sets
      uniqueCount <- 0
      containElement <- c()
      for (j in 1:length(newick.element)) {
        element <- newick.element[[j]]
        inter <-
          intersect(element, acluster)
        if (length(inter) == 0) {
          uniqueCount <- uniqueCount + 1
        }
        if (length(inter) == length(element)) {
          containElement <- c(containElement, j)
        }
      }

      #== Build Newick tree format
      if (uniqueCount == length(newick.element)) {
        newick.element[[j + 1]] <-  acluster
        acluster.text <-
          plyr::mapvalues(acluster,
                          cluster.all,
                          cluster.all.text,
                          warn_missing = FALSE)
        text <-
          paste0('(', paste(acluster.text, collapse = ','), '):', score)
        newick.text[[j + 1]] <- text
      }
      if (length(containElement) > 0) {
        elementList <- newick.element[containElement]
        adddiff <-
          setdiff(acluster, unlist(elementList))
        acluster.text <-
          plyr::mapvalues(adddiff,
                          cluster.all,
                          cluster.all.text,
                          warn_missing = FALSE)
        if (length(intersect(adddiff, unlist(newick.element))) == 0) {
          if (length(containElement) == 1) {
            newick.element[[containElement]] <- acluster
            text <-
              paste0('(',
                     newick.text[[containElement]],
                     ',',
                     paste(acluster.text, collapse = ','),
                     '):',
                     score)
            newick.text[[containElement]] <-
              text
          }
          else{
            newick.element <- newick.element[-containElement]
            newick.element[[length(newick.element) + 1]] <-
              acluster
            if (length(adddiff) == 0) {
              text <-
                paste0('(',
                       paste(newick.text[containElement], collapse = ','),
                       '):',
                       score)
            } else{
              text <-
                paste0(
                  '(',
                  paste(newick.text[containElement], collapse = ','),
                  ',',
                  paste(acluster.text, collapse = ','),
                  '):',
                  score
                )
            }
            newick.text <-
              newick.text[-containElement]
            newick.text[[length(newick.text) + 1]] <-
              text
          }
        }
      }
    }
  }
  #== End = identify supercluster

  #== Combine subtree (supercluster)
  single <-
    setdiff(cluster.all, unlist(newick.element))
  single.text <-
    plyr::mapvalues(single, cluster.all, cluster.all.text, warn_missing = FALSE)
  newick.text <-
    append(newick.text, single.text)
  newick.element <-
    append(newick.element, single)
  treeText <-
    paste0('(', paste(newick.text, collapse = ','), ');')
  t1 <- ape::read.tree(text = treeText)
  newick.score <-
    strsplit(unlist(newick.text), ':')
  newick.score <-
    unlist(lapply(newick.score, function(x) {
      x[length(x)]
    }))
  names(newick.element) <-
    paste(newick.score, 1:length(newick.score), sep = '_')
  newick.gene <-
    findOverlapGene(marker.table, newick.element)

  #== Build result object
  newickCS <-
    CS(
      score = newick.score,
      cluster = newick.element,
      gene = newick.gene,
      text = newick.text
    )
  reaValObj <- CSCoreClass(treeObject = t1,
                           newickCS = newickCS)
  return(reaValObj)
}
