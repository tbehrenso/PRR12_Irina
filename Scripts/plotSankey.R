
# starts from the gse object from RNAseq_Analysis

# --------------------------------
# With network D3 package
# --------------------------------
library(readxl)
library(networkD3)
library(htmlwidgets)
library(jsonlite)

goterms_sankey <- as.data.frame(read_xlsx('/Users/tbehr/Desktop/SankeyTable.xlsx', sheet = 'GSEA_NEU_UP'))

genes_source <- unlist(sapply(goterms_sankey$geneID, FUN=strsplit, split = '/'))
gene_counts <- sapply(sapply(goterms_sankey$geneID, FUN=strsplit, split = '/'), length)

goterms_target <- rep(goterms_sankey$Description, gene_counts)


# --- clean inputs ---
goterms_target <- trimws(goterms_target)
genes_source   <- trimws(genes_source)

# --- build links/nodes ---
links <- data.frame(
  source = genes_source,
  target = goterms_target,
  value  = 1,
  stringsAsFactors = FALSE
)

nodes <- data.frame(
  name = unique(c(links$source, links$target)),
  stringsAsFactors = FALSE
)

# each node gets its own group so genes stay uniquely colored by networkD3
nodes$group <- nodes$name

links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1
links$group <- nodes$group[links$IDtarget + 1]  # links inherit target group

# --- user-specified term colours (names must match term strings exactly) ---
term_colors <- c(
  "nuclear division" = "#e85d49",
  "double-strand break repair" = "#5ec2d9"
  # add more terms/colors as needed
)

# optional: check for missing term names
term_names_in_data <- sort(unique(nodes$name[nodes$name %in% goterms_target]))
missing <- setdiff(term_names_in_data, names(term_colors))
if(length(missing) > 0) {
  warning("You did not specify colours for these terms: ", paste(missing, collapse = ", "))
  # you can either supply colours for them or let them remain as networkD3 defaults
}

# --- build sankey (no fragile combined colourScale) ---
p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", NodeID = "name",
  NodeGroup = "group", LinkGroup = "group",
  fontSize = 12, nodeWidth = 30
)

# --- override term node colours and set links to target colours via onRender ---
term_colors_json <- toJSON(as.list(term_colors), auto_unbox = TRUE)

p <- onRender(p, sprintf("
  function(el, x) {
    var term_colors = %s;

    // override rect fills for nodes that are terms (if mapping exists)
    d3.select(el).selectAll('.node rect')
      .style('fill', function(d) {
        return term_colors[d.name] ? term_colors[d.name] : d3.select(this).style('fill');
      });

    // color links by their target node's name (if that target is in term_colors)
    d3.select(el).selectAll('.link')
      .style('stroke', function(d) {
        return term_colors[d.target.name] ? term_colors[d.target.name] : d3.select(this).style('stroke');
      })
      .style('stroke-opacity', 0.45);

    // optional: make node text black for readability
    d3.select(el).selectAll('.node text').style('fill', '#000');
  }
", term_colors_json))

p


# --------------------------------
# With ggplot2 + ggsankey packages
# --------------------------------
library(ggplot2)
library(ggsankey)
library(dplyr)



df_base <- data.frame(
  gene = c('geneA','geneB','geneB','geneC'),
  term = c('term1','term1','term2','term2')
)

df <- make_long(df_base, gene, term)


ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  scale_fill_discrete(drop=FALSE)




ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Sample Text")



# --------------------------------
# With ggalluvial
# --------------------------------
library(ggalluvial)

# Example dataset
df <- data.frame(
  id = 1:8,
  stage1 = c("A", "B", "C", "C", "D", "D", "E", "F"),
  stage2 = c("X", "X", "X", "Z", "Z", "X", "Z", "Z")
)

# Aggregate counts (needed for ggalluvial)
df_counts <- df %>%
  count(stage1, stage2)

# Plot
ggplot(df_counts,
       aes(axis1 = stage1, axis2 = stage2, y = n)) +
  geom_alluvium(aes(fill = stage2), width = 0.25, alpha = 0.8) +   # flows coloured by final node
  geom_stratum(width = 0.25, fill = "grey80", color = "black") +   # node blocks
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  # node labels
  scale_x_discrete(limits = c("Stage 1", "Stage 2"),
                   expand = c(.1, .05)) +
  labs(fill = "Final node") +
  theme_minimal()





