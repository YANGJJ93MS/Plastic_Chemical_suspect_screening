# Install required package if not already installed
if (!requireNamespace("networkD3", quietly = TRUE)) {
  install.packages("networkD3")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

library(networkD3)
library(RColorBrewer)

# Create a sample dataset based on the provided structure
filepath = 'D:/UCSF_postdoc_topic/REVEAL_topics/REVEAL_200samples_analysis/Manuscript/classification_classyfire/'
datapath = paste0(filepath,"112_chemicals_classification_for_sankygraph.csv")
data <- read.csv(datapath)

# Create a list of unique nodes
all_nodes <- c(
  data$Function_label, data$Polymer_label, data$Production_label,
  data$RTMSMS_level, data$Superclass
)

node_counts <- table(all_nodes)  # Count occurrences of each node
nodes <- names(node_counts)
nodes_with_counts <- paste0(nodes, " (", node_counts, ")")  # Append counts to node labels

# Map nodes to indices
node_map <- setNames(seq_along(nodes) - 1, nodes)

# Create links (source, target, value)
links <- data.frame(
  source = c(
    node_map[data$Function_label],
    node_map[data$Polymer_label],
    node_map[data$Production_label],
    node_map[data$RTMSMS_level]
  ),
  target = c(
    node_map[data$Polymer_label],
    node_map[data$Production_label],
    node_map[data$RTMSMS_level],
    node_map[data$Superclass]
  ),
  value = 1
)

# Define a custom color palette
palette <- brewer.pal(9, "Pastel1")  # Use a soft color palette
colourScale <- sprintf(
  'd3.scaleOrdinal().range(%s)',
  jsonlite::toJSON(palette, auto_unbox = TRUE)
)

# Create the Sankey diagram with custom colors
sankey <- sankeyNetwork(
  Links = links, Nodes = data.frame(name = nodes_with_counts),
  Source = "source", Target = "target", Value = "value",
  NodeID = "name", units = "chemicals", fontSize = 12, nodeWidth = 10,
  nodePadding = 5,# Adjust padding for better spacing
  colourScale = colourScale  # Apply custom color scale
)

# Display the Sankey diagram
sankey <- htmlwidgets::onRender(
  sankey,
  "
  function(el) {
    d3.select(el).selectAll('.node text')
      .style('font-family', 'Arial')
      .style('font-weight', 'bold');
  }
  "
)

# Display the Sankey diagram
sankey

