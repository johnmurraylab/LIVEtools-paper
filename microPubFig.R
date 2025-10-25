## Load and Configure Modules -----
library(reticulate)
#assumes that one lready have a python venv at this directory
reticulate::use_virtualenv("~/MurrayLab/LIVEtools_venv/")
reticulate::py_run_string("import sys")
library(LIVEtools)

## B -----
# Process a Single Embryo
#load a embryo file
embryo <- readEmbryoTable(file = "20241108_JIM767_25_L5-edit.zip")[[1]]
time = 139
### Draw the raw data into 3d plot -----
figO <- drawEmbLine(
  embryo, 139, cellSize = 4,
  xSize=0.08651, ySize=0.08651, zSize=0.5, #unprocessed data can have different x,y,z scales, specify here
  aligned = F)

saveEmbImg(
  figO[[1]], output_file = "microPubImg/B_original.png",
  viewPoint = list(x=-0,y=0,z=2.5),
  width = 900, height = 600
)
### Align embryo position/rotation and then plot -----
rotatedEmbryo <- totalRePosition(
  embryo, time = 139, indicatorP = "Cxp", indicatorD = "Cxaa", indicatorR = "MSap", indicatorL = "MSpp",
  xSize=0.08651, ySize=0.08651, zSize=0.5) #the xyz dimension is corrected during rotation
figRotated <- drawEmbLine(
  rotatedEmbryo, 139, xSize = 1, ySize = 1, zSize = 1, cellSize = 4,
  aligned = T)

saveEmbImg(
  figRotated[[1]], output_file = "microPubImg/B_rotated.png",
  viewPoint = list(x=-0,y=0,z=2.5),
  width = 900, height = 600
)

### Plot the position of reference lineage to demonstrate how embryo rotation is accomplished -----
figRotated_ref <- drawEmbLine(
  lineages = c("Cxp", "Cxaa", "MSap", "MSpp"),
  rotatedEmbryo, 139, xSize = 1, ySize = 1, zSize = 1, cellSize = 4,
  aligned = T)
MSap <- grepCells(rotatedEmbryo, lineages = "MSap", times = 139)[,c("x","y","z")]
MSpp <- grepCells(rotatedEmbryo, lineages = "MSpp", times = 139)[,c("x","y","z")]
Cxp <- grepCells(rotatedEmbryo, lineages = "Cxp", times = 139)[,c("x","y","z")]
Cxaa <- grepCells(rotatedEmbryo, lineages = "Cxaa", times = 139)[,c("x","y","z")]
color_s <- viridis::viridis_pal(option = "H")(4)
#### Figure out AP axis, find P side with Cxp lineage -----
figRotated_ref1 <- figRotated_ref[[1]] |>
  plotly::add_trace(
    name = "Cxp_avg", x = mean(Cxp$x), y = mean(Cxp$y), z = mean(Cxp$z),type = "scatter3d", mode="markers",
    marker = list(color = color_s[1], size = 10, symbol = 'x'))
saveEmbImg(
  figRotated_ref1, "microPubImg/B_ref_line_DV_view.png",
  viewPoint = list(x=-0,y=0,z=2.5),
  width = 900, height = 600
)
#### Figure out DV and RL axis, find D side with Cxaa, R side MSap, and L side with MSpp lineage -----
figRotated_ref2 <- figRotated_ref[[1]] |>
  plotly::add_trace(
    name = "Cxaa_avg", x = mean(Cxaa$x), y = mean(Cxaa$y), z = mean(Cxaa$z),type = "scatter3d", mode="markers",
    marker = list(color = color_s[2], size = 10, symbol = 'x'))
figRotated_ref2 <- figRotated_ref2 |>
  plotly::add_trace(
    name = "MSap_avg", x = mean(MSap$x), y = mean(MSap$y), z = mean(MSap$z),type = "scatter3d", mode="markers",
    marker = list(color = color_s[3], size = 10, symbol = 'x'))
figRotated_ref2 <- figRotated_ref2 |>
  plotly::add_trace(
    name = "MSpp_avg", x = mean(MSpp$x), y = mean(MSpp$y), z = mean(MSpp$z),type = "scatter3d", mode="markers",
    marker = list(color = color_s[4], size = 10, symbol = 'x'))
saveEmbImg(
  figRotated_ref2, "microPubImg/B_ref_line_AP_view.png",
  viewPoint = list(x=-2.5,y=0,z=0), up = list(x=0,y=0,z=1),
  width = 900, height = 600
)

## C, the expression value 3D plot -----
figExp <- drawEmbVal(
  rotatedEmbryo, time = 139, aligned = T,
  cellSize = 4, viewPoint = list(x=0,y=0,z=2.5)
)

saveEmbImg(
  figExp[[1]], "microPubImg/C.png",
  viewPoint = list(x=0,y=0,z=2.5),
  width = 900, height = 600
)

## F -----
### batch process -----
# Plot the reporter expression over time for a particular lineage or specified cells for all embryos.
embryos <- CD_In(directory = "JIM767_25/", prefix = "CD",
                 time.prefix = "TIME", time.suffix = ".csv",
                 info.prefix = "", info.suffix = "AuxInfo.csv")
### implement depth correction-----
model <- depthCorrectionParm( #learn model
  embryos, lineage = c("ABara", "ABalp", "E"),
  alignCell = "E", startT = 30, endT = 150, exc_zMin = 0.2, exc_zMax = 1, zMax = 67
)
embryos <- embryos |> #correct with model
  dataCorrection(lineage=c("ABara", "ABalp", "P1"), model = model, zMax = 67, exc_zMin = 0.2, exc_zMax = 1)
yrange = list(0,16000)
xrange = list(-50,200)
### palette for each embryo
embNames <- names(embryos[["CD"]])
library(viridis)
emb_palette <- turbo(length(embNames)) |> setNames(embNames)
### plot average expression at each timepoint within each embryo for selected lineages -----
lineExp <- c("ABalp", "ABara", "ABplpappp", "ABprpappp")
fig_blot_lines <- embryos|>plotBlotLine(
  title = 'average expression of ABalp, ABara, ABplpappp, ABprpappp lineages',
  aligningCell = lineExp[1], align_t = 0, lineages = lineExp,
  color_palette = emb_palette,
  xrange = xrange, yrange = yrange)

plotly::save_image(
  fig_blot_lines, "microPubImg/F.png",
  width = 1200, height = 800
)

## D, select an embryo and plot its lineage tree for MS lineage -----
embryo_sample = embryos[["CD"]][["20241108_JIM767_25_L5"]]
MS_tree <- CD_tree_plot(
  CD = embryo_sample,
  root = "MS", end_time = 240,
  min_gain = 0, max_gain = 5000
)
ggplot2::ggsave(
  filename = "microPubImg/D.png", plot = MS_tree,
  dpi = 300, width = 6, height = 6
)

## E, aligned lineage plot comparing reference vs experiment (RNAi) embryo sample -----
### example of lineage: lit-1i
sys1_lit1i <- readEmbryoTable(file = "CD20140203_sys-1_lit-1i_L1.csv")[[1]]
#fig_lit1i <- drawEmbLine(
#  sys1_lit1i, 139, cellSize = 4,
#  xSize=0.087, ySize=0.087, zSize=0.504, #unprocessed data can have different x,y,z scales, specify here
#  aligned = F)
sys1_lit1i_rotate <- totalRePosition(
  sys1_lit1i, time = 139,
  indicatorP = "Cxp", indicatorD = "Cxaa", indicatorR = "MSap", indicatorL = "MSpp",
  xSize=0.087, ySize=0.087, zSize=0.54) #the xyz dimension is corrected during rotation
fig_lit1i_rotate <- drawEmbLine(
  sys1_lit1i_rotate, 139,
  lineages = c("MS","E","C","D","P4"),
  xSize = 1, ySize = 1, zSize = 1, cellSize = 4,
  aligned = T)

saveEmbImg(
  fig_lit1i_rotate[[1]], "microPubImg/E_lit1i.png",
  viewPoint = list(x=-0,y=0,z=2.5),
  width = 900, height = 600
)

### reference embryo -----
reference <- readEmbryoTable(file = "CDReferenceModel.csv")[[1]]

reference_rotate <- totalRePosition(
  reference, time = 300,
  indicatorP = "Cxp", indicatorD = "Cxaa", indicatorR = "MSap", indicatorL = "MSpp",
  xSize=0.087, ySize=0.087, zSize=0.54) #the xyz dimension is corrected during rotation

fig_reference_rotate <- drawEmbLine(
  reference_rotate, 250, lineages = c("MS","E","C","D","P4"), xSize = 1, ySize = 1, zSize = 1, cellSize = 4,
  aligned = T)

saveEmbImg(
  fig_reference_rotate[[1]], "microPubImg/E_reference.png",
  viewPoint = list(x=-0,y=0,z=2.5),
  width = 900, height = 600
)
