load("spatial_data.RData")

mesh <- fmesher::fm_mesh_2d(
  loc = as.matrix(data[, c("lon", "lat")]),
  min.angle = 24,
  max.edge = c(0.1, 0.3),
  cutoff = 0.05
)
spde <- fmesher::fm_fem(mesh)

interpolator_data <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(data[, c("lon", "lat")])
)

interpolator_prediction <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(prediction_locations)
)

pdf("mesh.pdf", width = 4, height = 4)
plot(mesh)
dev.off()

save(mesh, spde, interpolator_data, interpolator_prediction, file = "spatial_helpers.RData")
