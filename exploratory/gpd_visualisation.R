# Load required package
library(evd)

# Define multiple sets of parameters (scale and shape) for the GPD
params <- list(
  list(scale = 2, shape = 0.5),
  list(scale = 2, shape = -0.5),
  list(scale = 2, shape = 0.3),
  list(scale = 2, shape = -0.3)
)

# Define the x-axis values over which to compute the density
x_vals <- seq(0.001, 10, length.out = 100)

# Set up the plot
plot(NULL, xlim = range(x_vals), ylim = c(0, 1),
     xlab = "x", ylab = "Density", main = "GPD Density Functions")

# Colors for different lines
colors <- c("blue", "red", "green", "purple")

# Plot each density with specified parameters and add to the legend
for (i in seq_along(params)) {
  scale <- params[[i]]$scale
  shape <- params[[i]]$shape
  density_vals <- dgpd(x_vals, loc = 0, scale = scale, shape = shape)
  
  # Add the density line to the plot
  lines(x_vals, density_vals, col = colors[i], lwd = 2)
}

# Add a legend with the scale and shape parameters
legend("topright", legend = sapply(params, function(p) paste0("scale=", p$scale, ", shape=", p$shape)),
       col = colors, lwd = 2, cex = 0.8)
