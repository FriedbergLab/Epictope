
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
output_filename <- gsub(".csv", ".tiff", filename)
score_data <- read.csv(filename)

# Mutate id and min_val to numeric
score_data$id <- as.numeric(score_data$position)
score_data$min_val <- as.numeric(score_data$min)

# Define moving average function
moving_average <- function(x, window_size) {
  n <- length(x)
  ma <- numeric(n - window_size + 1)
  pad <- rep(NA, floor(window_size/2))
  for (i in 1:(n - window_size + 1)) {
    ma[i] <- mean(x[i:(i + window_size - 1)])
  }
  res <- c(pad, ma, pad)
  return(res)
}

# reformat
plot_data <- score_data
plot_data$id <- as.numeric(plot_data$id)
plot_data <- plot_data[order(plot_data$id),]
plot_data$min_val_ma <- moving_average(plot_data$min_val, 7)
plot_data <- plot_data[, c("id", "min_val_ma")]

# Melt data to long format
plot_data <- reshape(plot_data,
                   varying = list(names(plot_data)[names(plot_data) != "id"]),
                   v.names = "value",
                   idvar = "id",
                   direction = "long")
# Convert value to numeric
plot_data$value <- as.numeric(plot_data$value)

# Calculate breaks for y-axis based on y_range
y_range <- round(range(as.numeric(plot_data$value), na.rm = TRUE), 2)
y_breaks <- seq(from = y_range[1], to = y_range[2], by = (y_range[2] - y_range[1])/5)
y_labels <- round(y_breaks, 2)
max_length <- max(nchar(y_labels))
mar_left <- max(4, max_length*0.2) 
min_id <- min(plot_data$id)
max_id <- max(plot_data$id)
num_ticks <- 6 

# function to generate x axis labels in increments of 100
generate_x_labels <- function(min_id, max_id) {
  num_steps <- (max_id - min_id) %/% 100
  result_vector <- seq(0, by = 100, length.out = num_steps + 1)[-1]
  if (tail(result_vector, 1) < max_id) {
    result_vector <- c(result_vector, max_id)
  }
  return(c(1, result_vector))
}
x_ticks <- generate_x_labels(min_id, max_id)
x_labels <- as.character(x_ticks)
par(mar=c(1,1,1,1))

# Open the TIFF device
tiff(paste0("outputs/", output_filename), width = 6, height = 4, units = "in", res = 300)
plot(plot_data$id, as.numeric(plot_data$value), type = "n", xlab = "Amino Acid Position", ylab = "Score",
     xaxt = 'n', yaxt = 'n', xlim = c(1, max_id), ylim = y_range)
lines(plot_data$id, as.numeric(plot_data$value), lty = "solid", col = "black")
axis(1, at = x_ticks, labels = x_labels)
axis(2, at = y_breaks, labels = round(y_breaks, 2))
box()

# Close the TIFF device
dev.off()