# Load required packages
library(tidyverse)
library(ggplot2)
library(gridExtra)

# Read in score files
score_folder <- "outputs/scores"
files = list.files(score_folder, "*.tsv")
score_data <- list()
for( i in 1:length(files)){
	score_data[[i]] = read_tsv(paste(score_folder, files[i], sep = "/"))
}
names(score_data) <- sapply(files, function(.x){strsplit(.x, "\\.")[[1]][1]})

# plot score files
fun_plot <- function(.x, title = ""){
	require(gridExtra)
	require(ggplot2)
	temp <- .x
	p1 <- ggplot(temp, aes(x = id, y = normalized_entropy)) + 
		geom_bar(stat = "identity", aes(fill = "Shannon Entropy")) + 
		ylab("") + xlab("") + ylim(0, 1.5) + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "", breaks = "Shannon Entropy", values = "darkolivegreen4") +
		theme(legend.position="top")

	p2 <- ggplot(temp, aes(x = id, y = struct_no)) + 
		geom_bar(stat = "identity", aes(fill = "Predicted Structure Interference"))  + 
		ylab("") + xlab("") + ylim(0, 1) + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "", breaks = "Predicted Structure Interference", values = "#E69F00") +
		theme(legend.position="top")

	p3 <- ggplot(temp, aes(x = id, y = soluble_sa)) + 
		geom_bar(stat = "identity", aes(fill = "Soluble Surface Area")) + 
		ylab("") + xlab("") + ylim(0, 1) + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "", breaks = "Soluble Surface Area", values = "#56B4E9") +
		theme(legend.position="top")
	
	p4 <- ggplot(temp, aes(x = id, y = inv_anchor2)) +
		geom_bar(stat = "identity", aes(fill = "Protein Binding Site")) + 
		ylab("") + xlab("") + ylim(0, 1) + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "", breaks = "Protein Binding Site", values = "darkorchid4") +
		theme(legend.position="top")

	p5 <- temp %>%
		mutate(smooth_score = zoo::rollmean(score, 11, fill = NA)) %>%
		ggplot(aes(x = id, y = score)) +
		geom_bar(stat = "identity", aes(fill = "Tagging Score")) + 
		geom_line(aes(y = smooth_score)) + 
		ylab("") + xlab("") + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "", breaks = "Tagging Score", values = "gray") +
		theme(legend.position="top")

	p6	<- temp %>%
		select(residue, normalized_entropy, struct_no, soluble_sa, inv_anchor2, score) %>%
		rowid_to_column("id") %>%
		group_by(id, residue) %>%
		pivot_longer(cols = c("normalized_entropy", "struct_no", "soluble_sa", "inv_anchor2", "score"), names_to = "var", values_to = "vals") %>%
		slice(which.min(vals)) %>%
		ungroup() %>%
		mutate(var = factor(var, levels = c("normalized_entropy", "struct_no", "soluble_sa", "inv_anchor2"))) %>%
		mutate(vals = ifelse(vals == 0, 0.05, vals)) %>%
		ggplot(aes(x = id, y = vals, fill = var)) +
		geom_bar(stat = "identity") +
		ylab("") + xlab("") + ylim(0, 1) + 
		scale_x_continuous(breaks=temp$id, labels=temp$residue) +
		scale_fill_manual(name = "Variable", values = c("darkolivegreen4", "#E69F00", "#56B4E9", "darkorchid4")) +
		theme(legend.position="top")

	res <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 6, top = title, bottom = "Higher values are better") 
	return(res)
}

for( i in 1:length(score_data)){
	title = names(score_data)[i]
	a <- fun_plot(score_data[[i]], title = title)
	pdf(paste("figures/", title, ".pdf", sep = ""), height = 12, width = 40)
	plot(a)
	dev.off()
}

# check correlations
all_scores <- bind_rows(score_data, .id = "tf") %>%
	select(c("soluble_sa", "iupred2", "anchor2", "entropy", "struct_no", score)) %>%
	rowwise() %>%
	mutate(inv_anchor2 = 1-anchor2) %>%
	rowwise() %>%
	select(-anchor2) %>%
	filter(!is.na(iupred2)) %>%
	select(inv_anchor2, everything()) 

# calculate kendall correlation
method = "spearman"
colnames(all_scores) <- c("Inverse Anchor2 Score", "Soluble Surface Area", "IUPred2", "Shannon Entropy", "Secondary Structure Score", "Tagging Score")
res <- cor(all_scores, method = method, use = "complete.obs")
ind <- upper.tri(res, diag = FALSE)
res[ind] <- 0


#a <- all_scores$`Shannon Entropy`
#b <- all_scores$`Inverse Anchor2 Score`
#cor(a, b, method = method, use = "complete.obs")
# reshape correlation data to long
res_long <- res %>%
	reshape2::melt() %>%
	filter(value != 0 & Var1 != "Tagging Score") %>%
	arrange(desc(value)) %>%
	rowwise() %>%
	mutate(value2 = signif(value, 2)) %>%
	rowwise()

p1 <- res_long %>%
	ggplot(aes(x = Var2, y = Var1, fill = value)) + 
	geom_tile() +
	geom_text(aes(label=value2)) +
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name=paste(method, "Correlation", sep = "\n")) + xlab("") + ylab("") +
 scale_x_discrete(limits = rev(levels(res_long$Var2))[2:6]) + 
 theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
 labs(title = "Correlation of Tagging Variables") 
pdf(paste("corr", method, ".pdf", sep = ""), height = 6, width = 8)
p1
dev.off()

#####
res_spearman <- cor(all_scores, method = "spearman")
ind_spearman <- upper.tri(res_spearman, diag = TRUE)
res_spearman[ind_spearman] <- 0

# reshape correlation data to long
res_long_spearman <- res_spearman %>%
	reshape2::melt() %>%
	filter(value != 0 & Var1 != "Tagging Score") %>%
	arrange(desc(value)) %>%
	rowwise() %>%
	mutate(value2 = signif(value, 2)) %>%
	rowwise()

p1_spearman <- res_long_spearman %>%
	ggplot(aes(x = Var1, y = Var2, fill = value)) + 
	geom_tile() +
	geom_text(aes(label=value2)) +
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation") + xlab("") + ylab("") +
 theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
 labs(title = "Correlation of Tagging Variables") 
pdf("corr_spearman.pdf", height = 6, width = 8)
p1_spearman
dev.off()

#####
p2 <- all_scores %>%
	reshape2::melt(id = "Shannon Entropy") %>%
	filter(variable != "Secondary Structure Score") %>%
	filter(variable != "Tagging Score") %>%
	ggplot(aes(x = `Shannon Entropy`, y = value)) + 
	geom_point(size = 0.5, alpha = 0.1) +
	xlab("Shannon Entropy") + ylab("") + 
	labs(title = "Scatterplot of Shannon Entropy with Variables") +
	facet_wrap(~variable, scales = "free")
pdf("entropy.pdf", height = 4, width = 12)
p2
dev.off()

p3 <- all_scores %>%
	reshape2::melt(id = "Soluble Surface Area") %>%
	filter(variable != "Secondary Structure Score") %>%
	filter(variable != "Tagging Score") %>%
	ggplot(aes(x = `Soluble Surface Area`, y = value)) + 
	geom_point(size = 0.5, alpha = 0.1) +
	xlab("Soluble Surface Area") + ylab("") + 
	labs(title = "Scatterplot of Soluble Surface Area with Variables") +
	facet_wrap(~variable, scales = "free")
pdf("solubility.pdf", height = 4, width = 12)
p3
dev.off()

p4 <- all_scores %>%
	reshape2::melt(id = "Inverse Anchor2 Score" ) %>%
	filter(variable != "Secondary Structure Score") %>%
	filter(variable != "Tagging Score") %>%
	ggplot(aes(x = `Inverse Anchor2 Score`, y = value)) + 
	geom_point(size = 0.5, alpha = 0.1) +
	xlab("Inverse Anchor2 Score" ) + ylab("") + 
	labs(title = "Scatterplot of Inverse Anchor2 Score with Variables") +
	facet_wrap(~variable, scales = "free")
pdf("inverse_anchor2.pdf", height = 4, width = 12)
p4
dev.off()

p5 <- all_scores %>%
	reshape2::melt(id = "IUPred2") %>%
	filter(variable != "Secondary Structure Score") %>%
	filter(variable != "Tagging Score") %>%
	ggplot(aes(x = `IUPred2`, y = value)) + 
	geom_point(size = 0.5, alpha = 0.1) +
	xlab("IUPred2") + ylab("") + 
	labs(title = "Scatterplot of IUPred2 with Variables") +
	facet_wrap(~variable, scales = "free")
pdf("iupred2.pdf", height = 4, width = 12)
p5
dev.off()

p6 <- all_scores %>%
	reshape2::melt(id = "Secondary Structure Score") %>%
	filter(variable != "Tagging Score") %>%
	mutate(`Secondary Structure Score` = as.factor(`Secondary Structure Score`)) %>%
	ggplot(aes(x = `Secondary Structure Score`, y = value)) + 
	geom_violin() + 
	geom_jitter(size = 0.1, alpha = 0.1) +
	xlab("Secondary Structure Score") + ylab("") + 
	labs(title = "Scatterplot of Secondary Structure Score with Variables") +
	facet_wrap(~variable, scales = "free")
pdf("sec_struct.pdf")
p6
dev.off()

#### 

