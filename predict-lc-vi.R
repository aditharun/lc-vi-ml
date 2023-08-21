library(tidyverse)
library(haven)
library(labelled)
library(broom)
library(caret)
library(gbm)
library(pROC)
library(scales)
library(xgboost)

set.seed(123)


options <- commandArgs(trailingOnly = TRUE)

papervarnames_file <- "lc-vi-paper-var-names.csv"

#matrix called "merged" with rows being patients and columns being symptoms
#one column called status that has the classification
merged <- read_csv(options[1])

#class1varname is the value of the class 1 as a character
class1_variablename <- options[2]

outdir <- options[3]

plots_only <- options[4]

if (!dir.exists(outdir)){
	dir.create(outdir)
}

#k (number) fold ; p (repeat) repeat cv
cv_number <- 5
cv_repeat <- 5

##########



getPerformance <- function(model_seq, merged, class1varname){
	fit <- 
	    model_seq$pred %>% 
	    arrange(rowIndex) %>%
	    group_by(rowIndex) %>%
	    dplyr::summarize(obs=obs[1],
	              prob=mean(!!sym(class1varname))) %>%
	    ungroup() %>% mutate(id = rowIndex) %>%
	    select(id, obs, prob) 

	dat <- fit %>% mutate(class = ifelse(obs == class1varname, 1, 0))

	rocobj <- roc(dat$class, dat$prob, direction = "<", ci=TRUE)

	 df.roc <- data.frame(spec=rocobj$specificities, sens=rocobj$sensitivities)

	  lb <- rocobj$ci[1] 
	  auc <- rocobj$ci[2] 
	  ub <- rocobj$ci[3]

	  lab.text <- paste0("AUC: ", auc %>% signif(2), " (", lb %>% signif(2), " - ", ub %>% signif(2) , ")")

	  ann.text <- data.frame(spec=0.35, sens=0.25, label=lab.text, n=length(rocobj$cases))

	overall.roc <- ggplot(df.roc, aes(spec, sens)) +
	    geom_line(size=1) +
	      scale_x_reverse(expand=c(0, 0.01),
	                      breaks=c(0, 0.25, 0.5, 0.75, 1),
	                      labels=as.character(
	                          c("0", ".25", ".50", ".75", "1.0"))) +
	      scale_y_continuous(expand=c(0, 0.01),
	                         labels=as.character(
	                             c("0", ".25", ".50", ".75", "1.0")))  +
	      theme_classic(base_size=17) +
	      theme(panel.background=element_blank(),
	            panel.grid=element_blank(),
	            legend.position="top",          axis.line.x.bottom=element_line(color="black"),
	            strip.background=element_blank(),
	            aspect.ratio=0.95) +
	      xlab("Specificity") + ylab("Sensitivity") + geom_text(data = ann.text,label = lab.text, size=5)

	list(roc_curve = overall.roc, auc = auc, lb = lb, ub = ub)
}


gbm_stats <- function(merged, response, class1varname, cv.number, cv.repeats){

	tg_gbm <- expand.grid(n.trees=c(50, 100, 150),
	                    interaction.depth=c(2,3,4),
	                    shrinkage=0.1,
	                    n.minobsinnode=10)

	ctrl <- trainControl(method="repeatedcv",
	                       number=cv.number,
	                       repeats=cv.repeats,
	                       verboseIter=TRUE,
	                       savePredictions="final",
	                       classProbs=TRUE,
	                       index=createMultiFolds(merged$status, cv.number, cv.repeats),
	                       summaryFunction=twoClassSummary)

	model_seq_gbm <- caret::train(status ~ .,
	                            data=merged,
	                            method="gbm",
	                            metric="ROC",
	                            tuneGrid=tg_gbm,
	                            trControl=ctrl)

	gbm_performance <- getPerformance(model_seq_gbm, merged, class1varname)

	stats <- data.frame(auc = gbm_performance$auc, lb = gbm_performance$lb, ub = gbm_performance$ub, n.vars = ncol(merged) - 1)

	return(list(stats = stats, model_seq = model_seq_gbm))
}

getvarimp <- function(model_seq){

	df <- varImp(model_seq, scale=FALSE, useModel = TRUE)$importance %>% mutate(var = rownames(.)) %>% as_tibble() %>% arrange(desc(Overall)) %>% dplyr::mutate(id = 1:n()) %>% mutate(s = Overall / sum(Overall)) %>% filter(s > 0) %>% mutate(s = round(s*100, 2))

	df$var <- gsub("`", "", df$var)

	df
}


getShap <- function(merged, class1varname, model_seq_xgb){

	label.merge <- merged %>% mutate(status = ifelse(status == class1varname, 1, 0)) %>% pull(status)

	fparams <- model_seq_xgb$finalModel$tuneValue %>% as_tibble() 

	ft_cols <- attr(model_seq_xgb$terms, "term.labels") %>% as.character() %>% gsub("`", "", .)

	merged <- merged %>% select(c(ft_cols, "status"))

	bst <- xgboost::xgboost(as.matrix( (merged %>% select(-status))),label.merge , nrounds = fparams$nrounds, eta = fparams$eta, max_depth = fparams$max_depth, subsample = fparams$subsample, gamma = fparams$gamma, colsample_bytree = fparams$colsample_bytree, min_child_weight = fparams$min_child_weight, objective = "binary:logistic", nthread = 2, verbose = 0)

	shap.df <- xgboost::xgb.plot.shap(data = as.matrix((merged %>% mutate(status = ifelse(status == class1varname, 1, 0)) %>% select(-status))), model = bst, top_n = 90, plot=FALSE)$shap_contr %>% as_tibble() %>% mutate_all(~ abs(.)) %>% colMeans()  %>% as.data.frame() 

	shap.df <- shap.df %>% mutate(var = rownames(shap.df)) %>% as_tibble() %>% magrittr::set_colnames(c("mean_abs_shap", "var")) %>% arrange(desc(mean_abs_shap)) %>% mutate(id = 1:n()) 

	return(shap.df)

}


xgb_stats <- function(merged, response, class1varname, cv.number, cv.repeats){

	tg_xgb <- expand.grid(nrounds=c(50, 100, 150),
                    eta=c(0.1, 0.3),
                    gamma = c(0, 5), 
                    max_depth = c(2, 4, 6),
                    colsample_bytree = 1, 
                    min_child_weight = 1, 
                    subsample = 1)

	ctrl <- trainControl(method="repeatedcv",
	                       number=cv.number,
	                       repeats=cv.repeats,
	                       verboseIter=TRUE,
	                       savePredictions="final",
	                       classProbs=TRUE,
	                       index=createMultiFolds(merged$status, cv.number, cv.repeats),
	                       summaryFunction=twoClassSummary)


	model_seq_xgb <- caret::train(status ~ .,
                            data=merged,
                            method="xgbTree",
                            metric="ROC",
                            tuneGrid=tg_xgb,
                            trControl=ctrl)

	xgb_performance <- getPerformance(model_seq_xgb, merged, class1varname)

	stats <- data.frame(auc = xgb_performance$auc, lb = xgb_performance$lb, ub = xgb_performance$ub, n.vars = ncol(merged) - 1)

	return(list(stats = stats, model_seq = model_seq_xgb))
}

if (plots_only != "TRUE"){
	###### GBM PERMUTE ##############

	gbmModel <- gbm_stats(merged, "status", class1_variablename, cv_number, cv_repeat)

	results_gbm <- gbmModel$stats 
	varimp_store_gbm <- list(getvarimp(gbmModel$model_seq))

	startvarlen <- floor(nrow(varimp_store_gbm[[1]])/5)*5 - 5

	while (startvarlen >= 1){

		fts <- c(varimp_store_gbm[[1]] %>% dplyr::slice(1: startvarlen) %>% pull(var), "status")

		fts <- fts %>% gsub("`", "", .)

		gs <- gbm_stats(merged %>% select(fts), "status", class1_variablename, cv_number, cv_repeat)

		results_gbm <- rbind(results_gbm, gs$stats)

		varimp_store_gbm[[length(varimp_store_gbm) + 1]] <- getvarimp(gs$model_seq)

		if (startvarlen <= 20){

			startvarlen <- startvarlen - 1
		} else{
				startvarlen <- startvarlen - 5
		}

	}

	####### XGB GAIN #############



	xgbModel <- xgb_stats(merged, "status", class1_variablename, cv_number, cv_repeat)

	results_xgb_gain <- xgbModel$stats 
	varimp_store_xgb_gain <- list(getvarimp(xgbModel$model_seq)) 

	startvarlen <- floor(nrow(varimp_store_xgb_gain[[1]])/5)*5 - 5

	while (startvarlen >= 1){

		fts <- c(varimp_store_xgb_gain[[1]] %>% dplyr::slice(1: startvarlen) %>% pull(var), "status")

		fts <- fts %>% gsub("`", "", .)

		gs <- xgb_stats(merged %>% select(fts), "status", class1_variablename, cv_number, cv_repeat)

		results_xgb_gain <- rbind(results_xgb_gain, gs$stats)

		if (startvarlen == 1){
			tmp.varimp <- data.frame(Overall = 1, var = fts[fts!="status"], id=1, s = 100)
		} else{
			tmp.varimp <- getvarimp(gs$model_seq)
		}

		varimp_store_xgb_gain[[length(varimp_store_xgb_gain) + 1]] <- tmp.varimp

		if (startvarlen <= 20){

			startvarlen <- startvarlen - 1
		} else{
				startvarlen <- startvarlen - 5
		}

	}


	####### XGB SHAP #############

	xgbModel <- xgb_stats(merged, "status", class1_variablename, cv_number, cv_repeat)

	results_xgb_shap <- xgbModel$stats 
	varimp_store_xgb_shap <- list(getShap(merged, class1_variablename, xgbModel$model_seq))

	startvarlen <- floor(nrow(varimp_store_xgb_shap[[1]])/5)*5 - 5

	while (startvarlen >= 1){

		fts <- c(varimp_store_xgb_shap[[1]] %>% dplyr::slice(1: startvarlen) %>% pull(var), "status")

		fts <- fts %>% gsub("`", "", .)

		gs <- xgb_stats((merged %>% select(fts)), "status", class1_variablename, cv_number, cv_repeat)

		results_xgb_shap <- rbind(results_xgb_shap, gs$stats)

		if (startvarlen == 1){
			tmp.varimp <- data.frame(mean_abs_shap = 1, var = fts[fts!="status"], id=1)
		} else{
			tmp.varimp <- getShap(merged, class1_variablename, gs$model_seq)
		}

		varimp_store_xgb_shap[[length(varimp_store_xgb_shap) + 1]] <- tmp.varimp

		if (startvarlen <= 20){

			startvarlen <- startvarlen - 1
		} else{
				startvarlen <- startvarlen - 5
		}

	}

saveRDS(file=file.path(outdir, "results.rds"), list(results_gbm=results_gbm, results_xgb_shap = results_xgb_shap, results_xgb_gain = results_xgb_gain, varimp_store_gbm = varimp_store_gbm, varimp_store_xgb_shap = varimp_store_xgb_shap, varimp_store_xgb_gain = varimp_store_xgb_gain))

} else{

	results <- file.path(outdir, "results.rds") %>% readRDS(.)
	results_gbm <- results$results_gbm
	results_xgb_shap <- results$results_xgb_shap
	results_xgb_gain <- results$results_xgb_gain
	varimp_store_gbm <- results$varimp_store_gbm
	varimp_store_xgb_gain <- results$varimp_store_xgb_gain
	varimp_store_xgb_shap <- results$varimp_store_xgb_shap

}


true_names <- papervarnames_file %>% read_csv() %>% select(2,3)

varimp_store_gbm <- lapply(varimp_store_gbm, function(x) x %>% left_join(true_names, by=c("var"="new_names")) %>% select(-var) %>% relocate(Overall, proper_names, id, s) %>% magrittr::set_colnames(c("Overall", "var", "id", "s")))

varimp_store_xgb_gain <- lapply(varimp_store_xgb_gain, function(x) x %>% left_join(true_names, by=c("var"="new_names")) %>% select(-var) %>% relocate(Overall, proper_names, id, s) %>% magrittr::set_colnames(c("Overall", "var", "id", "s")))


varimp_store_xgb_shap <- lapply(varimp_store_xgb_shap, function(x) x %>% left_join(true_names, by=c("var"="new_names")) %>% select(-var) %>% relocate(mean_abs_shap, proper_names, id) %>% magrittr::set_colnames(c("mean_abs_shap", "var", "id")))


library(scales)

#PLOTTING
sx_varimp_plot <- function(df, method = "importance"){

	cbb <- ggsci::pal_jama("default")(7)

	if (length(unique(df$colstatus)) <= 1){

		plt <- df %>% ggplot(aes(y=factor(id), x=s), color = "grey80", fill = "grey80") + geom_col() + theme_minimal() + theme(legend.position = "none")


	} else{

		plt <- df %>% ggplot(aes(y=factor(id), x=s, color=colstatus, fill = colstatus)) + geom_col() + theme_minimal() + scale_color_manual(name = "In Final Model", values = c("Yes" = cbb[1], "No" = cbb[2]), labels = c("No" = "No", "Yes" = "Yes")) + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12)) + scale_fill_manual(name = "In Final Model", values = c("Yes" = cbb[1], "No" = cbb[2]), labels = c("No" = "No", "Yes" = "Yes"))

	}

	if (method == "shap"){

			xlabstr <- "Mean Absolute SHAP"

		} else{

			xlabstr <- "Percent Importance (%)"
		}

 	plt	+ scale_y_discrete(breaks = df$id, labels = df$var) + xlab(xlabstr) + theme(panel.grid = element_blank()) + theme(panel.background  = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) +  scale_x_continuous(breaks = pretty_breaks(n = 10)) + ggtitle("") + theme(axis.title = element_text(hjust = 0.5, size = 12)) + theme(axis.text = element_text(size = 10), axis.ticks.x = element_line(color = "grey80"))


}

#relative change: >= 1.5% drop then call it
get_n_vars <- function(df){

	df$diff <- c(0, diff(df$auc))
	df <- df %>% mutate(percent_rel_change = (diff*100) / auc)

	res <- df %>% mutate(idx = 1:n()) %>% filter(round(percent_rel_change, 4) <= -1.5) %>% mutate(nlim = ifelse(n.vars <= 20, 1, 0))

	varn <- res %>% filter(nlim == 1) %>% dplyr::slice(1) %>% pull(idx)

	if (identical(varn, integer(0))){

			varn <- res$idx[1]

	} 

	df$n.vars[(varn - 1)]

}

final.n.vars.gbm <- get_n_vars(results_gbm)

final.var.imp_gbm <- results_gbm %>% mutate(idx = 1:n()) %>% filter(n.vars == final.n.vars.gbm) %>% pull(idx) %>% varimp_store_gbm[[.]]

gbm_roc_n <- results_gbm %>% ggplot(aes(x=n.vars, y=auc)) + geom_point(size = 2.25) + geom_line(size = 0.85) + geom_errorbar(aes(ymin=lb, ymax = ub), size = 0.45, width = 0.1, color="grey80") + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill = "transparent")) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18)) + ylab("AUC w/ 95% CI") + xlab("Number of variables\nin gbm method") + geom_vline(xintercept = final.n.vars.gbm, color="red", linetype = "dashed", size = 0.3) + theme(panel.grid.major.y = element_line(color = "grey80")) + theme(axis.ticks.x = element_line(color = "black"))

gbm_all_sx <- sx_varimp_plot(varimp_store_gbm[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.gbm, "Yes", "No")) )

gbm_final_sx <- sx_varimp_plot(final.var.imp_gbm) 



final.n.vars.gain <- get_n_vars(results_xgb_gain)

final.var.gain_xgb <- results_xgb_gain %>% mutate(idx = 1:n()) %>% filter(n.vars == final.n.vars.gain) %>% pull(idx) %>% varimp_store_xgb_gain[[.]]

xgb_roc_n_gain <- results_xgb_gain %>% ggplot(aes(x=n.vars, y=auc)) + geom_point(size = 2.25) + geom_line(size = 0.85) + geom_errorbar(aes(ymin=lb, ymax = ub), size = 0.45, width = 0.1, color="grey80") + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill = "transparent")) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18)) + ylab("AUC w/ 95% CI") + xlab("Number of variables\nin xgboost gain method") + geom_vline(xintercept = final.n.vars.gain, color="red", linetype = "dashed", size = 0.3) + theme(panel.grid.major.y = element_line(color = "grey80")) + theme(axis.ticks.x = element_line(color = "black"))


xgb_gain_all_sx <- sx_varimp_plot(varimp_store_xgb_gain[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.gain, "Yes", "No")) )

xgb_gain_final_sx <- sx_varimp_plot(final.var.gain_xgb)



final.n.vars.shap <- get_n_vars(results_xgb_shap)

final.var.shap_xgb <- results_xgb_shap %>% mutate(idx = 1:n()) %>% filter(n.vars == final.n.vars.shap) %>% pull(idx) %>% varimp_store_xgb_shap[[.]]

xgb_roc_n_shap <- results_xgb_shap %>% ggplot(aes(x=n.vars, y=auc)) + geom_point(size = 2.25) + geom_line(size = 0.85) + geom_errorbar(aes(ymin=lb, ymax = ub), size = 0.45, width = 0.1, color="grey80") + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill = "transparent")) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18)) + ylab("AUC w/ 95% CI") + xlab("Number of variables\nin xgboost shap method") + geom_vline(xintercept = final.n.vars.shap, color="red", linetype = "dashed", size = 0.3) + theme(panel.grid.major.y = element_line(color = "grey80")) + theme(axis.ticks.x = element_line(color = "black"))


xgb_shap_all_sx <- sx_varimp_plot(varimp_store_xgb_shap[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.shap, "Yes", "No")) %>% magrittr::set_colnames(c("s", "var", "id", "colstatus")), method = "shap")

xgb_shap_final_sx <- sx_varimp_plot(final.var.shap_xgb %>% magrittr::set_colnames(c("s", "var", "id")), method = "shap")






xgb_gain <- varimp_store_xgb_gain[[1]] %>% select(var, s) %>% mutate(type = "xgb gain") %>% arrange(desc(s)) %>% mutate(r = 1:n())
gbm_permute <- varimp_store_gbm[[1]] %>% select(var, s) %>% mutate(type = "gbm permute") %>% arrange(desc(s)) %>% mutate(r = 1:n())
xgb_shap <- varimp_store_xgb_shap[[1]] %>% magrittr::set_colnames(c("s", "var")) %>% select(s, var) %>% relocate(var)  %>% mutate(type = "xgb shap") %>% arrange(desc(s)) %>% mutate(r = 1:n()) 



get_cor_plot <- function(df1, df2){

	x <- left_join(df1, df2, by=c("var" = "var")) %>% filter(complete.cases(.))

	plt <- ggplot(x, aes(x=s.x, y=s.y)) + geom_point(size = 2) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) + xlab(paste0("Variable importance\nfrom ", x$type.x %>% unique(), " model")) + ylab(paste0("Variable importance\nfrom ", x$type.y %>% unique(), " model")) + theme(axis.title = element_text(size = 14), axis.text = element_blank()) + geom_smooth(method = "lm", se = FALSE) + geom_label(x=quantile(x$s.x, .90), y = median(x$s.y), label = paste0("Pearson Corr: ", round(cor(x$s.x, x$s.y, method = "pearson"), 3)), size = 6)

	print(cor.test(x$s.x, x$s.y, method = "pearson"))

	plt

}

a_cor <- get_cor_plot(xgb_gain, gbm_permute)
b_cor <- get_cor_plot(xgb_gain, xgb_shap)
c_cor <- get_cor_plot(xgb_shap, gbm_permute)


hmap <- rbind(varimp_store_gbm[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.gbm, "Yes", "No")) %>% mutate(model = "gbm permute"), 
	varimp_store_xgb_gain[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.gain, "Yes", "No")) %>% mutate(model = "xgb gain"),
	varimp_store_xgb_shap[[1]] %>% mutate(colstatus = ifelse(id <= final.n.vars.shap, "Yes", "No")) %>% magrittr::set_colnames(c("s", "var", "id", "colstatus")) %>% mutate(model = "xgb shap") %>% mutate(Overall = NA)
		)


hmap <- hmap %>% group_by(model) %>% mutate(q = cut(s,quantile(s),include.lowest=TRUE,labels=FALSE) %>% as.character())

#hmap <- hmap %>% mutate(val = ifelse(colstatus == "Yes", r , NA))


hmap <- expand.grid(model=hmap$model %>% unique(), var=hmap$var %>% unique()) %>% as_tibble() %>% left_join(., hmap, by=c("model"="model", "var"="var")) 

hmap$var <- factor(hmap$var, levels = c(final.var.imp_gbm$var, setdiff(hmap$var %>% unique(), final.var.imp_gbm$var)))

t <- hmap %>% mutate(idx = 1:n())

hmap.q <- t %>% group_by(model) %>% filter(colstatus == "Yes") %>% mutate(q = cut(s,quantile(s),include.lowest=TRUE,labels=FALSE) %>% as.character()) %>% ungroup() %>% select(idx, q) %>% right_join(t %>% select(-q), by=c("idx"="idx")) %>% mutate(q = ifelse(is.na(q), "Not Selected" ,q))

hmap.plt <- hmap.q %>% ggplot(aes(x=model, y=var, fill=factor(q, levels = c("Not Selected", "1", "2", "3", "4")))) + geom_tile(color="black") + theme_minimal() + theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent")) + theme(axis.text.y = element_text(size=10), axis.text.x = element_text(size=14), axis.title.x = element_blank(), axis.title.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + scale_fill_manual(name = "Quartile", values = c("Not Selected"="grey95", "1"="grey80", "2"="grey60", "3"="grey40", "4"="grey20")) + ylab("Variables with Non-Zero Importance") 


cor_plot <- cowplot::plot_grid(a_cor + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)), 
	c_cor + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)), 
	b_cor + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)),
	nrow = 3, ncol = 1, labels = c("B", "C", "D"), label_size = 25) 

cowplot::plot_grid(cowplot::plot_grid(hmap.plt, labels = "A", label_size = 25), cor_plot, nrow = 1, rel_widths = c(1.6, 1)) %>% ggsave(plot = ., device = cairo_pdf, filename = file.path(outdir, "supp-model-comparisons.pdf"), units = "in", height = 15, width = 14)



#main figure
cowplot::plot_grid(gbm_roc_n, 
		gbm_final_sx + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) + geom_col(position = position_dodge(width = 0.05), color = "grey40", fill = "grey80") + ylab("Symptom"),
		nrow = 2, labels = c("A", "B"), label_size = 25) %>% ggsave(filename = file.path(outdir, "main-model-fig.pdf"), device = cairo_pdf, units = "in", height = 10, width = 10, plot = .)



#supp figure xgb gain and xgb boost
cowplot::plot_grid(xgb_roc_n_gain, 
				xgb_roc_n_shap, 
				xgb_gain_final_sx + geom_col(position = position_dodge(width = 0.05), color = "grey40", fill = "grey80") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) + ylab("Symptom"), 
				xgb_shap_final_sx + geom_col(position = position_dodge(width = 0.05), color = "grey40", fill = "grey80") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) + ylab("Symptom"), 
				nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"), label_size = 25) %>% ggsave(filename = file.path(outdir, "supp-xgb-models.pdf"), device = cairo_pdf, units = "in", height = 12, width = 20, plot = .)



















