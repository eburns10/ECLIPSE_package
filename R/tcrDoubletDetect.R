#' @title Finds and removes doublets from Seurat objects
#'
#' @description Filters out cells with excessive numbers of TCR chains and
#' helps set thresholds for RNA counts based on TCR chain count.
#'
#' @param object Seurat object.
#' @param filter logical statement whether cells should be filtered or not.
#' Default is TRUE.
#' @param singleChainLimit number of alpha or beta chains that is allowed
#' before filtering occurs. Must be a integer, with a default of 2.
#' @param totalChainLimit number of total TCR chains that is allowed
#' before filtering occurs. Must be an integer. Default is no threshold.
#' @param showPlots logical statement that determines whether plots
#' are provided that help set additional filtering thresholds based on
#' RNA counts and TCR chain count. Default is TRUE.
#' @param verbose specifies whether notes should be printed. Default is TRUE.
#'
#' @return a Seurat object that is filtered if desired.
#'
#' @examples
#' data(test_seurat_with_tcr)
#'
#' cd8_tcr_filtered <- tcrDoubletDetect(
#'                             test_seurat_with_tcr,
#'                             singleChainLimit = 2)
#'
#' @export
tcrDoubletDetect <- function(object,
                                filter = TRUE,
                                singleChainLimit = 2,
                                totalChainLimit = -1,
                                showPlots = TRUE,
                                verbose = TRUE) {
    if (!"CTaa" %in% colnames(object@meta.data) |
        !"scR_CTaa" %in% colnames(object@meta.data)) {
        stop("The seurat object must have had ECLIPSE ran on it")
    }

    if (!is.logical(filter) | !is.logical(showPlots) | !is.logical(verbose)) {
        stop(
            "filter, showPlots, and verbose must be TRUE or FALSE. ",
            "Please try again."
        )
    }

    if (!is.double(singleChainLimit) | !is.double(totalChainLimit)) {
        stop(
            "singleChainLimit and totalChainLimit both must be ",
            "whole numbers of type double. Please try again."
        )
    }

    if (singleChainLimit %% 1 != 0 | totalChainLimit %% 1 != 0) {
        stop(
            "singleChainLimit and totalChainLimit both must be ",
            "whole numbers of type double. Please try again."
        )
    }

    if (singleChainLimit < 1 | (totalChainLimit < 2 & totalChainLimit != -1)) {
        stop(
            "singleChainLimit must be at least 2, and ",
            "totalChainLimit must be at least 2 or ",
            "-1 (to not filter based on total chain count). ",
            "Please try again."
        )
    }

    ### Makes column in meta data that has the number of alpha/beta
    ### chains both originally and after ECLIPSE.
    object@meta.data <- object@meta.data %>%
        separate(.data$highConfChains,
            remove = FALSE,
            sep = "_", into = c("scr_alpha", "scr_beta")
        )
    object@meta.data <- object@meta.data %>%
        separate(.data$CTaa,
            remove = FALSE, sep = "_",
            into = c("ctaa_alpha", "ctaa_beta")
        )

    object@meta.data <- object@meta.data %>%
        mutate(
            scr_alphaCount = case_when(
                .data$scr_alpha == "NA" | is.na(.data$scr_alpha) ~ 0,
                .data$scr_alpha != "NA" ~ str_count(.data$scr_alpha, ";") + 1
            ),
            scr_betaCount = case_when(
                .data$scr_beta == "NA" | is.na(.data$scr_beta) ~ 0,
                .data$scr_beta != "NA" ~ str_count(.data$scr_beta, ";") + 1
            )
        )


    ### This adjusts for if the clones are called by my pipeline as supposed
    ### to be having 2 of a chain.
    ### Still if something has 2 that doesn't mean it's 100% a doublet
    ### because not everything can be called by my pipeline
    object@meta.data <- object@meta.data %>%
        mutate(
            adj_scr_alphaCount = case_when(
                .data$scr_alpha == str_replace(.data$ctaa_alpha, ":::", ";") &
                    str_detect(.data$ctaa_alpha, ":::") &
                    .data$cloneCallSource == "VDJdive"
                ~ .data$scr_alphaCount - 1,
                TRUE ~ .data$scr_alphaCount
            ),
            adj_scr_betaCount = case_when(
                .data$scr_beta == str_replace(.data$ctaa_beta, ":::", ";") &
                    str_detect(.data$ctaa_beta, ":::") &
                    .data$cloneCallSource == "VDJdive"
                ~ .data$scr_betaCount - 1,
                TRUE ~ .data$scr_betaCount
            )
        )

    ### Anything with more than the desired number of chains is removed
    cell_count <- nrow(object@meta.data)

    if (filter == TRUE) {
        if (verbose == TRUE) {
            message(
                sum(object$scr_alphaCount > singleChainLimit |
                    object$scr_betaCount > singleChainLimit), " cells (",
                round(sum(object$scr_alphaCount > singleChainLimit |
                    object$scr_betaCount > singleChainLimit) /
                    cell_count * 100, 1),
                "% of total) were removed due to having >",
                singleChainLimit, " TCR\u03B1 or TCR\u03B2 chains originally"
            )
        }

        object <- object[, object$scr_alphaCount <= singleChainLimit &
            object$scr_betaCount <= singleChainLimit]

        if (singleChainLimit < 2) {
            warning(
                "Filtering cells with more than 1 \u03B1 or \u03B2 chain is ",
                "not recommended as some clones ",
                "express 2 \u03B1 and/or \u03B2  chains"
            )
        }
    }


    ### Finds the total amount of chains, both raw and adjusted for
    ### clones called as being 2 of a chain in my pipeline
    object@meta.data <- object@meta.data %>%
        mutate(
            scr_count = .data$scr_alphaCount + .data$scr_betaCount,
            adj_scr_count = .data$adj_scr_alphaCount + .data$adj_scr_betaCount
        )

    if (filter == TRUE) {
        if (totalChainLimit >= 0) {
            if (verbose == TRUE) {
                message(
                    sum(object$scr_count > totalChainLimit), " cells (",
                    round(sum(object$scr_count > totalChainLimit) /
                        cell_count * 100, 1),
                    "% of total) were removed due to having >",
                    totalChainLimit, " total chains originally"
                )
            }

            object <- object[, object$scr_count <= totalChainLimit]

        }

        if (totalChainLimit <= 2 & totalChainLimit > -1) {
            warning(
                "Filtering all cells with more than 2 total chains ",
                "is not recommended as some clones express 2 alpha ",
                "and/or beta chains"
            )
        }
    }


    if (showPlots == TRUE) {
        plot3 <- object@meta.data %>%
            group_by(.data$scr_count) %>%
            filter(n() > 1) %>%
            ggplot(aes(
                x = as.character(.data$scr_count),
                y = .data$nCount_RNA
            )) +
            geom_violin(aes(fill = as.character(.data$scr_count)),
                position = position_dodge(), width = 0.9
            ) +
            geom_boxplot(
                fill = "gray80", width = 0.15, linewidth = 0.6,
                outlier.size = -1
            ) +
            theme_classic() +
            xlab(paste("Number of TCR Chains",
                        "\noriginally detected", sep = "")) +
            ylab(paste("Number of RNA UMIs",
                        "\n(nCount_RNA)", sep = "")) +
            ggtitle("RNA Count vs. TCR Chain Count") +
            theme(
                axis.title = element_text(size = 16),
                legend.position = "none",
                plot.title = element_text(size = 16, face = "bold")
            )

        print(plot3)

        plot4 <- object@meta.data %>%
            group_by(.data$adj_scr_count) %>%
            filter(n() > 1) %>%
            ggplot(aes(
                x = as.character(.data$adj_scr_count),
                y = .data$nCount_RNA
            )) +
            geom_violin(aes(fill = as.character(.data$adj_scr_count)),
                position = position_dodge(), width = 0.9
            ) +
            geom_boxplot(
                fill = "gray80", width = 0.15, linewidth = 0.6,
                outlier.size = -1
            ) +
            theme_classic() +
            xlab(paste("Number of TCR Chains\n(adjusted to account for ",
                "clones that\nexpress 2 of one chain type)", sep = "")) +
            ylab(paste("Number of RNA UMIs\n",
                "(nCount_RNA)", sep = "")) +
            ggtitle("RNA Count vs. TCR Chain Count") +
            theme(
                axis.title = element_text(size = 16),
                legend.position = "none",
                plot.title = element_text(size = 16, face = "bold")
            )
        print(plot4)
    }

    object@meta.data <- object@meta.data %>% select(-all_of(c(
        "scr_alpha",
        "scr_beta",
        "ctaa_alpha",
        "ctaa_beta",
        "scr_alphaCount",
        "scr_betaCount",
        "adj_scr_alphaCount",
        "adj_scr_betaCount",
        "scr_count",
        "adj_scr_count"))
    )
    return(object)
}
