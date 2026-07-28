#' @title Finds clone with 3 total TCR chains and modifies contigs to address
#' this
#'
#' @description Analyzes chain pairings across a sample to find clones with
#' 2 TCR alpha or 2 TCR beta chains. Then, the contigs are modified so that
#' these cells only have 1 alpha and 1 beta contig, with ::: separating each
#' CDR3 sequence. These contigs are then returned.
#'
#'
#' @param contigs tibble of contigs from cellranger all or
#' filtered contig files.
#' @param chain "alpha" or "beta". This tells it whether clones with two of
#' which chain should be looked for.
#' @param mode "standard" is the default, where modified contigs are returned.
#' If "show" is provided, two tables are returned. The first
#' shows the final set of clones with a double chain of one type, and the
#' second shows how/why some of those clones were called as such. The second
#' is pre-filtering, so some clones may not have been called as having 3 total
#' chains in the end.
#' @param printStats specifies whether notes/stats
#' should be printed. Default is TRUE.
#'
#'
#'
#' @return a tibble of modified contigs,
#' except if mode show is selected.
#'
#' @examples
#' data(post_qc_contigs)
#'
#' modified_contigs <- findThreeChainClones(
#'                         contigs = post_qc_contigs,
#'                         chain = "a")
#'
#' @export
findThreeChainClones <- function(contigs,
                                    chain,
                                    mode = "standard",
                                    printStats = TRUE) {
    if (sum(c("barcode", "high_confidence", "chain", "cdr3", "umis") %in%
        colnames(contigs)) != 5) {
        stop(
            "Contigs should be in the all_contig_annotations or ",
            "filtered_contig_annotations format. ",
            "Please try again."
        )
    }

    if (chain %in% c(
        "alpha", "Alpha", "a", "A", "ALPHA",
        "beta", "Beta", "b", "B", "BETA"
    ) == FALSE) {
        stop(
            "Chain should be alpha or beta. ",
            "Please try again."
        )
    }

    if (mode %in% c("standard", "show") == FALSE) {
        stop(
            "Mode must be standard or show ",
            "Please try again."
        )
    }

    if (!is.logical(printStats)) {
        stop(
            "printStats must be TRUE or FALSE. ",
            "Please try again."
        )
    }

    if (chain %in% c("alpha", "Alpha", "a", "A", "ALPHA")) {
        ### THE FIRST PART HERE IS FINDING CLONES WITH 2 ALPHA CHAINS AND
        ### A SHARED BETA CHAIN

        ### Uses scRepertoire to call clones at baseline
        tcrs <- combineTCR(contigs)[[1]]
        tcrs <- tcrs %>% separate(
            col = "CTaa",
            into = c("alpha", "beta"),
            sep = "_", remove = FALSE
        )
        tcrs <- tcrs %>% filter(.data$beta != "NA" & .data$alpha != "NA")

        ### Makes data frame called combos that has all cells with 3+ chains
        ### and lists the different alphas with each beta
        combos_orig <- tcrs %>%
            filter(str_count(.data$alpha, ";") < 2 &
                !str_detect(.data$beta, ";"))
        combos_full <- combos_orig %>%
            group_by(.data$beta, .data$alpha) %>%
            summarise(count = n(), .groups = "drop") %>%
            ungroup() %>%
            group_by(.data$beta) %>%
            filter(!((.data$count / sum(.data$count)) < 0.1 &
                str_detect(.data$alpha, ";"))) %>%
            ungroup()
        combos <- combos_full %>%
            group_by(.data$beta) %>%
            arrange(desc(.data$count)) %>%
            summarise(
                alphas = paste(.data$alpha, collapse = " "),
                counts = paste(.data$count, collapse = " "),
                total = sum(.data$count),
                chains = n()
            )
        if (nrow(combos) > 0) {
            combos <- combos %>%
                filter(.data$chains >= 3) %>%
                select(all_of(c(
                    "chains",
                    "counts",
                    "total",
                    "beta",
                    "alphas"))
                ) %>%
                separate(.data$alphas,
                    into = paste("A", seq(1, max(combos$chains)), sep = ""),
                    sep = " ", remove = FALSE, fill = "right"
                )
        } else if (nrow(combos) == 0) {
            return(contigs)
        }

        combos <- as.data.frame(combos)
        combos <- combos %>% mutate(combo = NA)

        if (nrow(combos) > 0) {
            ### Scans if there are any clones that have multiple chains
            ### and assigns all those cells to 1 clone
            ### If there is 1+ each of A1A2-B, A1-B, and
            ### A2-B then it is called as a double. Or also 3+ of A1A2-B.
            for (i in seq_len(nrow(combos))) {
                end <- combos[i, 1] + 5
                for (j in 6:end - 1) {
                    if (str_detect(combos[i, j], ";") == FALSE) {
                        c1 <- combos[i, j]
                        c2 <- NA
                        for (k in (j + 1):end) {
                            if (str_detect(combos[i, k], ";") == FALSE) {
                                c2 <- combos[i, k]
                                chains <- as.character(combos[i, 6:end])

                                combo1 <- paste(c1, c2, sep = ";")
                                combo2 <- paste(c2, c1, sep = ";")

                                if (combo1 %in% chains &
                                    is.na(combos[i, "combo"])) {
                                    combos[i, "combo"] <- combo1
                                } else if (combo2 %in% chains &
                                    is.na(combos[i, "combo"])) {
                                    combos[i, "combo"] <- combo2
                                }
                            }
                        }
                    }
                }
            }
        }
        if (mode == "show") {
            combo_table <- combos
        }
        combos <- combos %>% select(all_of(c("beta", "combo")))
        combos <- combos %>% filter(!is.na(.data$combo))

        ### Uses combos_full to find more clones that may not have all 3
        ### conditions met,
        ### but still have cells with 2 alpha that are >= 10% and
        ### also 3+ in count
        combos2 <- combos_full %>%
            filter(str_detect(.data$alpha, ";") & .data$count >= 3)
        combos2 <- combos2 %>%
            arrange(desc(.data$count)) %>%
            group_by(.data$beta) %>%
            filter(row_number() == 1) %>%
            ungroup()
        combos2 <- combos2 %>%
            select(all_of("beta"), "combo" = all_of("alpha")) %>%
            filter(!.data$beta %in% combos$beta)

        combos <- rbind(combos, combos2)


        ### This finds the beta chain that is most common
        ### for each alpha chain.
        ### Note: this is not just 3 TCR clones, this is for all clones
        ### We want to be able to filter so that for each 2A 1B clone,
        ### each A must be most commonly seen with that B.
        pairings <- combos_orig %>%
            group_by(.data$alpha, .data$beta) %>%
            summarise(count = n(), .groups = "drop")
        pairings <- pairings %>% separate(.data$alpha,
            into = c("a1", "a2"),
            sep = ";", fill = "right"
        )
        pairings <- pairings %>%
            pivot_longer(
                cols = all_of(c("a1", "a2")),
                names_to = "original", values_to = "alpha"
            ) %>%
            select(-all_of("original")) %>%
            filter(!is.na(.data$alpha))
        pairings <- pairings %>%
            group_by(.data$alpha, .data$beta) %>%
            summarise(
                Count = sum(.data$count),
                .groups = "drop"
            )
        pairings <- pairings %>%
            arrange(desc(.data$Count)) %>%
            group_by(.data$alpha) %>%
            mutate(rank = min_rank(desc(.data$Count))) %>%
            group_by(.data$alpha, .data$rank) %>%
            filter(sum(.data$rank) == 1) %>%
            ungroup() %>%
            select(-all_of("rank"))

        ### This part compares the most common beta for each alpha and
        ### sees if for every 3 TCR clone both alphas are
        ### most commonly with the beta
        ### If this isn't the case, they are filtered out here
        pairings2 <- pairings %>% left_join(combos, by = "beta")
        pairings2 <- pairings2 %>%
            filter(str_detect(
                .data$combo,
                paste("^", .data$alpha, ";|;", .data$alpha, "$", sep = "")
            ))
        pairings2 <- pairings2 %>%
            group_by(.data$beta) %>%
            mutate(combos = n_distinct(.data$combo)) %>%
            group_by(.data$beta, .data$combo) %>%
            summarise(
                count = n(), combos = first(.data$combos),
                .groups = "drop"
            ) %>%
            filter(.data$count == 2)

        if (nrow(pairings2) > 0) {
            if (max(pairings2$combos) > 1) {
                warning("Muliple combos with a beta remain")
            }
        }
        combos <- pairings2 %>%
            select(all_of(c("beta", "combo"))) %>%
            ungroup()
        if (mode == "show") {
            combos3 <- combos
        }

        ### Makes a vector special that has all of the
        ### alpha CDR3s of chains that have 2 alphas
        special <- as.character((combos %>%
            filter(!is.na(.data$combo)) %>%
            select(all_of("combo")))[[1]])

        ### Makes a data frame tcr_mod where combos
        ### (beta CDR3 + combo alpha CDR3s) are
        ### joined to tcrs based off of shared beta CDR3
        ### This allows us to know which barcodes need to change.
        ### Then, if the current alpha CDR3 is detected
        ### within the combo alpha CDR3
        ### that has the same beta CDR3 then it will be assigned
        ### as then new alpha CDR3
        ### Lastly, only the barcode and
        ### combo alpha CDR3 are retained and cells
        ### that aren't getting assigned with a combo alpha CDR3 are removed
        tcr_mod <- tcrs %>% left_join(combos, by = "beta")

        tcr_mod <- tcr_mod %>%
            mutate(alpha_mod = case_when(
                str_detect(
                    .data$combo,
                    paste("^", .data$alpha, ";", "|", ";",
                        .data$alpha, "$",
                        sep = ""
                    )
                ) |
                    .data$combo == .data$alpha ~ .data$combo,
                TRUE ~ NA
            ))

        tcr_mod <- tcr_mod %>%
            select(all_of(c("barcode", "alpha_mod"))) %>%
            filter(!is.na(.data$alpha_mod))

        ### A df called new_contigs is made that takes the
        ### original contigs data and left joins in tcr_mod by shared barcode
        ### Then, if the row is an alpha chain and
        ### has a new alpha combo it is called as being a combo
        new_contigs <- contigs %>% left_join(tcr_mod, by = "barcode")
        new_contigs <- new_contigs %>%
            mutate(new_cdr3 = case_when(
                !is.na(.data$alpha_mod) & .data$chain == "TRA" ~
                    .data$alpha_mod
            ))

        ### new_contigs is then filtered by reads/UMI and then,
        ### for newly assigned clones only,
        ### only the chain with the most amount of UMIs/reads for a
        ### given chain/barcode pair is retained.
        ### This filtering is necessary because essentially each row = 1 chain,
        ### so if you assign 1 row to a combo (2 chains),
        ### then you would be duplicating the data.
        ### Removing duplicates rows ensures that downstream the
        ### clone will not get duplicated.
        new_contigs <- new_contigs %>%
            group_by(.data$barcode, .data$chain) %>%
            arrange(desc(.data$umis), desc(.data$reads)) %>%
            mutate(row_num = row_number()) %>%
            ungroup()
        new_contigs <- new_contigs %>%
            filter(!(!is.na(.data$new_cdr3) & .data$row_num > 1))
        new_contigs <- new_contigs %>%
            mutate(cdr3 = case_when(
                !is.na(.data$new_cdr3) ~ .data$new_cdr3,
                is.na(.data$new_cdr3) ~ .data$cdr3
            ))

        ### Makes a new df ab with the modified contigs.
        ### Feeds the first 32 columns of ab into scRepertoire so
        ### show find new clone calls.
        ### This makes a new df called post
        ab <- new_contigs

        post <- combineTCR(ab)[[1]]

        if (mode == "show") {

            combo_table <- combo_table %>%
                select(-all_of(c("chains", "alphas"))) %>%
                rename("cell_counts_per_alpha_pairing" = all_of("counts"),
                    "total_cells" = all_of("total"),
                    "predicted_double_alpha_pairing" = all_of("combo"),
                    "beta_chain" = all_of("beta"))

            combo_table <- combo_table %>%
                select("total_cells",
                    "cell_counts_per_alpha_pairing",
                    "beta_chain",
                    "predicted_double_alpha_pairing",
                    everything()) %>%
                arrange(desc(.data$total_cells))

            colnames(combo_table)[5:ncol(combo_table)] <-
                str_replace(colnames(combo_table)[5:ncol(combo_table)],
                    "A", "observed_alpha_pairing")


            message("Returning a list. The first table shows all ",
                    "beta chains predicted to have a double alpha ",
                    "pairing. The second table shows how/why many ",
                    "of these clones were called as such. Please note,",
                    " the second table is pre-filtering, so not all ",
                    "clones may been in the first table since some may ",
                    "have been filtered out."
                    )
            colnames(combos3) <- c("beta_chain",
                "predicted_double_alpha_pairing")

            return(list(combos3, combo_table))

        }





        ### THIS SECOND PART IS FINDING ALL THE CELLS WITH AN ALPHA AND
        ### A MISSING BETA THAT ARE IN CLONES WITH COMBO ALPHA CHAINS
        ### THEN THESE CELLS ARE ASSIGNED TO THEIR COMBO CLONE SO
        ### THEY ARE CALLED DOWNSTREAM AS 1 CLONE

        if (length(special) > 0) {
            ### Makes a df vdj which summarizes the new clone sizes
            ### Then, if the alpha isn't missing,
            ### calls the number of combo CDR3s it could be from.
            vdj <- post %>%
                group_by(.data$CTaa) %>%
                summarise(count = n()) %>%
                arrange(desc(.data$count)) %>%
                ungroup()
            vdj <- vdj %>%
                separate(.data$CTaa,
                    remove = FALSE,
                    into = c("alpha", "beta"), sep = "_"
                )
            vdj <- vdj %>%
                rowwise() %>%
                mutate(special_count = case_when(
                    .data$alpha != "NA"
                    ~ sum(str_detect(
                            special,
                            paste("^", .data$alpha, ";", "|", ";",
                                .data$alpha, "$",
                                sep = ""
                            )
                        ) | special == .data$alpha),
                    .data$alpha == "NA" ~ 0
                ))


            ### Makes a new column called special_tcr that has the combo alpha
            ### that the alpha chain belongs to.
            vdj$special_tcr <- NA
            for (i in seq_len(nrow(vdj))) {
                if (vdj$special_count[i] == 1) {
                    vdj$special_tcr[i] <-
                        special[which(str_detect(
                            special,
                            paste("^", vdj$alpha[i], ";", "|", ";",
                                vdj$alpha[i], "$",
                                sep = ""
                            )
                        ) | special == vdj$alpha[i])]
                }
            }

            vdj_bad <- vdj %>% filter(.data$special_count > 1 &
                .data$beta == "NA")
            if (nrow(vdj_bad) > 0) {
                message(
                    "There are some alpha chains that are in ",
                    "multiple combo clones. ",
                    "Because of this, ",
                    sum(vdj_bad$count),
                    " cells with these alphas and no beta won't be",
                    "able to be assinged to a combo clone"
                )
            }

            ### Makes a df call special_cells that calls the cells with beta NA
            ### and that belong to combo alpha clones
            ### Begins by only retaining cells that could be called
            ### and then calculates the percentage of cells with each beta
            ### and the same combo alpha.
            ### You can think of these percentages as
            ### the percentage likelihood of
            ### each beta for cells with a NA beta and the same alpha
            special_cells <- vdj %>%
                filter(!is.na(.data$special_tcr)) %>%
                mutate(mod_count = case_when(
                    .data$beta == "NA" ~ 0,
                    .data$beta != "NA" ~ .data$count
                ))
            special_cells <- special_cells %>%
                group_by(.data$special_tcr) %>%
                mutate(group_count = sum(.data$mod_count)) %>%
                ungroup()
            special_cells <- special_cells %>%
                mutate(group_pct = 100 * .data$mod_count / .data$group_count)

            ### Next we try to assign for NA beta cells.
            ### This starts with making a column new_alpha which
            ### is only there for rows where the beta is 80%+ likely.
            ### We then add that same alpha (that pairs with the likely beta)
            ### to a new column new_alpha2 for all cells
            ### in the same alpha group.
            ### Lastly, makes a new column new_alpha3 that,
            ### for cells with no beta and a paired alpha,
            ### has the combo alpha chain CDR3
            special_cells <- special_cells %>%
                mutate(new_alpha = case_when(
                    .data$group_pct >= 80 ~ .data$alpha,
                    .data$group_pct < 80 ~ NA
                ))
            special_cells <- special_cells %>%
                arrange(desc(.data$group_pct)) %>%
                group_by(.data$special_tcr) %>%
                mutate(new_alpha2 = dplyr::first(.data$new_alpha)) %>%
                ungroup()
            special_cells <- special_cells %>%
                rowwise() %>%
                mutate(new_alpha3 = case_when(
                    is.na(.data$new_alpha2)
                    ~ NA,
                    .data$beta == "NA" & (str_detect(
                        .data$new_alpha2,
                        paste("^", .data$alpha, ";", "|", ";",
                            .data$alpha, "$",
                            sep = ""
                        )
                    ) |
                        .data$new_alpha2 == .data$alpha) ~ .data$new_alpha2,
                    TRUE ~ NA
                )) %>%
                ungroup()

            ### Makes a df special_cells2 that selects only the CDR3s
            ### and the new alpha combo CDR3 and then
            ### filters rows with no new assignment
            special_cells2 <- special_cells %>%
                dplyr::select(all_of(c("CTaa", "new_alpha3"))) %>%
                filter(!is.na(.data$new_alpha3))

            ### Joins special_cells with post (the scRepertoire output with
            ### combo alpha chains) by CDR3s to make vdj2
            ### This allow us to know which barcodes need to change.
            ### Barcodes that don't need to change are then removed
            vdj2 <- post %>% left_join(special_cells2, by = "CTaa")
            vdj2 <- vdj2 %>%
                select(all_of(c("barcode", "new_alpha3"))) %>%
                filter(!is.na(.data$new_alpha3))

            ### ab (contigs after joining alpha chains) is
            ### joined to vdj2 to make ab2
            ### ab2 is filtered to retain only the row of that barcode/chain
            ### combo with the most UMIs/reads.
            ### This is only for alpha chains that are changing
            ### This is necessary because, again, not doing so would result
            ### in some duplications down stream
            ### Lastly, the alpha cdr3 for these combo alpha - missing beta
            ### cells are changed
            ab2 <- ab %>% left_join(vdj2, by = "barcode")
            ab2 <- ab2 %>%
                group_by(.data$barcode, .data$chain) %>%
                arrange(desc(.data$umis), desc(.data$reads)) %>%
                mutate(row_num = row_number()) %>%
                ungroup()
            ab2 <- ab2 %>%
                filter(!(!is.na(.data$new_alpha3) & .data$row_num > 1 &
                    .data$chain == "TRA"))
            ab2 <- ab2 %>%
                mutate(cdr3 = case_when(
                    .data$chain == "TRA" & !is.na(.data$new_alpha3) ~
                        .data$new_alpha3,
                    TRUE ~ .data$cdr3
                ))
        } else if (length(special) == 0) {
            ab2 <- ab
        }

        if (printStats == TRUE) {
            message(
                length(special), " clones with an extra TCR\u03B1 chain ",
                "were detected, comprising ",
                round(sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRA") /
                n_distinct(ab2$barcode) * 100, 1), "% of cells (",
                sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRA"),
                " total cells)"
            )
        }

        ab2 <- ab2 %>% select(-all_of(c(
            "alpha_mod",
            "new_cdr3",
            "row_num"))
        )
        return(ab2)

    } else if (chain %in% c("beta", "Beta", "b", "B", "BETA")) {
        ### THE FIRST PART HERE IS FINDING CLONES WITH 2 BETA CHAINS AND
        ### A SHARED ALPHA CHAIN

        ### Uses scRepertoire to call clones at baseline
        tcrs <- combineTCR(contigs)[[1]]
        tcrs <- tcrs %>%
            separate(
                col = "CTaa", into = c("alpha", "beta"),
                sep = "_", remove = FALSE
            )
        tcrs <- tcrs %>% filter(.data$beta != "NA" & .data$alpha != "NA")

        ### Makes data frame called combos that has all cells with 3+ chains
        ### and lists the different betas with each alpha
        combos_orig <- tcrs %>%
            filter(str_count(.data$beta, ";") < 2 &
                !str_detect(.data$alpha, ";"))
        combos_full <- combos_orig %>%
            group_by(.data$beta, .data$alpha) %>%
            summarise(count = n(), .groups = "drop") %>%
            ungroup() %>%
            group_by(.data$alpha) %>%
            filter(!((.data$count / sum(.data$count)) < 0.1 &
                str_detect(.data$beta, ";"))) %>%
            ungroup()
        combos <- combos_full %>%
            group_by(.data$alpha) %>%
            arrange(desc(.data$count)) %>%
            summarise(
                betas = paste(.data$beta, collapse = " "),
                counts = paste(.data$count, collapse = " "),
                total = sum(.data$count),
                chains = n()
            )

        if (nrow(combos) > 0) {
            combos <- combos %>%
                filter(.data$chains >= 3) %>%
                select(all_of(c(
                    "chains",
                    "counts",
                    "total",
                    "alpha",
                    "betas"))
                ) %>%
                separate(.data$betas,
                    into = paste("A", seq(1, max(combos$chains)),
                        sep = ""
                    ),
                    sep = " ", remove = FALSE, fill = "right"
                )
        } else if (nrow(combos) == 0) {
            return(contigs)
        }

        combos <- as.data.frame(combos)
        combos <- combos %>% mutate(combo = NA)


        if (nrow(combos) > 0) {
            ### Scans if there are any clones that have multiple chains and
            ### assigns all those cells to 1 clone
            ### If there is 1+ each of A-B1B2, A-B1, and A-B2 then it is
            ### called as a double. Or also 3+ of A-B1B2.
            for (i in seq_len(nrow(combos))) {
                end <- combos[i, 1] + 5
                for (j in 6:end - 1) {
                    if (str_detect(combos[i, j], ";") == FALSE) {
                        c1 <- combos[i, j]
                        c2 <- NA
                        for (k in (j + 1):end) {
                            if (str_detect(combos[i, k], ";") == FALSE) {
                                c2 <- combos[i, k]
                                chains <- as.character(combos[i, 6:end])

                                combo1 <- paste(c1, c2, sep = ";")
                                combo2 <- paste(c2, c1, sep = ";")

                                if (combo1 %in% chains &
                                    is.na(combos[i, "combo"])) {
                                    combos[i, "combo"] <- combo1
                                } else if (combo2 %in% chains &
                                    is.na(combos[i, "combo"])) {
                                    combos[i, "combo"] <- combo2
                                }
                            }
                        }
                    }
                }
            }
        }
        if (mode == "show") {
            combo_table <- combos
        }
        combos <- combos %>% select(all_of(c("alpha", "combo")))
        combos <- combos %>% filter(!is.na(.data$combo))

        ### Uses combos_full to find more clones that may not
        ### have all 3 conditions met,
        ### but still have cells with 2 beta that are >= 10% and
        ### also 3+ in count
        combos2 <- combos_full %>%
            filter(str_detect(.data$beta, ";") & .data$count >= 3)
        combos2 <- combos2 %>%
            arrange(desc(.data$count)) %>%
            group_by(.data$alpha) %>%
            filter(row_number() == 1) %>%
            ungroup()
        combos2 <- combos2 %>%
            select(all_of("alpha"), "combo" = all_of("beta")) %>%
            filter(!.data$alpha %in% combos$alpha)

        combos <- rbind(combos, combos2)

        #### This finds the alpha chain that is most common for
        ### each beta chain.
        ### Note: this is not just 3 TCR clones, this is for all clones
        #### We want to be able to filter so that for each 2B 1A clone,
        ### each B must be most commonly seen with that A.
        pairings <- combos_orig %>%
            group_by(.data$alpha, .data$beta) %>%
            summarise(count = n(), .groups = "drop")
        pairings <- pairings %>%
            separate(.data$beta,
                into = c("b1", "b2"),
                sep = ";", fill = "right"
            )
        pairings <- pairings %>%
            pivot_longer(
                cols = all_of(c("b1", "b2")),
                names_to = "original", values_to = "beta"
            ) %>%
            select(-all_of("original")) %>%
            filter(!is.na(.data$beta))
        pairings <- pairings %>%
            group_by(.data$alpha, .data$beta) %>%
            summarise(
                Count = sum(.data$count),
                .groups = "drop"
            )
        pairings <- pairings %>%
            arrange(desc(.data$Count)) %>%
            group_by(.data$beta) %>%
            mutate(rank = min_rank(desc(.data$Count))) %>%
            group_by(.data$beta, .data$rank) %>%
            filter(sum(.data$rank) == 1) %>%
            ungroup() %>%
            select(-all_of("rank"))

        ### This part compares the most common alpha for each beta and
        ### sees if for every 3 TCR clone both betas are
        ### most commonly with the alpha
        ### If this isn't the case, they are filtered out here
        pairings2 <- pairings %>% left_join(combos, by = "alpha")
        pairings2 <- pairings2 %>%
            filter(str_detect(
                .data$combo,
                paste("^", .data$beta, ";|;", .data$beta, "$", sep = "")
            ))
        pairings2 <- pairings2 %>%
            group_by(.data$alpha) %>%
            mutate(combos = n_distinct(.data$combo)) %>%
            group_by(.data$alpha, .data$combo) %>%
            summarise(
                count = n(), combos = first(.data$combos),
                .groups = "drop"
            ) %>%
            filter(.data$count == 2)

        if (nrow(pairings2) > 0) {
            if (max(pairings2$combos) > 1) {
                warning("Muliple combos with an alpha remain")
            }
        }
        combos <- pairings2 %>%
            select(all_of(c("alpha", "combo"))) %>%
            ungroup()
        if (mode == "show") {
            combos3 <- combos
        }

        ### Makes a vector special that has all of the beta CDR3s of
        ### chains that have 2 betas
        special <- as.character((combos %>% filter(!is.na(.data$combo)) %>%
            select(all_of("combo")))[[1]])

        ### Makes a data frame tcr_mod where combos
        ### (alpha CDR3 + combo beta CDR3s)
        ### is joined to to tcrs based off of shared alpha CDR3
        ### This allows us to know which barcodes need to change.
        ### Then, if the current beta CDR3 is detected within
        ### the combo beta CDR3 that
        ### has the same alpha CDR3 then it will be assigned as
        ### then new beta CDR3
        ### Lastly, only the barcode and combo beta CDR3 are retained
        ### and cells that aren't getting assigned with a combo
        ### beta CDR3 are removed
        tcr_mod <- tcrs %>% left_join(combos, by = "alpha")
        tcr_mod <- tcr_mod %>%
            mutate(beta_mod = case_when(
                str_detect(
                    .data$combo,
                    paste("^", .data$beta, ";", "|", ";",
                        .data$beta, "$",
                        sep = ""
                    )
                ) |
                    .data$combo == .data$beta ~ .data$combo,
                TRUE ~ NA
            ))

        tcr_mod <- tcr_mod %>%
            select(all_of(c("barcode", "beta_mod"))) %>%
            filter(!is.na(.data$beta_mod))

        ### A df called new_contigs is made that takes
        ### the original contigs data
        ### and left joins in tcr_mod by shared barcode
        ### Then, if the row is a beta chain and has a new beta combo it is
        ### called as being a combo
        new_contigs <- contigs %>% left_join(tcr_mod, by = "barcode")
        new_contigs <- new_contigs %>%
            mutate(new_cdr3 = case_when(
                !is.na(.data$beta_mod) & .data$chain == "TRB" ~ .data$beta_mod
            ))

        ### new_contigs is then filtered by reads/UMI and then,
        ### for newly assigned clones only,
        ### only the chain with the most amount of UMIs/reads for a
        ### given chain/barcode pair is retained.
        ### This filtering is necessary because essentially each row = 1 chain,
        ### so if you assign 1 row to a combo (2 chains),
        ### then you would be duplicating the data.
        ### Removing duplicates rows ensures that downstream
        ### the clone will not get duplicated.
        new_contigs <- new_contigs %>%
            group_by(.data$barcode, .data$chain) %>%
            arrange(desc(.data$umis), desc(.data$reads)) %>%
            mutate(row_num = row_number()) %>%
            ungroup()
        new_contigs <- new_contigs %>%
            filter(!(!is.na(.data$new_cdr3) & .data$row_num > 1))
        new_contigs <- new_contigs %>%
            mutate(cdr3 = case_when(
                !is.na(.data$new_cdr3) ~ .data$new_cdr3,
                is.na(.data$new_cdr3) ~ .data$cdr3
            ))

        ### Makes a new df ab with the modified contigs.
        ### Feeds the first 32 columns of ab into scRepertoire so
        ### show find new clone calls.
        ### This makes a new df called post
        ab <- new_contigs

        post <- combineTCR(ab)[[1]]

        if (mode == "show") {

            combo_table <- combo_table %>%
                select(-all_of(c("chains", "betas"))) %>%
                rename("cell_counts_per_beta_pairing" = all_of("counts"),
                        "total_cells" = all_of("total"),
                        "predicted_double_beta_pairing" = all_of("combo"),
                        "alpha_chain" = all_of("alpha"))

            combo_table <- combo_table %>%
                select("total_cells",
                        "cell_counts_per_beta_pairing",
                        "alpha_chain",
                        "predicted_double_beta_pairing",
                        everything()) %>%
                arrange(desc(.data$total_cells))

            colnames(combo_table)[5:ncol(combo_table)] <-
                str_replace(colnames(combo_table)[5:ncol(combo_table)],
                    "A", "observed_beta_pairing")


            message("Returning a list. The first table shows all ",
                    "alpha chains predicted to have a double beta ",
                    "pairing. The second table shows how/why many ",
                    "of these clones were called as such. Please note,",
                    " the second table is pre-filtering, so not all ",
                    "clones may been in the first table since some may ",
                    "have been filtered out."
                    )
            colnames(combos3) <- c("alpha_chain",
                "predicted_double_beta_pairing")

            return(list(combos3, combo_table))

        }



        ### THIS SECOND PART IS FINDING ALL THE CELLS WITH AN BETA AND
        ### A MISSING ALPHA THAT ARE IN CLONES WITH COMBO BETA CHAINS
        ### THEN THESE CELLS ARE ASSIGNED TO THEIR COMBO CLONE SO
        ### THEY ARE CALLED DOWNSTREAM AS 1 CLONE
        if (length(special) > 0) {
            ### Makes a df vdj which summarizes the new clone sizes
            ### Then, if the beta isn't missing,
            ### calls the number of combo CDR3s it could be from.
            vdj <- post %>%
                group_by(.data$CTaa) %>%
                summarise(count = n()) %>%
                arrange(desc(.data$count)) %>%
                ungroup()
            vdj <- vdj %>% separate(.data$CTaa,
                remove = FALSE,
                into = c("alpha", "beta"), sep = "_"
            )
            vdj <- vdj %>%
                rowwise() %>%
                mutate(special_count = case_when(
                    .data$beta != "NA"
                    ~ sum(str_detect(
                            special,
                            paste("^", .data$beta, ";", "|", ";", .data$beta,
                                "$",
                                sep = ""
                            )
                        ) |
                            special == .data$beta),
                    .data$beta == "NA"
                    ~ 0
                ))

            ### Makes a new column called special_tcr that has
            ### the combo beta that
            ### the beta chain belongs to.
            vdj$special_tcr <- NA
            for (i in seq_len(nrow(vdj))) {
                if (vdj$special_count[i] == 1) {
                    vdj$special_tcr[i] <- special[which(str_detect(
                        special,
                        paste("^", vdj$beta[i], ";", "|", ";",
                            vdj$beta[i], "$",
                            sep = ""
                        )
                    ) |
                        special == vdj$beta[i])]
                }
            }

            vdj_bad <- vdj %>% filter(.data$special_count > 1 &
                .data$alpha == "NA")
            if (nrow(vdj_bad) > 0) {
                message(
                    "There are some beta chains that are ",
                    "in multiple combo clones. Because of this, ",
                    sum(vdj_bad$count),
                    " cells with these betas and no alpha won't be able ",
                    "to be assinged to a combo clone"
                )
            }

            ### Makes a df call special_cells that calls
            ### the cells with alpha NA
            ### and that belong to combo beta clones
            ### Begins by only retaining cells that could be called and
            ### then calculates the percentage of cells with each
            ### alpha and the same combo beta
            ### You can think of these percentages as
            ### the percentage likelihood of
            ### each alpha for cells with a NA alpha and the same beta
            special_cells <- vdj %>%
                filter(!is.na(.data$special_tcr)) %>%
                mutate(mod_count = case_when(
                    .data$alpha == "NA" ~ 0,
                    .data$alpha != "NA" ~ .data$count
                ))
            special_cells <- special_cells %>%
                group_by(.data$special_tcr) %>%
                mutate(group_count = sum(.data$mod_count)) %>%
                ungroup()
            special_cells <- special_cells %>%
                mutate(group_pct = 100 * .data$mod_count / .data$group_count)

            ### Next we try to assign for NA alpha cells.
            ### This starts with making a column new_beta which is
            ### only there for rows where the alpha is 80%+ likely.
            ### We then add that same beta (that pairs with the likely alpha)
            ### to a new column new_beta2 for all cells in the same beta group.
            ### Lastly, makes a new column new_beta3 that,
            ### for cells with no alpha and a paired beta,
            ### has the combo beta chain CDR3
            special_cells <- special_cells %>%
                mutate(new_beta = case_when(
                    .data$group_pct >= 80 ~ .data$beta,
                    .data$group_pct < 80 ~ NA
                ))
            special_cells <- special_cells %>%
                arrange(desc(.data$group_pct)) %>%
                group_by(.data$special_tcr) %>%
                mutate(new_beta2 = dplyr::first(.data$new_beta)) %>%
                ungroup()
            special_cells <- special_cells %>%
                rowwise() %>%
                mutate(new_beta3 = case_when(
                    is.na(.data$new_beta2) ~ NA,
                    .data$alpha == "NA" & (str_detect(
                        .data$new_beta2,
                        paste("^", .data$beta, ";", "|", ";",
                            .data$beta, "$",
                            sep = ""
                        )
                    ) |
                        .data$new_beta2 == .data$beta) ~ .data$new_beta2,
                    TRUE ~ NA
                )) %>%
                ungroup()

            ### Makes a df special_cells2 that selects only the CDR3s and
            ### the new beta combo CDR3 and then filters rows
            ### with no new assignment
            special_cells2 <- special_cells %>%
                dplyr::select(all_of(c("CTaa", "new_beta3"))) %>%
                filter(!is.na(.data$new_beta3))

            ### Joins special_cells with post (the scRepertoire output
            ### with combo beta chains) by CDR3s to make vdj2
            ### This allow us to know which barcodes need to change.
            ### Barcodes that don't need to change are then removed
            vdj2 <- post %>% left_join(special_cells2, by = "CTaa")
            vdj2 <- vdj2 %>%
                select(all_of(c("barcode", "new_beta3"))) %>%
                filter(!is.na(.data$new_beta3))

            ### ab (contigs after joining beta chains) is
            ### joined to vdj2 to make ab2
            ### ab2 is filtered to retain only the row of that barcode/chain
            ### combo with the most UMIs/reads.
            ### This is only for beta chains that are changing
            ### This is necessary because, again, not doing so would
            ### result in some duplications down stream
            ### Lastly, the beta cdr3 for these combo beta -
            ### missing alpha cells are changed
            ab2 <- ab %>% left_join(vdj2, by = "barcode")
            ab2 <- ab2 %>%
                group_by(.data$barcode, .data$chain) %>%
                arrange(desc(.data$umis), desc(.data$reads)) %>%
                mutate(row_num = row_number()) %>%
                ungroup()
            ab2 <- ab2 %>%
                filter(!(!is.na(.data$new_beta3) & .data$row_num > 1 &
                    .data$chain == "TRB"))
            ab2 <- ab2 %>%
                mutate(cdr3 = case_when(
                    .data$chain == "TRB" & !is.na(.data$new_beta3) ~
                        .data$new_beta3,
                    TRUE ~ .data$cdr3
                ))
        } else if (length(special) == 0) {
            ab2 <- ab
        }

        if (printStats == TRUE) {
            message(
                length(special), " clones with an extra TCR\u03B2 chain ",
                "were detected, comprising ",
                round(sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRB") /
                n_distinct(ab2$barcode) * 100, 1), "% of cells (",
                sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRB"),
                " total cells)"
            )
        }

        ab2 <- ab2 %>% select(-all_of(c(
            "beta_mod",
            "new_cdr3",
            "row_num"))
        )

        return(ab2)
    }
}
