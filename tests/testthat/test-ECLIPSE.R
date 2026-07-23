data(test_seurat)
data(post_qc_contigs)

files <- system.file("extdata/examples", package = "ECLIPSE")
tcr_paths <- list.files(files,
    pattern = "\\.csv.gz$",
    full.names = TRUE
)

test_seurat_with_tcr <- ECLIPSE(
    contig_file_paths = tcr_paths,
    contig_type = "manual",
    seurat_object = test_seurat,
    batch = "index",
    donor = "donor",
    group = "donor",
    original_barcode = "Cell_Barcode",
    detailedEclipseStats = FALSE,
    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
    allowDoubleChains = "both",
    verbose = FALSE
)

out_contigs <- findThreeChainClones(post_qc_contigs,
    chain = "b"
)

filtered_seurat <- tcrDoubletDetect(test_seurat_with_tcr)





test_that("ECLIPSE accepts only valid inputs", {
    expect_error(
        ECLIPSE(
            contig_folder_paths = list(tcr_paths),
            contig_type = "all",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "contig_folder_paths must be a vector"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = list(tcr_paths),
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "contig_file_paths must be a vector"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = list(),
            contig_type = "all",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Please provide folders to contig_folder_paths that have"
    )

    expect_error(
        ECLIPSE(
            contig_type= "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "No contig file paths provided"
    )


    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Batch column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "fake",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Donor column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "fake",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Group column name provided is not in the Seurat object"
    )


    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "bc",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Original barcode column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "fake",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Incorrect contig type"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = 1,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "TRUE or FALSE not provided to detailedEclipseStats"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            format = "fake"
        ),
        "None, none, or Blank not provided to format"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(1, 2, 3, 4),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Incorrect input to combineTcrArgs"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "fake",
            verbose = FALSE
        ),
        "Alpha, Beta, Both, or Neither not provided to allowDoubleChains"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = "fale"
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = 1
        ),
        "TRUE or FALSE not provided to returnFilteredContigs"
    )
})




test_that("ECLIPSE returns a correctly formatted Seurat object", {
    ### Expected columns present in meta.data
    expect_equal(
        c(
            "CTaa", "scR_CTaa", "clonalFrequency",
            "orig.ident", "nCount_RNA"
        ) %in%
            colnames(test_seurat_with_tcr@meta.data),
        rep(TRUE, 5)
    )

    ### Meta data retained
    expect_true(nrow(test_seurat_with_tcr@meta.data) > 1 &
        ncol(test_seurat_with_tcr@meta.data) > 1)

    ### clonalProportion adds up to 1 for each group
    sum_freq_per_donor <- test_seurat_with_tcr@meta.data %>%
        filter(!is.na(clonalProportion)) %>%
        group_by(donor, CTaa) %>%
        summarise(prop = first(clonalProportion)) %>%
        group_by(donor) %>%
        summarise(total = sum(prop))

    expect_equal(mean(sum_freq_per_donor$total), 1)

    ### clonalFrequency is whole numbers
    expect_equal(
        sum(!is.na(test_seurat_with_tcr@meta.data$clonalFrequency) &
              test_seurat_with_tcr$clonalFrequency %% 1 != 0),
        0
    )

    ### Clone names are formatted correctly (i.e no -)
    expect_equal(
        length(grep("-", test_seurat_with_tcr$CTaa)),
        0
    )

    ### Columns didn't get removed
    expect_equal(
        sum(as.character(seq(1, nrow(test_seurat_with_tcr@meta.data)))
        == rownames(test_seurat_with_tcr@meta.data)),
        0
    )

    ### Number of cells annotated in scR_CTaa = CTaa
    expect_equal(
        sum(is.na(test_seurat_with_tcr$scR_CTaa)),
        sum(is.na(test_seurat_with_tcr$CTaa))
    )
})




test_that("findThreeChainClones accepts only valid inputs", {
    expect_error(
        findThreeChainClones(post_qc_contigs[, 30:37],
            chain = "a",
            mode = "standard",
            printStats = TRUE
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "fake",
            mode = "standard",
            printStats = TRUE
        ),
        "Chain should be alpha or beta"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "fake",
            printStats = TRUE
        ),
        "Mode must be standard or show"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "standard",
            printStats = "fake"
        ),
        "printStats must be TRUE or FALSE"
    )
})




test_that("findThreeChainClones returns contigs", {
    ### Returns data.frame
    expect_true(is.data.frame(out_contigs))

    ### This should keep contigs as is or join them together
    ### leading to the same or less rows
    expect_lte(
        nrow(out_contigs),
        nrow(post_qc_contigs)
    )

    ### findThreeChainClones shouldn't add any new columns
    expect_equal(
        ncol(out_contigs),
        ncol(post_qc_contigs)
    )

    ### Expected columns present
    expect_equal(
        sum(c("barcode", "high_confidence", "chain", "cdr3", "umis")
        %in% colnames(out_contigs)),
        5
    )
})




test_that("tcrDoubletDetect accepts only Seurat objects
          after ECLIPSE was ran", {
    expect_error(
        tcrDoubletDetect(test_seurat),
        "The seurat object must have had ECLIPSE ran on it"
    )

    expect_error(
        tcrDoubletDetect(test_seurat_with_tcr,
            filter = "fake",
            singleChainLimit = 2,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "filter, showPlots, and verbose must be TRUE or FALSE"
    )

    expect_error(
        tcrDoubletDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 2.5,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "singleChainLimit and totalChainLimit both must be"
    )

    expect_error(
        tcrDoubletDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 0,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "singleChainLimit must be at least 2"
    )
})



test_that("tcrDoubletDetect returns a filtered Seurat object", {
    ### No columns are added
    expect_equal(
        ncol(test_seurat_with_tcr@meta.data),
        ncol(filtered_seurat@meta.data)
    )

    ### Rows are the same or decreased
    expect_lte(
        nrow(filtered_seurat@meta.data),
        nrow(test_seurat_with_tcr@meta.data)
    )
})
