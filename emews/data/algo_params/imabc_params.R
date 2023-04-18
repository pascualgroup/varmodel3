algo.params <- list(
    extras = list(priors.path = ""),
    imabc.args = list(
        target_fun = NULL,
        # previous_results_dir = '/lus/eagle/projects/ModelCRC/cisnet_calib/miscan/experiments/imabc_2.4.1.nt',
        N_start = 8000,
        # N_start = 26,
        seed = 12345,
        latinHypercube = TRUE,
        N_centers = 10,
        Center_n = 1000,
        N_post = 1000,
        max_iter = 2,
        N_cov_points = 250,
        sample_inflate = 2,
        verbose = TRUE,
        output_tag = "timestamp",
        improve_method = "percentile"
    )
)
