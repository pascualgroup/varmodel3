resolve_algo_params() {
    PARAMS_FILE=$1
    NEW_PARAMS_FILE=$2

    PREV_RESULTS="# previous_results_dir = "
    if [ ! -z ${CFG_PREV_RESULTS+x} ]; then
        # CFG_PREV_RESULTS is set, check for empty
        if [ ! -z ${CFG_PREV_RESULTS} ]; then
            # CFG_PREV_RESULTS is not emtpy
            PREV_RESULTS="previous_results_dir = '$CFG_PREV_RESULTS'"
        fi
    fi

    sed "s|\${PREV_RESULTS_DIR}.*|$PREV_RESULTS,|g" $PARAMS_FILE > $NEW_PARAMS_FILE
}