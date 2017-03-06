## epik.cmd: command line interface to epik package

This package provides commands for running specific analysis from [epik](https://github.com/jokergoo/epik) package 
under command-line mode.

To see commands which are supported:

```
~/project\> Rscript -e "epik.cmd::epik()"
Usage: Rscript -e "epik.cmd::epik()" cmd [options]

Available cmd:

  chromatin_states_transitions
  correlated_regions
  cr_add_fdr
  cr_enriched_heatmap
  cr_genic_stat
  cr_global_vis
  cr_gviz
  cr_reduce
  cr_scatterplot
  cr_sig_heatmap
  diff_meth_in_cgi_and_shore
  diff_meth_in_genomic_features
  general_methylation
  hc_difference

For each cmd, use 'Rscript -e "epik.cmd::epik()" cmd --help' to get help.
```

For each command, `--help` option prints the help message. E.g:

```
~/project\> Rscript -e "epik.cmd::epik()" chromatin_states_transitions --help
Visualize chromatin states transitions in two subgroups

Usage: Rscript -e "epik.cmd::epik()" chromatin_states_transitions [options]

  --config character
    A configuration R script. Check the help page of `load_config()` function to
    find out how to properly set one.

  --min1 numeric
    If there are multiple samples in the group 1, it is possible that a segment has
    more than one states asigned to it. If the recurrency of each state is
    relatively low, it means there is no one dominant state for this segment and it
    should be removed. This argument controls the minimal value for the
    recurrency of states in a given segment. The value is a percent. (0.5 by
    default)

  --min2 numeric
    Same as `min1`, but for samples in group 2. (0.5 by default)

  --methdiff numeric
    The segments for which the methylation difference between two groups is less
    than this value are removed. (0 by default)

  --equalscale
    When there are more than one comparisons, should all chord diagrams have same
    scale? (TRUE by default)

  --help
    Print help message and exit.

  --version
    Print version information and exit.

To successfully run it, `chipseq_hooks$chromHMM()` must be set beforehand in
the configuration file.

The orders of chromatin states in the plot can be set by `CHROMATIN_STATES_ORDER`
(a vector of state names) and the color can be set by `CHROMATIN_STATES_COLOR` (a
vector of colors for which names are state names) in the configuration file.

```

## License

GPL (>= 2)
