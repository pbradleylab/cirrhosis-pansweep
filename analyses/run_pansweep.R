# Non-lachnos only
PanSweep::PanSweep_Analysis("cirrhosis_3.json", verbose=TRUE)

PanSweep::PanSweep_Analysis("cirrhosis_2.json")


PanSweep::PanSweep_Analysis("cirrhosis.json")


PanSweep::PanSweep_Analysis("cirrhosis.json", correlation_function = PanSweep::spearman_nonzero_wrapper)
