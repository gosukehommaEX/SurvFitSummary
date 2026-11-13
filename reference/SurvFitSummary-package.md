# SurvFitSummary: Survival Extrapolation with Standard Parametric and Spline Models

Prepares time-to-event data with
[`processing_dataset`](https://gosukehommaEX.github.io/SurvFitSummary/reference/processing_dataset.md),
summarises them with Kaplan-Meier estimates
([`summarise_data`](https://gosukehommaEX.github.io/SurvFitSummary/reference/summarise_data.md),
[`km_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/km_table.md)),
Cox proportional hazards models
([`estimate_hazard_ratio`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_hazard_ratio.md)),
the Grambsch-Therneau test
([`test_proportional_hazards`](https://gosukehommaEX.github.io/SurvFitSummary/reference/test_proportional_hazards.md))
and smoothed empirical hazards
([`estimate_empirical_hazard`](https://gosukehommaEX.github.io/SurvFitSummary/reference/estimate_empirical_hazard.md)).
Test data can be simulated with
[`generate_dummy_survival_data`](https://gosukehommaEX.github.io/SurvFitSummary/reference/generate_dummy_survival_data.md).

Standard parametric and Royston-Parmar spline models are fitted with
[`fit_survival_models`](https://gosukehommaEX.github.io/SurvFitSummary/reference/fit_survival_models.md),
and their results are reported by
[`parameter_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/parameter_table.md),
[`covariance_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/covariance_table.md),
[`gof_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/gof_table.md),
[`predict_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/predict_survival.md),
[`survival_summary`](https://gosukehommaEX.github.io/SurvFitSummary/reference/survival_summary.md)
and
[`spline_knot_table`](https://gosukehommaEX.github.io/SurvFitSummary/reference/spline_knot_table.md).

The plots are
[`plot_km`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_km.md),
[`plot_empirical_hazard`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_empirical_hazard.md),
[`plot_survival_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_survival_fits.md),
[`plot_hazard_fits`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_hazard_fits.md),
[`plot_distribution_diagnostics`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_distribution_diagnostics.md),
[`plot_schoenfeld`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_schoenfeld.md)
and
[`plot_qq_survival`](https://gosukehommaEX.github.io/SurvFitSummary/reference/plot_qq_survival.md).

The results are written to an Excel workbook by
[`write_survival_workbook`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_workbook.md)
and to a Word report by
[`write_survival_report`](https://gosukehommaEX.github.io/SurvFitSummary/reference/write_survival_report.md).
The vignettes describe the analysis
([`vignette("SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/SurvFitSummary.md)),
the models
([`vignette("models", package = "SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/models.md))
and the validation of the package
([`vignette("validation", package = "SurvFitSummary")`](https://gosukehommaEX.github.io/SurvFitSummary/articles/validation.md)).

## See also

Useful links:

- <https://gosukehommaEX.github.io/SurvFitSummary/>

- <https://github.com/gosukehommaEX/SurvFitSummary>

- Report bugs at
  <https://github.com/gosukehommaEX/SurvFitSummary/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>
([ORCID](https://orcid.org/0000-0002-6854-5627))

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6854-5627))
