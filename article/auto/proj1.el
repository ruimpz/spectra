(TeX-add-style-hook
 "proj1"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "aa"
    "aa10"
    "graphicx"
    "txfonts"
    "amsmath")
   (LaTeX-add-labels
    "eq:W_estimation"
    "fig:line_fit"
    "eq:log_W_relation"
    "eq:T_exc"
    "fig:temp_estimation_sun"
    "eq:exp_profile_params")
   (LaTeX-add-bibliographies
    "spectra"))
 :latex)

