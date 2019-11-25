(TeX-add-style-hook
 "proj1"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "aa"
    "aa10"
    "graphicx"
    "txfonts"
    "amsmath"
    "subfigure")
   (LaTeX-add-labels
    "eq:W_estimation"
    "fig:line_fit"
    "sec:temp_estimation"
    "eq:log_W_relation"
    "eq:T_exc"
    "fig:temp_estimation_sun"
    "eq:exp_profile_params"
    "sec:k"
    "fig:bad_synth_lines")
   (LaTeX-add-bibliographies
    "spectra"))
 :latex)

