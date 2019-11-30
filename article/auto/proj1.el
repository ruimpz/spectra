(TeX-add-style-hook
 "proj1"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "aa"
    "aa10"
    "graphicx"
    "txfonts"
    "amsmath"
    "subfigure"
    "booktabs"
    "hyperref")
   (LaTeX-add-labels
    "eq:W_estimation"
    "fig:line_fit"
    "sec:temp_estimation"
    "eq:log_W_relation"
    "eq:T_exc"
    "fig:temp_estimation_sun"
    "sec:database_lookup"
    "eq:exp_profile_params"
    "eq:vsinI"
    "fig:line_fft"
    "eq:rot_profile"
    "fig:rot_profile"
    "sec:k"
    "fig:bad_synth_lines"
    "sec:continuum_estimation"
    "sec:norm1"
    "fig:additive_parameter"
    "sec:norm2"
    "sec:multiplet_selection"
    "sec:star1"
    "fig:star1_temp_estimate"
    "fig:star1_best_fit"
    "tab:star1_params"
    "sec:star2"
    "fig:star2"
    "tab:star2_params"
    "sec:improvements")
   (LaTeX-add-bibliographies
    "spectra"))
 :latex)

