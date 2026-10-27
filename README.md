# Christophorou_EarliestColony

- **Developed by:** Anneliese Jarman
- **Based at:** Imaging Facility, Babraham Institute
- **Completed on:** 9/9/26
- **Developed for:** Noah Shriever and Maria Christophorou

## What this macro does

This macro opens a whole folder of Incucyte data one by one (it splits filenames on hyphens, not underscores, for parsing and looping) and turns each into a time series. It then finds the earliest timepoint at which any GFP+ cells appear in any colony detected in the final frame that is larger than 20 pixels. Finally, it saves the ROIs, the results, and the stack of masks for the full time series, and closes everything.

## Glossary

- **Blowout:** A rapid increase in the auto-thresholded "foreground" area that occurs when the signal to background contrast drops because of illumination inconsistency (coupled with generally low fluorescence signal).