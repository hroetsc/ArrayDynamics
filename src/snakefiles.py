
outfol = expand("results/SIMresults/{lattice}_J{J}_r{r}/",
  lattice = features["lattice"],
  r = features["r"],
  J = features["J"] if lattice=="Kagome" else features["Jsq"])

# ------------------------------ execute for every combination of parameters ------------------------------
rule DynamicMC:
  output:
    SIMresults = expand(join(outfol, "_c{concentration}_met-{methylation}_rep{replicate}.RData",
      concentration = features["concentration"],
      methylation = features["methylation"],
      replicate = features["replicate"])),
    MC = touch(expand(outfol+"DynamicMC_c{concentration}_met-{methylation}_rep{replicate}.done",
      concentration = features["concentration"],
      methylation = features["methylation"],
      replicate = features["replicate"]))
  shell:
    "Rscript src/0_DynamicMC.R --met {methylation} --lattice {lattice} --c {concentration} --J {J} --r0 {r} --rep {replicate}"



# ------------------------------ selected/aggregated replicates ------------------------------
rule plotGrids:
  input:
    SIMresults = expand(outfol+"_c{concentration}_met-{methylation}_rep1.RData",
      concentration = features["concentration"],
      methylation = features["methylation"])
  output:
    grids = touch(expand(outfol+"GIFs/gifs_c{concentration}_met-{methylation}.done",
      concentration = features["concentration"],
      methylation = features["methylation"]))
    shell:
      "Rscript src/0_visualisation.R --met {methylation} --lattice {lattice} --c {concentration} --J {J} --r0 {r} --rep 1"


rule plotFluct:
  input:
    MC = touch(expand(outfol+"DynamicMC_c{concentration}_met-{methylation}_rep{replicate}.done",
      concentration = features["concentration"],
      methylation = features["methylation"],
      replicate = features["replicate"]))
  output:
    ovA = expand(join(outfol,"fluct/A_c{concentration}_met-{methylation}.RData",
      concentration = features["concentration"],
      methylation = features["methylation"])),
    ovM = expand(outfol+"fluct/M_c{concentration}_met-{methylation}.RData",
      concentration = features["concentration"],
      methylation = features["methylation"]),
    fluct = touch(expand(outfol+"fluct/Fluct_c{concentration}_met-{methylation}.done",
      concentration = features["concentration"],
      methylation = features["methylation"]))
  shell:
    "Rscript src/0_visualisation2.R --met {methylation} --lattice {lattice} --c {concentration} --J {J} --r0 {r}"



# ------------------------------ concentration: 0 ------------------------------
rule computePSD:
  input:
    ovA = expand(outfol+"fluct/A_c0_met-{methylation}.RData",
      methylation = features["methylation"])
  output:
    PSD = expand(outfol+"PSD/A_c0_met-{methylation}.RData",
      methylation = features["methylation"])
  shell:
    "Rscript src/1_PSD.R --met {methylation} --lattice {lattice} --J {J} --r0 {r}"



# ------------------------------ range of concentrations ------------------------------
rule dose_response:
  input:
    fluct = touch(expand(outfol+"fluct/Fluct_c{concentration}_met-{methylation}.done",
      concentration = features["concentration"],
      methylation = features["methylation"]))
  output:
    dr = expand(outfol+"dose-response/A_dr_met-{methylation}.RData",
      methylation = features["methylation"])
  shell:
    "Rscript src/1_dose-response.R --met {methylation} --lattice {lattice} --J {J} --r0 {r}"

