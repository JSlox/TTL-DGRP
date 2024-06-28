# TTL-DGRP
```
Scripts to analyse thermal tolerance landscape of 100 lines of the DGRP as shown in Soto et al. 2024\

KO_4temp_100DGRP.xlsx ---> Reaction_norms.R\

Data file **KO_4temp_100DGRP.xlsx** contains all data from the thermal assays performed on 100 lines of the DGRP over 4 temperatures, indicating cell of the fly (celda), experimental block (bloque), date (fecha), sex, genotype ID (ID), code of genotype (Geno), experimental temperature (temp) and knockdown time on seconds (KO). This is the input file for 3 scripts:\

  - (1) **Reaction_norms.R**, that makes the reaction norms figures (1 and 2)\
  - (2) **TDT_PLOT.R**, that make the plots of the TDT lines (Figure 3)\
  - (3) **LMMs.R**, that performs the linear mixed models on Knockdown data and CTmax and z data (Tables 1,2 and 3)\

Script **Thermal_landscape.R** contains a function to calculate thermal tolerance landscape parameters (CTmax and Z)\



```


