**allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim.rds** is the list of all genes in the transcriptome, made from the realigned transcriptome data to the charr genome by Alexander.

**DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds** is the list of all DE genes identified by Alexander. However, 146 of these DE genes are not in the full gene list.
So I do not like using this file.

**DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR.rds** is the list of all DE genes identified by Alexander, minus the 146 genes that are missing from the full gene list. I prefer using this file.

Morph specific datasets:

**DEgenesFullInfo_NConly_uniq_LBonly_20k_GR.rds**, **DEgenesFullInfo_NConly_uniq_SBonly_20k_GR.rds** and **DEgenesFullInfo_NConly_uniq_PLonly_20k_GR.rds** are the list of DE genes separated by morphs. But these contain some of the 146 genes that are not in the full gene list, so I would recommend against using them.

**DEgenesFullInfo_NConly_uniq_LBonly_nomissinggene_20k_GR.rds**, **DEgenesFullInfo_NConly_uniq_SBonly_nomissinggene_20k_GR.rds** and **DEgenesFullInfo_NConly_uniq_PLonly_nomissinggene_20k_GR.rds** are the same as the three files above, but without any of the 146 genes that are not in the full gene list. I prefer using these.
