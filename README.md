# atdp
BioWardrobe ATDP
[BioWardrobe atdp](https://github.com/cincinnati-childrens-hospital/biowardrobe/tree/master/src/atdp "BioWardrobe")

The software is used to calculate average tag density profile around all annotated TSS.
Such data can be used to estimate the success of ChIP-Seq type experiments for some histone
modifications (e.g. H3K4me)

## Usage
`atdp [options] --in=pathToFile --a=pathtoFile --out=pathToFile`

|Option          |Type    |Info                                                |
|----------------|--------|----------------------------------------------------|
|--in            |Required| Input BAM file, bam                                |
|--out           |Required| Base output file name, tsv                         |
|--a             |Required| Annotation file, tsv                               |
|--log           |Optional| Log filename [./logfile_def.log]           |
|--index         |Optional| Input index BAI file, bai                          |
|--sam_twicechr  |Optional| Which chromosome to double, str                    |
|--sam_ignorechr |Optional| Which chromosome to ignore, str                    |
|--f             |Optional| Fragmentsize, int [150]                            |
|--avd_window    |Optional| Average tag density window, int [5000]             |
|--avd_smooth    |Optional| Average smooth window (odd), int [0]               |
