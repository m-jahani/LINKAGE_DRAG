library(dplyr)
library(tidyr)
library(tibble)
library(data.table)


options(scipen = 99999) # avoid scientific number

# args = commandArgs(trailingOnly = TRUE)

CHR_LENGTH_FILE <- args[1]
PC_ADMIX <- args[2]
SAMPLE_LIST <- args[3]
SAVE_DIR <- args[4]

# CHR_LENGTH_FILE <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/introg_genotyping/HA412v2_chromosome.txt"
# PC_ADMIX <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/introg_genotyping/pcadmix.xxx.regions.txt"
# SAMPLE_LIST <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/introg_genotyping/SAM_LIST"
# SAVE_DIR <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/introg_genotyping"

# read length of referenc genome for each chromosome
fread(CHR_LENGTH_FILE) %>%
    select(
        CHR = V1,
        LENGTH = V2
    ) -> CHR_LENGTH

# make 1kb intervals across genome
RANGES <- NULL
for (i in 1:nrow(CHR_LENGTH)) {
    data.frame(start = seq(1, as.numeric(CHR_LENGTH[i, 2]), 1000)) %>%
        mutate(end = as.numeric(start) + 999) %>%
        mutate(end = ifelse(end > as.numeric(CHR_LENGTH[i, 2]), as.numeric(CHR_LENGTH[i, 2]), end)) %>%
        mutate(chr = as.character(CHR_LENGTH[i, 1])) %>%
        mutate(
            chr_num = as.numeric(gsub("Ha412HOChr", "", chr)),
            cM = 0,
            ID = paste0(chr, ":", start, "-", end)
        ) %>%
        rbind(., RANGES) -> RANGES
}

rm(CHR_LENGTH)

RANGES %>%
    rbind(
        .,
        data.frame(
            start = 1, # dummy interval for avoid stopping script by samples without retrogression
            end = 1000,
            chr = "dummy",
            chr_num = 18,
            cM = 0,
            ID = "dummy:1-1000"
        )
    ) -> RANGES

# sort and write 1kb intervals
RANGES %>%
    arrange(
        as.numeric(chr_num),
        as.numeric(start)
    ) %>%
    select(
        chr,
        start,
        end
    ) %>%
    fwrite(paste0(SAVE_DIR, "/HA412_1kb_windows.bed"),
        sep = "\t",
        col.names = F,
        quote = F
    )

# read PC admix result
fread(PC_ADMIX) %>%
    separate(hap,
        into = c("sample", "haplotype"),
        sep = "_",
        remove = T
    ) %>%
    mutate(donor = ifelse(anc == 1, "ANNUUS", "2nd_GERMPLASM")) %>%
    select(-anc) -> INTROG_SAMPLE

# DONORS vector
INTROG_SAMPLE %>%
    distinct(donor) %>%
    pull(donor) -> DONORS

# HAPLOTYPES vector
INTROG_SAMPLE %>%
    distinct(haplotype) %>%
    pull(haplotype) -> HAPLOTYPES

# SAMPLES vector
fread(SAMPLE_LIST,
    header = F
) %>%
    distinct(V1) %>%
    pull(V1) -> SAMPLES


for (DONOR in 1:length(DONORS)) { # donor loops start

    # write tfam file
    fread(SAMPLE_LIST,
        header = F
    ) %>%
        mutate(
            V2 = V1,
            V3 = 0,
            V4 = 0,
            V5 = 0,
            V6 = -9
        ) %>%
        select(
            V1,
            V2,
            V3,
            V4,
            V5,
            V6
        ) %>%
        fwrite(paste0(
            SAVE_DIR,
            "/SAM_introgression_donor_",
            as.character(DONORS[DONOR]),
            ".tfam"
        ),
        sep = " ",
        col.names = F,
        quote = F
        )

    # four first column of the final result
    RANGES %>%
        select(
            chr_num,
            ID,
            cM,
            start
        ) -> RESULT

    # write map file
    RESULT %>%
        arrange(
            as.numeric(chr_num),
            as.numeric(start)
        ) %>%
        filter(ID != "dummy:1-1000") %>%
        fwrite(paste0(
            SAVE_DIR,
            "/SAM_introgression_donor_",
            as.character(DONORS[DONOR]),
            ".map"
        ),
        sep = " ",
        col.names = F,
        quote = F
        )

    for (SAMPL in 1:length(SAMPLES)) { # samples loop start
        for (HAPLO in 1:length(HAPLOTYPES)) { # haplotype loop start
            # introgressions for each donor, sample and haplotype
            INTROG_SAMPLE %>%
                filter(
                    sample == as.character(SAMPLES[SAMPL]),
                    haplotype == as.character(HAPLOTYPES[HAPLO]),
                    donor == as.character(DONORS[DONOR])
                ) %>%
                mutate(chr_num = as.numeric(gsub("Ha412HOChr", "", CHR))) %>%
                arrange(
                    as.numeric(chr_num),
                    as.numeric(start)
                ) %>%
                select(
                    CHR,
                    start,
                    end
                ) %>%
                rbind(
                    .,
                    data.frame(
                        CHR = "dummy",
                        start = 1, # dummy interval for avoid stopping script by samples without retrogression
                        end = 1000
                    )
                ) %>%
                fwrite(paste0(
                    SAVE_DIR,
                    "/",
                    SAMPLES[SAMPL],
                    "_HAP",
                    HAPLOTYPES[HAPLO],
                    "_",
                    DONORS[DONOR],
                    ".introgression.bed"
                ),
                sep = "\t",
                col.names = F,
                quote = F
                )
            # find overlap between 1kb windows and introgressions of each donor, sample and haplotype
            system(paste0(
                "bedtools intersect -wa -a ",
                SAVE_DIR,
                "/HA412_1kb_windows.bed -b ",
                SAVE_DIR,
                "/",
                SAMPLES[SAMPL],
                "_HAP",
                HAPLOTYPES[HAPLO],
                "_",
                DONORS[DONOR],
                ".introgression.bed -sorted > ",
                SAVE_DIR,
                "/",
                SAMPLES[SAMPL],
                "_HAP",
                HAPLOTYPES[HAPLO],
                "_",
                DONORS[DONOR],
                "_overlap.bed"
            ))
            # remove introgression of each loop
            system(paste0(
                "rm ",
                SAVE_DIR,
                "/",
                SAMPLES[SAMPL],
                "_HAP",
                HAPLOTYPES[HAPLO],
                "_",
                DONORS[DONOR],
                ".introgression.bed"
            ))
            # merge the overlap windows as allele 2 for making tped file
            fread(paste0(
                SAVE_DIR,
                "/",
                SAMPLES[SAMPL],
                "_HAP",
                HAPLOTYPES[HAPLO],
                "_",
                DONORS[DONOR],
                "_overlap.bed"
            )) %>%
                mutate(ID = paste0(V1, ":", V2, "-", V3)) %>%
                distinct(ID) %>%
                mutate(!!paste0(SAMPLES[SAMPL], "_", HAPLOTYPES[HAPLO]) := 1) %>% # !!paste0(SAMPLES[SAMPL],"_",HAPLOTYPES[HAPLO]) :
                # mutate(sample_hap = paste0(SAMPLES[SAMPL],"_",HAPLOTYPES[HAPLO])) %>%
                full_join(., RESULT) -> RESULT

            # remove overlaps of each loop
            system(paste0(
                "rm ",
                SAVE_DIR,
                "/",
                SAMPLES[SAMPL],
                "_HAP",
                HAPLOTYPES[HAPLO],
                "_",
                DONORS[DONOR],
                "_overlap.bed"
            ))
        } # haplotype loop end
    } # samples loop end

    # replace NA with zero.
    RESULT[is.na(RESULT)] <- 2

    RESULT %>%
        select(
            -chr_num,
            -ID,
            -cM,
            -start
        ) %>%
        colnames() %>%
        enframe() %>%
        select(value) %>%
        separate(value,
            into = c("sample", "haplotype"),
            sep = "_",
            remove = F
        ) %>%
        mutate(sample = as.numeric(gsub("SAM", "", sample))) %>%
        arrange(
            sample,
            haplotype
        ) %>%
        select(value) %>%
        pull(value) -> column_list

    RESULT %>%
        filter(ID != "dummy:1-1000") %>%
        select(
            chr_num,
            ID,
            cM,
            start,
            all_of(column_list)
        ) %>%
        arrange(
            as.numeric(chr_num),
            as.numeric(start)
        ) %>%
        fwrite(paste0(
            SAVE_DIR,
            "/SAM_introgression_donor_",
            as.character(DONORS[DONOR]),
            ".tped"
        ),
        sep = " ",
        col.names = F,
        quote = F,
        na = "NA"
        )

    rm(RESULT, column_list)
} # donor loops end

system(paste0("rm ", SAVE_DIR, "/HA412_1kb_windows.bed"))