#!/usr/bin/env Rscript

# dependency

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

mixdist   <- function(Estimate, diffusion){

  if(!diffusion){
    if(length(Estimate) == 4){

      X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")
    } else
      if(length(Estimate) == 7){

        X <-  data.table(dist = seq(0, 4e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")
      } else
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 4e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  }

  if(diffusion){
    if(length(Estimate) == 4){

      X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")


    } else
      if(length(Estimate) == 7){

        X <- data.table(dist = seq(0, 10e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 10e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")

      } else
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 10e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 10e5) %>% {. / sum(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  }


  return(X[,variable := factor(variable, c("M", "E","M0", "U"))][])



}
frollconf <- function(X, Y, N){
  Ypad <- c(rep(NA, round(1.5 * N)), Y, rep(NA, round(1.5 * N)))

  data.table(ymeanpad = frollmean(Ypad, n = 2 * N + 1, align = "center", na.rm = T, hasNA = T)[],
             y95pad = frollapply(Ypad, n = 2 * N + 1, FUN = sd, align = "center", na.rm = T) * 2
  )[,.(x = X,
       ymean = ymeanpad[!is.na(Ypad)],
       y95 = y95pad[!is.na(Ypad)],
       yse = c(y95pad / sqrt(frollsum(!is.na(Ypad) * 1, n = 2 * N + 1, align = "center")))[!is.na(Ypad)])
  ]

}

dependencies <- setNames(c("_filterer",    "_overlapscore", "_filterer", "_pc",      "_pc",    "_pc", "_filterer", "_hmmfit"),
                         c("overlapscore", "sphaseorder", "pc",       "psmashc", "mmfit", "BICfit", "hmmfit", "singlefork"))

# arguments

arguments <- commandArgs(trailingOnly=TRUE)

print(arguments)

exp_pattern <- "RPE"

base_folder <-"~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/"

base_folder142 <- "~/Desktop/post-doc/Experiments/exp142-CAF1-15E-60WO-15E-F/"


base_folder136 <- "~/Desktop/post-doc/Experiments/exp136-CAF1-15E-60WO-15E-F/"



#
#all_plots <- list()

#for(i in 1:26) all_plots[[paste0(LETTERS, LETTERS)[i]]] <- "empty"

##### PLOT files
#cat("\n1. PLOT files \n")
#try({
AA <- data.table(file = list.files(base_folder, pattern = exp_pattern, recursive = T, include.dirs = F)
)[, exp := gsub(".*/|_.*|-pl.*","", file)
][, type := sub("_|-", "", gsub(paste0(exp, "|.tsv.gz|.RDS"), "", sub(".*/", "", file))), .(file)
][, size := file.size(paste0(base_folder, file)) / 1e6
][str_detect(type, "memstats|plotted", T) & str_detect(exp, "plot", T)] %>%
  ggplot() + geom_raster(aes(x = type, y = exp, fill = (size))) +
  facet_grid(cols = vars(str_detect(type, "pl")), scales = "free")
#})


##### PLOT memory usage
cat("\n2. PLOT memory usage \n")
#try({all_plots$
BB <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*memstats"), recursive = T),
                                   function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[, job := sub("_.*", "", JobID)
][][, jobtype := sub(".*\\.", "", JobName)
][][,.(mem_used = sum(as.numeric(sub("K", "", MaxRSS)), na.rm = T) / 1e6,
     time_req = Timelimit[1],
     jobtype = jobtype[1],
     time_elapsed = Elapsed[1]), .(exp, job, AllocCPUS, mem_req = as.numeric(sub("Gn", "", ReqMem)))
][][, dependency := grep(dependencies[jobtype],
                       list.files(base_folder, pattern = exp, full.names = T, recursive = T), value = T)[1], .(exp, job)
][][, mem_depend := file.size(dependency) / 1e9
][!is.na(mem_depend), depend_ratio := coef(lm(mem_used ~ mem_depend, data = .SD))[2], .(jobtype)
][] %>%
  ggplot() +
  #geom_point(aes(x = as.POSIXct(time_req), y = as.POSIXct(time_elapsed), col = jobtype)) +
  geom_point(aes(x = mem_depend, y = mem_used, col = jobtype)) +
  geom_smooth(aes(x = mem_depend, y = mem_used, col = jobtype), method = "lm", se = F) +
  geom_label(aes(x = 0.6, y = as.numeric(as.factor(jobtype)) * 2, label = paste(jobtype, round(depend_ratio, 2)), col = jobtype)) +
  #geom_point(aes(x = mem_req, y = mem_used, col = jobtype)) +
  geom_abline(slope = 1, intercept = 0)


#})


##### PLOT poisson test
CC <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*poissontest"), recursive = T),
                                   function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
) %>%
  ggplot() +
  geom_point(aes(x = log(mean), y = poisson_score, col = exp),
             size = 0.1) +
  facet_wrap(vars(exp), scales = "free") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.1) +
  geom_vline(xintercept = c(-2.5, 0))


##### PLOT posterior per experiment
DD <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T),
                                   function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[,.(N = as.double(.N)), .(exp, posterior)
][,N := N / max(N), .(exp)] %>%
  ggplot() + geom_tile(aes(x = 1, y = posterior, height = 1, width = N)) +
  facet_wrap(vars(exp))
DD

##### PLOT shpase order QCs
EE <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T),
                                   function(x){readRDS(paste0(base_folder, x))$bootstrap[, exp := gsub(".*/|_.*", "", x)]})
) %>%
  ggplot() + geom_point(aes(x = rank, y = value, col = N), size = 0.2) +
  facet_wrap(vars(exp)) +
  scale_color_viridis_c()
EE

FF <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T),
                                   function(x){as.data.table(readRDS(paste0(base_folder, x))$final_order)[, exp := gsub(".*/|_.*", "", x)]})
) %>% ggplot(aes(y = X1.1, x = X2.1)) +
  geom_segment(aes(yend = X2, xend = Y2), size = 0.05) +
  geom_point(aes(col = rank), size = 1) +
  coord_fixed() +
  facet_wrap(vars(exp)) +
  scale_color_viridis_c(name = "Tour index")
FF
####### plotting cut tracks


RepTiming = list(
  fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/4D_early_8R_UC.tsv.gz")%>% setNames(c("chr", "bp"))%>% mutate(round = round(bp / 10e4) * 10e4,
                                                                                                               RT= "early") %>%
    group_by(chr, round, RT) %>%
    summarize(count = dplyr::n()),

  fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/4D_late_3O_JH.tsv.gz")%>% setNames(c("chr", "bp"))%>% mutate(round = round(bp / 10e4) * 10e4,
                                                                                                              RT= "late") %>%
    group_by(chr, round, RT) %>%
    summarize(count = dplyr::n())
)%>% rbindlist(idcol = F)

RT1=RepTiming %>% dplyr::filter(chr == "chr2", round>2e7 & round<6e7)%>%
  pivot_wider(names_from = RT, values_from = count, values_fill = 0)%>% mutate(total = early+late)%>%
  dplyr::filter(total>10)%>%
  mutate(reptiming= log2(early/late))%>% ggplot()+geom_raster(aes(x=round, y=0, fill=reptiming))+theme_bw()+
  scale_fill_distiller(palette = "RdYlBu", name= "Replication Timing")+coord_cartesian(expand = F)

RT2=fread("~/Desktop/post-doc/Experiments/datasets_RPE1/Halazonetis_Nature_RPE1_origins/RPE1.tsv.gz") %>%
  rename(chr= V1, bp = V2) %>% mutate(RT = "origins") %>% group_by(chr) %>%
  mutate(round = round(bp / 5e4) * 5e4) %>%
  group_by(chr, round) %>%
  summarize(count = dplyr::n())%>% filter(chr==2, round>2e7, round<6e7)%>%
  mutate(count = (count - quantile(count,0.05)) / (quantile(count, 0.95) - quantile(count,0.05)))%>%
  ggplot()+geom_col(aes(x=round, y=count), fill="red", width = 10e4)+theme_bw()+
  coord_cartesian(expand = F)+ylim(0,3)


momics        <- fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/scEU_rpf_chic_pb_gw_10kb_class.tsv.gz")
gene_location <- fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/hg38_gene_location.tsv")


momics_track <- momics[chr == "2" & round < 6e7 & round > 2e7 & (filterbin) & cluster2 == 3] %>% select(-rpf) %>%
  pivot_longer(cols = c(k27me3, k9me3, k36me3, nascent), names_to = "sample") %>%
  ggplot()+
  #geom_segment(data = gene_location[chr == "2" & to < 5e7],
  #             aes(x = from, xend = to, y = 0, yend = 0), size = 1.5) +
  geom_area(aes(x = round, y = value, fill = sample))+
  theme_bw()+
  facet_grid(rows = (vars(sample)),scales = "free")+
  scale_fill_manual(values=c( "#D674BA", "#86BF88", "#7A416A" , "cyan3"))+
  ylab("Z-score")+
  coord_cartesian(expand=F)+
  xlab("Chromosome 2 (bp)")+
  theme(legend.position = "none")


GG <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T),
                                   function(x){
                                     fread(paste0(base_folder, x)
                                     )[chr == 2 & posterior > 700
                                     ][, exp := gsub(".*/|_.*", "", x)
                                     ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                                     ][, .(N = as.double(.N)), .(exp, rank, bp = round(bp, -5))
                                     ][, N := (N) / quantile(N, 0.95), .(rank)
                                     ][bp > 2e7 & bp < 6e7
                                     ][N > 1, N := 1]
                                   })
)[str_detect(exp, "DMSO|24"), rank := abs(rank - 1)  #"DMSO-2hr|V1-2hr"
]%>% #dplyr::filter(str_detect(exp, "DMSO|2hr"))%>%
  ggplot() +
  geom_raster(aes(y = rank, x = bp / 1e6, fill = (N)), interpolate = T) +
  scale_y_continuous(trans = "reverse") + ylab("S-phase Progression") + xlab("Chromosome 2 (Mbp)") +
  facet_wrap(vars(fct_rev(exp)), scales = "free", ncol = 1) +
  scale_fill_gradientn(colours = c("white", "blue", "blue4")) +
  theme_bw() +
  coord_cartesian(expand = F)

GG
(momics_track/RT1/RT2/GG)+plot_layout(heights = c(6,1,1, 14))

####### plotting genome pile up
HH <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T),
                                   function(x){
                                     fread(paste0(base_folder, x)
                                     )[chr %in% 1:22, .(N = as.double(.N), exp = gsub(".*/|_.*", "", x)), .(cell, chr, bp = round(bp / 1e6))
                                     ][, N := N / sum(N), .(cell)
                                     ][, .(N = sum(N)), .(exp, chr, bp)
                                     ][, N := N / quantile(N, 0.97), .(chr)
                                     ][N > 1, N := 1][]
                                   })) %>%
  ggplot() +
  geom_col(aes(y = N, x = bp, fill = log10(N)) ) +
  facet_grid(rows = vars(exp), cols = vars(chr), scales = "free") +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F)
HH

####### plotting model fits
II <- merge.data.table(
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                function(x){fread(paste0(base_folder, x)
                )[, exp := gsub(".*/|_.*", "", x)
                ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                ]})
  )[((!diffusion) & model != "quadruple") | (diffusion)
  ][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
  ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
  ][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
  ][parameter != 11, mixdist(Estimate, diffusion = diffusion[1]), .(exp, cell, rank, diffusion)],
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*psmashc"), recursive = T),
                function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})),
  by = c("exp", "cell", "dist")
)[dist > 0
][,  N := as.double(N)
][, N := N / sum(N), .(exp, cell, diffusion, variable)
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
  #][rank > 0.4 & rank < 0.6
][,.SD[cell %in% sample(unique(cell), 3)], .(exp)
][dist > 0 & dist < 4e5] %>%
  ggplot() +
  geom_col(aes(x = dist, y = value, fill = variable)) +
  geom_line(aes(x = dist , y = N)) +
  facet_wrap(vars(exp, rank), scales = "free", nrow = 4)
#})

II

##### PLOT speeds
JJ <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                                   function(x){fread(paste0(base_folder, x)
                                   )[, exp := gsub(".*/|_.*", "", x)
                                   ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                                   )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                                   ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.8
  #][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>%
  {
    .[!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>%
      ggplot() +
      geom_point(data = ., aes(x = rank, y = Estimate )) +
      geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95), alpha = 0.3) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse), alpha = 0.3) +
      geom_line(aes(x = x, y = ymean)) +
      facet_grid(cols = vars(exp)) + theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")

  }

JJ

KK <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                                   function(x){fread(paste0(base_folder, x)
                                   )[, exp := gsub(".*/|_.*", "", x)
                                   ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                                   )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                                   ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.8

][!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>%
  ggplot() +
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, col = exp), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill = exp), alpha = 0.3) +
  geom_line(aes(x = x, y = ymean, col = exp)) +
  theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")

KK

pc_sim <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/pc_sim_super_all.tsv") %>%
  dplyr::filter(U !="U")%>%
  mutate_if(is.character, as.numeric)%>% as.data.table()

str(pc_sim)
pc_sim[, c("u_R", "s_R")       := .(abs(u - mean(u)) / sd(u), abs(s - mean(s)) / sd(s)), .(U,S)
][, c("u_rank", "s_rank") := .(rank(u_R), rank(s_R)), .(U,S)]

U_model <- loess(U ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])


LL1 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                                   function(x){fread(paste0(base_folder, x)
                                   )[, exp := gsub(".*/|_.*", "", x)
                                   ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                                   )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                                   ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|2hr"))%>%
  {
    .[!is.na(rank)][order(rank), frollconf(rank, u, 25), .(exp)]  %>%
      ggplot() +
      geom_point(data = ., aes(x = rank, y = u, col =exp), size=1.75, alpha=0.5) +
      #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
      geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
      #facet_grid(cols = vars(fct_rev(exp))) +
      theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")+
      scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
      scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
      theme(legend.position = "none")+ylim(0,2)#+coord_cartesian()

  }



LL2 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|2hr"))%>%
  ggplot() +
  geom_histogram(aes( y = u, fill =exp), position = position_identity(), alpha=0.5) +
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  #facet_grid(cols = vars(fct_rev(exp))) +
  theme_bw() + ylab("") + xlab("")+
  scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
  scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
  theme(legend.position = "none")+ylim(0,2)#+coord_cartesian()

library(patchwork)
(LL1|LL2)+plot_layout(widths = c(4,1.5))

LL3 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|10m")) %>%
  mutate(rankbin = (round(rank*30))/30)%>%
  group_by(exp,rankbin)%>%
  summarise(meanu = mean(u))%>% pivot_wider(names_from = exp, values_from = meanu) %>%
  mutate(speed_diff = (`JvB-140-RPE-CAF1-2hr-dTAG-15E-60WO-15E-F` - `JvB-140-RPE-CAF1-DMSO-15E-60WO-15E-F`))%>%
  ggplot()+geom_point(aes(x=rankbin, y=speed_diff), col="#FF7B15")+
  geom_smooth(aes(x=rankbin, y=speed_diff), span=0.4, level=0.75, col="#FF7B15", fill="#FF7B15", alpha=0.5)+
  xlim(c(-.05,0.8))+
  ylim(c(-1, 0))+
  theme_bw()+
  ylab("Speed difference (kb/min)")+xlab('')


((LL3|LL3)/(LL1|LL2))+plot_layout(widths = c(4,1.5), heights = c(1.5,6))

LLa <- (LL3/LL1)+plot_layout(widths = c(4,1.5), heights = c(1.25,6))
LLb <- (LL2/LL2)+plot_layout(widths = c(4,1.5), heights = c(1.25,6))
(LLa|LLb)+plot_layout(widths = c(4,1.5))


GG

base_folder142 <- "~/Desktop/post-doc/Experiments/exp142-CAF1-15E-60WO-15E-F/"

GG|GG142


GG142 <- rbindlist(map(list.files(base_folder142, pattern = paste0(exp_pattern, ".*filterer"), recursive = T),
                    function(x){
                      fread(paste0(base_folder142, x)
                      )[chr == 2 & posterior > 700
                      ][, exp := gsub(".*/|_.*", "", x)
                      ][, rank := as.data.table(readRDS(list.files(base_folder142, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                      )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                      ][, .(N = as.double(.N)), .(exp, rank, bp = round(bp, -5))
                      ][, N := (N) / quantile(N, 0.95), .(rank)
                      ][bp > 2e7 & bp < 6e7
                      ][N > 1, N := 1]
                    })
)[str_detect(exp, "RPE-CAF1-dTAG-15E|RPE-CAF1-DMSO-15E"), rank := abs(rank - 1)
]%>% dplyr::filter(str_detect(exp, "CAF1-dTAG-15E|CAF1-DMSO-15E"))%>%
  ggplot() +
  geom_raster(aes(y = rank, x = bp / 1e6, fill = (N)), interpolate = T) +
  scale_y_continuous(trans = "reverse") + ylab("S-phase Progression") + xlab("Chromosome 2 (Mbp)") +
  facet_wrap(vars((exp)), scales = "free", ncol = 1) +
  scale_fill_gradientn(colours = c("white", "blue", "blue2", "blue4")) +
  theme_bw() +
  coord_cartesian(expand = F)


GG142

GG|GG142


LL1 <- rbindlist(map(list.files(base_folder142, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                     function(x){fread(paste0(base_folder142, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder142, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "RPE-CAF1-dTAG-15E|RPE-CAF1-DMSO-15E"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "CAF1-dTAG-15E|CAF1-DMSO-15E"))%>%
  {
    .[!is.na(rank)][order(rank), frollconf(rank, u, 25), .(exp)]  %>%
      ggplot() +
      geom_point(data = ., aes(x = rank, y = u, col =fct_rev(exp)), size=1.75, alpha=0.5) +
      #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=fct_rev(exp)), alpha = 0.2) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=fct_rev(exp)), alpha = 0.2) +
      geom_line(aes(x = x, y = ymean, col=fct_rev(exp)), size = 1) +
      #facet_grid(cols = vars(fct_rev(exp))) +
      theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")+
      scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
      scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
      theme(legend.position = "none")+ylim(0,2.5)#+coord_cartesian()

  }



LL2 <-rbindlist(map(list.files(base_folder142, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                    function(x){fread(paste0(base_folder142, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder142, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "RPE-CAF1-dTAG-15E|RPE-CAF1-DMSO-15E"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "CAF1-dTAG-15E|CAF1-DMSO-15E"))%>%
  ggplot() +
  geom_histogram(aes( y = u, fill =fct_rev(exp)), position = position_identity(), alpha=0.5) +
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  #facet_grid(cols = vars(fct_rev(exp))) +
  theme_bw() + ylab("") + xlab("")+
  scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
  scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
  theme(legend.position = "none")+ylim(0,2.5)#+coord_cartesian()



LL3 <- rbindlist(map(list.files(base_folder142, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                     function(x){fread(paste0(base_folder142, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder142, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "RPE-CAF1-dTAG-15E|RPE-CAF1-DMSO-15E"), rank := abs(rank - 1)
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "CAF1-dTAG-15E|CAF1-DMSO-15E"))%>%
  mutate(rankbin = (round(rank*30))/30)%>%
  group_by(exp,rankbin)%>%
  summarise(meanu = mean(u))%>% pivot_wider(names_from = exp, values_from = meanu) %>%
  mutate(speed_diff = (`JvB-142-RPE-CAF1-dTAG-15E-60WO-15E-F` - `JvB-142-RPE-CAF1-DMSO-15E-60WO-15E-F`))%>%
  ggplot()+geom_point(aes(x=rankbin, y=speed_diff), col="#FF7B15")+
  geom_smooth(aes(x=rankbin, y=speed_diff), span=0.5, level=0.75, col="#FF7B15", fill="#FF7B15", alpha=0.5)+
  xlim(c(-.05,0.8))+
  ylim(c(-1, 0.1))+
  theme_bw()+
  ylab("Speed difference (kb/min)")+xlab('')


((LL3|LL3)/(LL1|LL2))+plot_layout(widths = c(4,1.5), heights = c(1.5,6))

LLa <- (LL3/LL1)+plot_layout(widths = c(4,1.5), heights = c(1.25,6))
LLb <- (LL2/LL2)+plot_layout(widths = c(4,1.5), heights = c(1.25,6))
(LLa|LLb)+plot_layout(widths = c(4,1.5))






base_folder136 <- "~/Desktop/post-doc/Experiments/exp136-DRB_PARP-CAF1-15E-60WO-15E-F/"






forks_DMSO <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-140-RPE-CAF1-DMSO-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-140-RPE-CAF1-DMSO-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'DMSO'
]%>% mutate(rank = abs(rank-1))%>% as.data.table()

forks_10 <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-10m-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-10m-dTAG-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-10m-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-10m-dTAG-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := '10min'
]


forks_2 <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-2hr-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-2hr-dTAG-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-2hr-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-2hr-dTAG-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := '2hr'
]

forks_DMSO[, bin := round(center / 2e5) * 2e5]
forks_DMSO[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_DMSO[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_DMSO[order(rank), cumreads := cumsum(N), .(chr, bin)]


forks_10[, bin := round(center / 2e5) * 2e5]
forks_10[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_10[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_10[order(rank), cumreads := cumsum(N), .(chr, bin)]

forks_2[, bin := round(center / 2e5) * 2e5]
forks_2[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_2[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_2[order(rank), cumreads := cumsum(N), .(chr, bin)]




forks_30 <- fread("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-24hr-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-24hr-dTAG-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp140_CAF1_10m-2hr-30hr-15E-60WO-15E-F/scEdU-seq/JvB-140-RPE-CAF1-24hr-dTAG-15E-60WO-15E-F/JvB-140-RPE-CAF1-24hr-dTAG-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'dTAG'
]%>% mutate(rank = abs(rank-1))%>% as.data.table()



forks_30[, bin := round(center / 2e5) * 2e5]
forks_30[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_30[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_30[order(rank), cumreads := cumsum(N), .(chr, bin)]



#forks_plus8 <- bind_rows(forks_DMSO, forks_30)
forks_plus8 <- bind_rows(forks_DMSO, forks_10)

fork_8_10 <- forks_plus8[chr == 2 & bin > 1e7 & bin < 6e7 & width_est <2e6 & width_est >2.5e4 & N>2] %>%
  ggplot() +
  geom_segment(aes(x = (center - width_est / 2) / 1e6, xend = (center + width_est / 2) / 1e6, y = abs(rank-1), yend = abs(rank-1)),
                   color ="blue2", size =0.2, alpha =0.8) +
  ylab("S-phase progression") +
  xlab("Genomic location [Mb]") +
  coord_cartesian(expand = F) + scale_fill_manual("black")+
  theme_bw()+facet_grid(rows = vars(fct_rev(treatment)))#+ggtitle('Clone 8')

(RT1/RT2/fork_8_10)+plot_layout(heights = c(1,1, 9))


library(patchwork)
gnome

forks_plus8 %>% #mutate(rank = rank/max(rank))%>%
  dplyr::filter(rank< 0.8)%>%
  group_by(treatment, rank)%>%
  summarize(n= n())%>%
  select(treatment, n) %>%
  dplyr::filter(n < 4500)%>%
  ungroup() %>%
  rstatix::pairwise_t_test(n ~ treatment)




forks_6b3_DMSO <- fread("~/Desktop/post-doc/Experiments/exp136-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp136-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'DMSO'
]






forks_6b3_dTAG <- fread("~/Desktop/post-doc/Experiments/exp136-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp136-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'dTAG'
]%>% mutate(rank = abs(rank))%>% as.data.table()



forks_6b3_DMSO[, bin := round(center / 2e5) * 2e5]
forks_6b3_DMSO[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_6b3_DMSO[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_6b3_DMSO[order(rank), cumreads := cumsum(N), .(chr, bin)]

forks_6b3_dTAG[, bin := round(center / 2e5) * 2e5]
forks_6b3_dTAG[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_6b3_dTAG[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_6b3_dTAG[order(rank), cumreads := cumsum(N), .(chr, bin)]

forks_plus6b3 <- bind_rows(forks_6b3_DMSO, forks_6b3_dTAG)



fork_6b3 <- forks_plus6b3[chr == 2 & bin > 2e7 & bin < 8.5e7 & width <2e6] %>%
  ggplot() +
  geom_segment(aes(x = (center - width_est / 2) / 1e6, xend = (center + width_est / 2) / 1e6, y = abs(rank-1), yend = abs(rank-1)),
               color ="grey10", size =0.25) +
  ylab("S-phase progression") +
  xlab("Genomic location [Mb]") +
  coord_cartesian(expand = F) + scale_fill_manual("black")+
  theme_bw()+facet_grid(rows = vars(treatment))+ggtitle('Clone 6B3')


fork_8|fork_6b3

forks_plus8 <- bind_rows(forks_DMSO, forks_30)

fork_8 <- forks_plus8[chr == 2 & bin > 2e7 & bin < 8.5e7 & width <2e6] %>%
  ggplot() +
  geom_segment(aes(x = (center - width_est / 2) / 1e6, xend = (center + width_est / 2) / 1e6, y = abs(rank-1), yend = abs(rank-1)),
               color ="grey10", size =0.25) +
  ylab("S-phase progression") +
  xlab("Genomic location [Mb]") +
  coord_cartesian(expand = F) + scale_fill_manual("black")+
  theme_bw()+facet_grid(rows = vars(treatment))+ggtitle('Clone 8')



min10a <- forks_plus %>% dplyr::filter(rank < 0.33) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/247,
         `10min` = `10min`/249)%>% select(loc, DMSO, `10min`)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`10min`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  xlim(0,0.5)+ylim(0,0.5)+ggtitle('Early')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+
  xlab('Norm. fork coverage \n per bin (200kb, DMSO)')+ ylab('Norm. fork coverage \n per bin (200kb, dTAG)')+
  coord_fixed(expand = F, ratio =1)+theme(legend.position = 'none')

min10b <- forks_plus %>% dplyr::filter(rank >0.33 & rank <0.66) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/249,
         `10min` = `10min`/248)%>% select(loc, DMSO, `10min`)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`10min`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+
  xlim(0,0.6)+ylim(0,0.6)+ggtitle('Mid')+
  xlab('Norm. fork coverage \n per bin (200kb, DMSO)')+ ylab('Norm. fork coverage \n per bin (200kb, dTAG)')+
  theme_bw()+coord_fixed(expand = F, ratio =1)+theme(legend.position = 'none')

min10c <- forks_plus %>% dplyr::filter(rank >0.66) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/255,
         `10min` = `10min`/256)%>% select(loc, DMSO, `10min`)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`10min`), bins = 30)+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+
  xlim(0,0.4)+ylim(0,0.4)+ ggtitle('Late')+
  xlab('Norm. fork coverage \n per bin (200kb, DMSO)')+ ylab('Norm. fork coverage \n per bin (200kb, dTAG)')+
  theme_bw()+coord_fixed(expand = F, ratio =1)+theme(legend.position = 'right')

(min10a|min10b|min10c)
corrplot

Late <- forks_plus %>% dplyr::filter(rank >0.66) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/250,
         `10min` = `10min`/250)%>% select(loc, DMSO, `10min`)%>%
  setNames(c("loc", "Late_DMSO", "Late_10min"))%>% as.data.table()


Mid <- forks_plus %>% dplyr::filter(rank < 0.66 & rank > 0.33 ) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/250,
         `10min` = `10min`/250)%>% select(loc, DMSO, `10min`)%>%
  setNames(c("loc", "Mid_DMSO", "Mid_10min"))%>% as.data.table()

Early <- forks_plus %>% dplyr::filter(rank < 0.33 ) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/250,
         `10min` = `10min`/250)%>% select(loc, DMSO, `10min`)%>%
  setNames(c("loc", "Early_DMSO", "Early_10min")) %>% as.data.table()


EM <- merge(Early, Mid, by = "loc" )

All <- merge(EM, Late, by = "loc"  )

corrplot <- All %>% ungroup() %>% select(-loc)%>% drop_na()%>%
  cor(method="pearson") %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  #mutate(rowname = factor(rowname,levels = c("ref_early","DOT_early","ref_late","DOT_late")),
  #       name = factor(name,levels = rev(c("ref_early","DOT_early","ref_late","DOT_late")))) %>%
  ggplot(aes(x=factor(rowname, c("Early_DMSO",
                              "Early_10min",
                              "Mid_DMSO",
                              "Mid_10min",
                              "Late_DMSO",
                              "Late_10min")),
             y= factor(name, c("Early_DMSO",
                              "Early_10min",
                              "Mid_DMSO",
                              "Mid_10min",
                              "Late_DMSO",
                              "Late_10min"))))+
  geom_tile(aes(fill=value),color="black", size=0.35)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "Pearson \nCorrelation") +
  geom_text(aes(label=round(value,2)), size=5)+coord_fixed()+
  #guides(shape = guide_legend(override.aes = list(size = 1.5)),
  #       color = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.position = 'right')+ggtitle('Replication Timing')



forks_plus %>% dplyr::filter(rank<0.33 | rank>0.66) %>%
  mutate(early = rank < 0.33) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment, early) %>%
  summarize(n= n()) %>%
  group_by(early)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = early, values_from = n)%>%
  ggplot()+
  geom_bin2d(aes(x=`TRUE`, y=`FALSE`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+
  #xlim(0,0.4)+ylim(0,0.4)+
  ggtitle('Early vs Late Correlation')+
  xlab('Early')+ ylab('late')+
  theme_bw()+coord_fixed(expand = F, ratio = 1)+theme(legend.position = 'none')

OO2 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                     function(x){fread(paste0(base_folder, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 5e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|10m")) %>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp), y=s, fill=fct_rev(exp)), outlier.shape = NA, alpha=0.25)+
  geom_jitter(aes(x=fct_rev(exp), y=s, col=fct_rev(exp)), size =0.2, alpha=0.5)+
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  #facet_grid(cols = vars(fct_rev(exp))) +
  theme_bw() + ylab("Variance (kb/min)") + xlab('S-phase Progression')+
  scale_fill_manual(values = c("#1F97FF", "#FF7B15"))+
  scale_color_manual(values = c("#1F97FF", "#FF7B15"))+
  theme(legend.position = "none")+
  #xlim(0,0.8)+
  ylim(0.15,0.45)


rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
              function(x){fread(paste0(base_folder, x)
              )[, exp := gsub(".*/|_.*", "", x)
              ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
              )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
              ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 5e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|10m")) %>%
  rstatix::pairwise_t_test(s ~ exp)


OO3 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                     function(x){fread(paste0(base_folder, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 5e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% dplyr::filter(str_detect(exp, "DMSO|10m")) %>%
  ggplot()+
  geom_point(aes(x=rank, y=s, col=exp), alpha=0.5)+
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
  facet_grid(cols = vars(fct_rev(exp))) +
  theme_bw() + ylab("Variance (kb/min)") + xlab('S-phase Progression')+
  scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
  scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
  theme(legend.position = "none")+
  xlim(0,0.8)+
  ylim(0.15,0.45)

OO4 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                     function(x){fread(paste0(base_folder, x)
                     )[, exp := gsub(".*/|_.*", "", x)
                     ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                     )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                     ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 5e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
  ] %>% dplyr::filter(str_detect(exp, "DMSO|10m")) %>%
  ggplot()+
  geom_point(aes(x=u, y=s, col=exp), alpha=0.5)+
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, fill=exp), alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill=exp), alpha = 0.2) +
  #geom_line(aes(x = x, y = ymean, col=exp), size = 1) +
#  facet_grid(cols = vars(fct_rev(exp))) +
  theme_bw() + ylab("Variance (kb/min)") + xlab('DNA replication speed (kb/min)')+
  scale_fill_manual(values = c("#FF7B15", "#1F97FF"))+
  scale_color_manual(values = c("#FF7B15",  "#1F97FF"))+
  theme(legend.position = "none")+
  #xlim(0,0.8)+
  ylim(0.15,0.45)

(OO2|OO4|OO3)+plot_layout(width =c(1.5,4,6))

forks_plus %>% dplyr::filter(rank < 0.33) %>%  group_by(treatment) %>% distinct(rank) %>% summarize(n=n())
#min10

hr2 <- forks_plus2 %>% dplyr::filter(rank < 0.3) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/114,
         `2hr` = `2hr`/105)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`2hr`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+coord_cartesian(expand = F)




hr24 <- forks_plus24 %>% dplyr::filter(rank < 0.3) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/227,
         `24hr` = `24hr`/211)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`24hr`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+coord_cartesian(expand = F)

min10|hr2|hr24

hr2LATE <- forks_plus2 %>% dplyr::filter(rank > 0.7) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/114,
         `2hr` = `2hr`/105)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`2hr`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+coord_cartesian(expand = F)+theme(legend.position = 'none')


min10LATE <- forks_plus %>% dplyr::filter(rank < 0.7) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/114,
         `10min` = `10min`/112)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`10min`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+theme(legend.position = 'none')+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+coord_cartesian(expand = F)+theme(legend.position = 'none')

hr24LATE <- forks_plus24 %>% dplyr::filter(rank < 0.7) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  group_by(treatment)%>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n)%>%
  mutate(DMSO = DMSO/151,
         `24hr` = `24hr`/141)%>%
  ggplot()+
  geom_bin2d(aes(x=DMSO, y=`24hr`))+
  geom_abline(slope = 1, intercept = 0, size=1, linetype=2, col= "blue")+
  scale_fill_viridis_c(option="B", trans = "log10", limits = c(1, 150), oob = scales::squish)+theme_bw()+coord_cartesian(expand = F)

(min10|hr24)/(min10LATE|hr24LATE)


forks_corr <- forks_plus24 %>% dplyr::filter(rank < 0.66) %>%
  group_by(chr, treatment) %>%
  mutate(bin  = round(center / 2e5) * 2e5)%>%
  group_by(chr, bin, treatment) %>%
  summarize(n= n()) %>%
  mutate(loc = paste0(chr, "-", bin))%>%
  pivot_wider(names_from = treatment, values_from = n, values_fill =  0)%>%
  mutate(DMSO = DMSO/151,
         `24hr` = `24hr`/141)%>%
  ungroup()%>%  select(DMSO, `24hr`)

cor(forks_corr$DMSO, forks_corr$`24hr`, method = "spearman")


meth_cuts_reads <- readRDS("~/Desktop/post-doc/ChiC/RPE-1_rawdata/cuts_reads.RDS")#[[1]][, max(bp) - min(bp), .(V1, chr, bin)]

# which bin in which mod
mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
  x[, .(N = as.double(.N)), .(V1, chr, bin = round(bp / 1.5e4)*1.5e4)][, N := N / sum(N)]})
)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = sum(N)), .(mod, chr, bin)]


mod_bin2 <- mod_bin %>% group_by(mod) %>% mutate(score = scale(N))%>%
  ungroup()%>% group_by(chr,bin)%>%
                                                      mutate(max = max(score),
                                                      maxN = score/max)%>%
  dplyr::filter(maxN == 1) %>% select(chr, bin, mod) %>%
  mutate(loc= paste0(chr, "-", bin))%>%
  as.data.table()


mod_bin2 %>% group_by(mod) %>%
  summarize(n=n())%>%
  ungroup()%>%
  mutate(max = sum(n),
         Nnorm = n/max)%>%
  ggplot()+geom_col(aes(x=mod, y=Nnorm*100, fill=factor(mod, c('k9me3','k27me3','k36me3'))))+
  #geom_text(aes(y = ypos, label = mod), color = "white", size=6) +
  scale_fill_manual(values= c(   "#7A416A", "#D674BA",  "#86BF88"))+
  scale_color_manual(values= c(   "#D674BA","#86BF88", "#7A416A"  ))+theme_bw()+
  theme(legend.position = 'none')+
  ylab("Percentage of 15kb bins")+
  xlab('')



forks_DMSO2 <- forks_DMSO[, bin := round(center / 2e5) * 2e5] %>%
  group_by(chr, bin) %>% summarise(medrank = mean(rank)) %>%
  mutate(loc= paste0(chr, "-", bin))%>%
  as.data.table()


merge.data.table(mod_bin2, forks_DMSO2, by = "loc")%>%
  ggplot()+
  geom_jitter(aes(y=factor(mod, c('k9me3', 'k27me3','k36me3')), x=medrank, col = factor(mod, c('k9me3', 'k27me3','k36me3'))), alpha=0.5, size =0.20)+
  geom_violin(aes(y=factor(mod, c('k9me3', 'k27me3','k36me3')), x=medrank, fill= factor(mod, c('k9me3', 'k27me3','k36me3'))),draw_quantiles = T, alpha=0.25, scale = "area")+
  coord_cartesian()+
  theme_bw()+theme(legend.position = 'none')+
  scale_fill_manual(values= c(  "#7A416A",  "#D674BA", "#86BF88" ))+
  scale_color_manual(values= c( "#7A416A", "#D674BA","#86BF88" ))+
  ylab('')+
  xlab('S-phase Progression')



d1 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)%>% dplyr::filter(!is.na(rank))%>%
  dplyr::filter(str_detect(exp, "DMSO|10m"))%>%
  distinct(cell, exp)%>% mutate(plate = ceiling(cell/384)) %>% group_by(plate, exp) %>% summarize(rank=n())

d2 <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                    function(x){fread(paste0(base_folder, x)
                    )[, exp := gsub(".*/|_.*", "", x)
                    ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                    )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                    ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2.5 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
]%>%  dplyr::filter(str_detect(exp, "DMSO|10m"))%>% mutate(plate = ceiling(cell/384)) %>% distinct(plate, cell, exp)%>% group_by(plate, exp) %>% summarize(speed=n())


OO <- left_join(d1, d2, by=c("exp", "plate"))%>%
  mutate(ratio = speed/rank, ratioEdU = rank/1520)%>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp),y=ratio, fill= fct_rev(exp)),alpha=0.5)+
  geom_jitter(aes(x=fct_rev(exp),y=ratio, col= fct_rev(exp)))+
  theme_bw()+theme(legend.position = "none")+ ylim(0,1)+
  scale_fill_manual(values = c("#1F97FF", "#FF7B45" ))+
  scale_color_manual(values = c( "#1F97FF", "#FF7B45"))+
  ggtitle("Ratio of RFS+/EdU+")+xlab('')+ylab('Percentage of cells')


OOO <- left_join(d1, d2, by=c("exp", "plate"))%>%
  mutate(ratio = speed/rank, ratioEdU = rank/376)%>%
  ggplot()+
  geom_boxplot(aes(x=fct_rev(exp),y=ratioEdU, fill= fct_rev(exp)),alpha=0.5)+
  geom_jitter(aes(x=fct_rev(exp),y=ratioEdU, col= fct_rev(exp)))+
  theme_bw()+theme(legend.position = "none")+ ylim(0,1)+
  scale_fill_manual(values = c("#1F97FF", "#FF7B45" ))+
  scale_color_manual(values = c( "#1F97FF", "#FF7B45"))+
  ggtitle("QC PASS")+xlab('')+ylab('Percentage of cells')


PP<- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T),
                   function(x){fread(paste0(base_folder, x)
                   )[, exp := gsub(".*/|_.*", "", x)
                   ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                   )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                   ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4
][str_detect(exp, "DMSO|24"), rank := abs(rank - 1)
][rank < 0.9
][, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
][] %>% dplyr::filter(str_detect(exp, "DMSO|10m"))%>%
  mutate(rankbin = (round(rank*100))/100)%>%
  group_by(exp)%>%
  arrange(rank)%>%
  mutate(cumulative = cumsum(rep(1,n()))/n())%>%
  ggplot()+geom_line(aes(x=rank, y=cumulative,
                         col=factor(exp)),
                     size=1.5)+theme_bw()+
  ylab('Cumulative RFS+ cells')+xlab('S-phase Progression')+
  #scale_fill_manual(values = c("#6D6E71", "#007D20"))+
  #scale_color_manual(values = c("#6D6E71", "#007D20"))+coord_cartesian()+
  theme(legend.position = 'bottom')+coord_equal()+
  theme_bw()+
  scale_fill_manual(values = c("#1F97FF", "#FF7B45" ))+
  scale_color_manual(values = c( "#1F97FF", "#FF7B45"))+ theme(legend.position = "none")


OOO|OO|PP




forks_DMSO <- fread("~/Desktop/post-doc/Experiments/exp136-DRB_PARP-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp136-DRB_PARP-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F/JvB-136-RPE-CAF1-DMSO-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := 'DMSO'
]%>% mutate(rank = abs(rank-1))%>% as.data.table()

forks_10 <- fread("~/Desktop/post-doc/Experiments/exp136-DRB_PARP-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F_singlefork.tsv.gz"
)[state == 2
][, width_est := width + width / N
][order(center)
][, rank := as.data.table(readRDS("~/Desktop/post-doc/Experiments/exp136-DRB_PARP-CAF1-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F/JvB-136-RPE-CAF1-dTAG-15E-60WO-15E-F_sphaseorder.RDS"
)$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]][, treatment := '10min'
]




forks_DMSO[, bin := round(center / 2e5) * 2e5]
forks_DMSO[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_DMSO[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_DMSO[order(rank), cumreads := cumsum(N), .(chr, bin)]


forks_10[, bin := round(center / 2e5) * 2e5]
forks_10[order(rank), cumfire := cumsum(rep(1, .N)), .(chr, bin)]
forks_10[, cumfire := cumfire / uniqueN(forks_plus$cell)]
forks_10[order(rank), cumreads := cumsum(N), .(chr, bin)]

forks_plus <- bind_rows(forks_DMSO, forks_10)


forksa <- forks_plus %>% #mutate(rank = rank/max(rank))%>%
  #dplyr::filter(rank< 0.8)%>%
  group_by(treatment, rank)%>%
  summarize(n= n())%>%
  ggplot()+geom_violin(aes(x = fct_rev(treatment),
                           y=n,
                           fill= fct_rev(treatment)), alpha =0.3)+
  geom_jitter(aes(x = fct_rev(treatment),
                  y=n,
                  col= fct_rev(treatment)), alpha=0.5, size=0.75)+
  ylim(0,4500)+ xlab("DNA replication forks per cell")+ylab("DNA replication forks per cell")+
  theme_bw()+
  scale_fill_manual(values = c("#1F97FF", "#FF7B45" ))+
  scale_color_manual(values = c( "#1F97FF", "#FF7B45"))+ theme(legend.position = "none")

forksb <- forks_plus %>% #mutate(rank = rank/max(rank))%>%
  #dplyr::filter(rank< 0.9)%>%
  group_by(treatment, rank)%>%
  summarize(n= n())%>%
  ggplot()+
  geom_smooth(aes(x = rank, y=n, col= fct_rev(treatment), fill=fct_rev(treatment)), method = 'gam')+
  geom_point(aes(x = rank, y=n, col=fct_rev(treatment)), size=1, alpha = 0.6)+
  ylim(0,4500)+ xlab("S-phase Progression")+ ylab("")+
  theme_bw()+xlab("")+
  theme_bw()+
  scale_fill_manual(values = c("#1F97FF", "#FF7B45" ))+
  scale_color_manual(values = c( "#1F97FF", "#FF7B45"))+ theme(legend.position = "right")+
  facet_wrap(vars(fct_rev(treatment)),ncol = 2)+
  theme(legend.position = "none")


(forksa|forksb)+plot_layout(widths =  c(2,8))
