library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(parallel)
library(readxl)
library(fuzzyjoin)
library(GGally)
library(forcats) 
library(job)
library(magrittr)
library(igraph)
library(tidyr)
library(readr)
library(ggbeeswarm)
library(ggpubr)
library(Biostrings)
#library(edgeR)
#library(tidyverse)
#library(patchwork)
#library(multtest)
#library(stringdist)
#library(rstatix)
#library(ggbeeswarm)
#library(ggnetwork)
#library(devtools)
#library(tcrgrapher)


#metadata----
cumulative_tol_dose<-fread('/data/vskatova/vskatova/single_cell_poised_cd4/POISED_CTD_XH.csv')
cd4_all_inf<-fread('/data/vskatova/vskatova/single_cell_poised_cd4/POISED_peanut_reactive_CD4_scTCR_cluster14_XH.csv') #изначальная таблица

#peanut reacive cells----
cd4_inf<-cd4_all_inf %>% # удобная таблица
  select(V1,seurat_clusters,Group,timepoint,Group2,cluster14,Cell_Group_pea,CTgene,
         CTnt,CTaa,CTstrict,Frequency,cloneType,Cell_Group) %>% 
  separate(col='CTnt',into = c('cdr3nt_tra','cdr3nt_trb'),sep=('_')) %>% 
  separate(col='CTaa',into=c('cdr3aa_tra','cdr3aa_trb'),sep=('_')) %>% 
  separate(col='CTgene',into=c('vj_tra','vj_trb'),sep=('_'),remove =FALSE) %>% 
  mutate_at(c('vj_tra','vj_trb'),
            ~str_replace_all(string=., pattern='\\.',replacement = '_')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove_all(string=.,pattern='\\*\\d+')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove_all(string=.,pattern='\\_TR.C.*')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove(string=.,pattern='None_')) %>% 
  mutate(vj_tra = str_remove_all(vj_tra,'/')) %>% 
  mutate(vj_trb=str_remove_all(vj_trb,'\\TRBD\\d+'))

cd4_inf %<>% mutate_all(., ~ str_replace_na(.,replacement='NA'))  #меняю на NA
cd4_inf %<>%  mutate_all(.,~ str_replace_all(.,pattern='__','_'))

#cd4_inf %>% filter(vj_tra=='NA'& vj_trb=='NA') %>% nrow() # проверить, для скольких сиквенсов нет альфа и бета
#cd4_inf %>% filter(vj_tra=='NA'& vj_trb!='NA') %>% nrow() # проверить, для скольких сиквенсов нет альфа 
#cd4_inf %>% filter(vj_trb=='NA'& vj_tra!='NA') %>% nrow() # проверить, для скольких сиквенсов нет беты 

cd4_inf %<>% filter(!(vj_tra=='NA'& vj_trb=='NA'))# убрать тех реюят, у кого оба сиквенса отсутствуют

##adding F/Q in the cdr3aa----
cd4_seq_corrected_only_allergic<- cd4_inf %>% 
  mutate(cdr3aa_tra=case_when(cdr3aa_tra!='NA'~ paste0('C',cdr3aa_tra,'F'),
                              TRUE~cdr3aa_tra),
         cdr3aa_trb=case_when(cdr3aa_trb!='NA'~ paste0('C',cdr3aa_trb,'F'),
                              TRUE~cdr3aa_trb)) %>% 
  select(-seurat_clusters,-CTgene,-CTstrict,-cloneType) %>% 
  filter(Group!='HC') %>% 
  mutate(Participant_PPID=str_extract(V1,pattern='P\\d+')) %>% #extract names of POISED participant
  filter(Participant_PPID!='NA')

cd4_seq_corrected_only_healthy <- cd4_inf %>% 
  mutate(cdr3aa_tra=case_when(cdr3aa_tra!='NA'~ paste0('C',cdr3aa_tra,'F'),
                              TRUE~cdr3aa_tra),
         cdr3aa_trb=case_when(cdr3aa_trb!='NA'~ paste0('C',cdr3aa_trb,'F'),
                              TRUE~cdr3aa_trb)) %>% 
  select(-seurat_clusters,-CTgene,-CTstrict,-cloneType) %>% 
  filter(Group=='HC') %>% 
  mutate(Participant_PPID=str_extract(V1,pattern='C\\d+')) %>% #extract names of HEALTHY participants
  filter(Participant_PPID!='NA')
  
cd4_seq_corrected<-cd4_seq_corrected_only_allergic %>% 
  bind_rows(cd4_seq_corrected_only_healthy)

##identical cdr3aa in tissue samples----
###tissues----
final_metadata<-fread('/data/vskatova/vskatova/poised_trb/analysis/final_metadata.tsv') #read bulk poised metadata

final_metadata %<>%
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr/',sample))

tissue_samples_bulk_trb<-read_tsv(file=final_metadata %>% #bulk
                                    filter(tissue!='PBMC') %>% 
  filter(Participant_PPID %in% cd4_seq_corrected$Participant_PPID) %>% 
  pull(real_path),id="real_path") %>% 
  dplyr::rename(Count_bulk=cloneCount,
                Freq_bulk=cloneFraction,
                cdr3nt_trb=nSeqCDR3,
                cdr3aa_trb=aaSeqCDR3) %>% 
  mutate(vj_trb=paste0(bestVGene,'_',bestJGene)) |> 
  #filter(Count_bulk>1) %>%  #filter clones with count 1
  group_by(cdr3aa_trb,vj_trb,real_path) %>% # group by cdr3aa and vj to get rid of duplicated because of pcr mistakes 
    summarise(Count_bulk=sum(Count_bulk),
              Freq_bulk=sum(Freq_bulk),
              n_cdr3nt_variants_bulk=n_distinct(cdr3nt_trb)) %>% 
  inner_join(final_metadata %>% select(tissue,real_path, week,Participant_PPID
                      ),by='real_path')

cd4_seq_corrected_only_trb<- cd4_seq_corrected %>% # 
  filter(vj_trb!='NA') %>% # filter ones without trb chain
  dplyr::rename(Count_single_cell=Frequency,
                timepoint_single_cell=timepoint,
                Group_dose_outcome=Group,
                cell_subset=cluster14,
                Cell_type=Cell_Group,
                OIT=Group2) %>% 
  select(Group_dose_outcome,timepoint_single_cell,cell_subset,Count_single_cell,
         vj_trb,cdr3aa_trb,OIT,cdr3nt_trb,Cell_type,Participant_PPID) %>% 
  mutate(Count_single_cell=as.numeric(Count_single_cell)) %>% 
    group_by(Group_dose_outcome,timepoint_single_cell,cell_subset,
             vj_trb,cdr3aa_trb,OIT,Cell_type,Participant_PPID) %>%  #могут быть идентичные aa клоны,но 
     summarise(n_cdr3nt_variants_sc=n_distinct(cdr3nt_trb), #разные по нуклеотидному сиквенсу,объединяем их в один
               Count_single_cell=sum(Count_single_cell))                                                                  

nClones_tissues<-read_tsv(file=final_metadata %>% #bulk для нормализации нужно знать кол-во клонов
           filter(tissue!='PBMC') %>% 
           filter(Participant_PPID %in% cd4_seq_corrected$Participant_PPID) %>% 
           pull(real_path),id="real_path") %>% 
  group_by(real_path) |> 
  summarise(cloneClount=n()) |> 
  inner_join(final_metadata %>% select(tissue,real_path, week,Participant_PPID),by='real_path')


shared_tissue_clns<-inner_join(tissue_samples_bulk_trb, #take bulk
                          cd4_seq_corrected_only_trb,# take antigen spec data
                          by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") #grouping variables

any_tp_sc_clones_in_w0_w52_tissues<-shared_tissue_clns |> 
  mutate(cell_subset=ifelse(cell_subset=='PR_TNa_2','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_3','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_4','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_Tfh13','PR_Th2conv',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Treg_act2','PR_Treg_act1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_act2','PR_Teff_Me_act1',cell_subset)) |> 
  group_by(real_path,tissue,week,Participant_PPID,Group_dose_outcome,cell_subset,OIT) |> 
  summarise(sum_freq_by_tissue_week=sum(Freq_bulk)) |> 
  right_join(nClones_tissues,by=c('real_path','tissue','week','Participant_PPID')) |> 
  mutate(normalized_summed_freq=sum_freq_by_tissue_week/cloneClount) |> 
  mutate(normalized_summed_freq=ifelse(is.na(normalized_summed_freq),0,normalized_summed_freq)) |> 
  filter(week=='W000' |week=='W104') |> 
  filter(OIT=='Active') |> 
  ungroup() |> 
  select(Participant_PPID,tissue,cell_subset,week,normalized_summed_freq) |> 
  complete(Participant_PPID,tissue,cell_subset,week,fill = list(normalized_summed_freq=0)) |> 
 filter(cell_subset=='PR_Th2conv'|cell_subset=='PR_Teff_Me_Th1_CTL' |cell_subset=='PR_Teff_Me_Th2a' |
           cell_subset=='PR_Treg_act1') |> 
  group_by(cell_subset,tissue,Participant_PPID) |> 
  mutate(n_timepoints=n_distinct(week)) |> 
  filter(n_timepoints==2) |> 
  ungroup()
 
  
tissues_n_clones_per_subset<-shared_tissue_clns |> 
  mutate(cell_subset=ifelse(cell_subset=='PR_TNa_2','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_3','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_4','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_Tfh13','PR_Th2conv',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Treg_act2','PR_Treg_act1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_act2','PR_Teff_Me_act1',cell_subset)) |> 
  filter(Participant_PPID %in% any_tp_sc_clones_in_w0_w52_tissues$Participant_PPID) |> 
  group_by(Participant_PPID,cell_subset) |> count()

ay_tp_sc_clones_in_tissues_w0_w52<-tissues_n_clones_per_subset |> 
  filter(cell_subset=='PR_Th2conv'|cell_subset=='PR_Teff_Me_Th1_CTL' |cell_subset=='PR_Teff_Me_Th2a' |
           cell_subset=='PR_Treg_act1') |> 
  left_join(any_tp_sc_clones_in_w0_w52_tissues,by=c('cell_subset','Participant_PPID')) |>
  left_join(nClones_tissues,by=c('tissue','week','Participant_PPID'))
 
ay_tp_sc_clones_in_tissues_w0_w52 |> 
  filter(Participant_PPID!='P118') |> 
  filter( tissue=='stomach' | tissue=='duodenum') |> 
  mutate(tissue=factor(x=tissue,levels=c('stomach','duodenum'))) %>%
  ggplot(aes(x=week,y=normalized_summed_freq,color=tissue))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_point(size=2,alpha=0.5,aes(color=tissue))+
  geom_line(aes(group=Participant_PPID))+
  scale_color_manual(values = c('#e7298a','#19A519'))+
  #scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  facet_wrap(~tissue+cell_subset,scales='free',nrow = 2)+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

ggsave(filename ='/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/all_PR_clones_w0_w104_stomach_duodenum_normalized.png',
       width = 10,height = 6,device = 'png',dpi = 400)

#пересекаются ли они между пациентами?
shared_tissue_clns |> 
  group_by(cdr3aa_trb,vj_trb) |> 
  summarise(Part_label=paste(unique(Participant_PPID),collapse = ','),
            n_patients=n_distinct(unique(Participant_PPID))) |> View()


  group_by(cdr3aa_trb,vj_trb,Participant_PPID,Group_dose_outcome) %>% 
  summarise(weeks_bulk=paste(unique(week),collapse =',' ),
            weeks_sc=paste(unique(timepoint_single_cell),collapse=','),
            #Part_label=paste(unique(Participant_PPID),collapse = ','),#there are no common seq for several patients
            tissues=paste(unique(tissue),collapse =',' ),
            Count_cells=n(),
            Count_sc_2=paste(Count_single_cell,collapse =',' ),
            Count_bulk_sum=sum(Count_bulk),
            cells_subsets=paste(cell_subset,collapse=','))
            #n_pat_labels=n_distinct(Participant_PPID),#there are no such seqs
            #Group_dose_outcome=paste(unique(Group_dose_outcome),collapse =',' ))

make.sequence.graph<-function (.data, .name = '',max_errs=1) {
  G <- graph.empty(n = length(.data), directed=F)
  tmp<-stringdistmatrix(.data,.data,method="hamming") #считает рёбра, возвращает файл типа dist
  print(tmp)
  G <- add.edges(G, t(which(tmp<=max_errs,arr.ind=T))) #добавляет в граф рёбра
  print(G)
  G <- igraph::simplify(G) # To eliminate multiple edges from a graph
  G <- set.vertex.attribute(G, 'label', V(G), .data)
  print(G)
  G
}
###graph---- 
table_for_graph<-shared_tissue_clns |>
  mutate(cdr3aa_trb_vj=paste0(cdr3aa_trb,'_',vj_trb)) |> 
  group_by(cdr3aa_trb_vj) |> 
  summarise(Part_label=paste(unique(Participant_PPID),collapse = ','),
            n_patients=n_distinct(unique(Participant_PPID)),
            tissues=paste(unique(tissue),collapse = ','),
            n_tissues=n_distinct(unique(tissue)))


shared_clones_tissues_graph<- make.sequence.graph(table_for_graph |> 
                                                  pull(cdr3aa_trb_vj),
                                      max_errs = 1)

shared_clones_tissues_graph<-set.vertex.attribute(shared_clones_tissues_graph, 
                                                  'ppid', 
                                                  V(shared_clones_tissues_graph), 
                                                  table_for_graph$Part_label)

shared_clones_tissues_graph<-set.vertex.attribute(shared_clones_tissues_graph,
                                                  'tissues',
                                                  V(shared_clones_tissues_graph), 
                                                  table_for_graph$tissues)

shared_clones_tissues_graph<-set.vertex.attribute(shared_clones_tissues_graph, 
                                                  'n_tissues', 
                                                  V(shared_clones_tissues_graph), 
                                                  table_for_graph$n_tissues)
write_graph(shared_clones_tissues_graph,
            '/data/vskatova/vskatova/poised_trb/analysis/PR_clones_found_in_tissues.gml',format = 'gml')

##to identify subsets----
shared_tissue_clns %<>% 
  mutate(cdr3aa_trb_vj=paste0(cdr3aa_trb,'_',vj_trb))

number_of_cells_that_were_sorted<-cd4_seq_corrected_only_trb %>% 
  filter(Participant_PPID %in% shared_tissue_clns$Participant_PPID) %>% 
  group_by(Participant_PPID,Group_dose_outcome,cell_subset,timepoint_single_cell) %>% 
  count()

mi_dark<-c("#198020",
           "#5F31CC",
           "#C26A27",
           "#068A94",
           "#A324B2",
           "#659406",
           "#105BCC",
           "#AD3757",
           "#5E5E70")
           
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",  "#5F31CC",
                        "#AD3757", "#198020","#A324B2",  "#105BCC" )
#a bar plot           
cells_from_bulk_found_in_sc<-cd4_seq_corrected_only_trb %>% 
  mutate(cdr3aa_trb_vj=paste0(cdr3aa_trb,'_',vj_trb)) %>% 
  filter(Participant_PPID %in% shared_tissue_clns$Participant_PPID) %>%  # отсортировать только тех пациентов, у которых нашлись клоны
  filter(cdr3aa_trb_vj %in% shared_tissue_clns$cdr3aa_trb_vj )

#отсортировать только те клоны, которые есть в тканях
cells_from_bulk_found_in_sc %>% 
  group_by(timepoint_single_cell,Participant_PPID,Group_dose_outcome) %>% 
  mutate(n_cells_one_pat_1timepoint=n()) %>% 
  ungroup() %>% 
  group_by(cell_subset,Participant_PPID,Group_dose_outcome,timepoint_single_cell) %>% 
  summarise(normalized_cell_count=n()/dplyr::first(n_cells_one_pat_1timepoint)) %>%  #может надо поделить на общее количество клеток,которые были найдены
  filter(timepoint_single_cell=='BL' | timepoint_single_cell=='W104') %>%
  #group_by(Participant_PPID,Group_dose_outcome,timepoint_single_cell) %>% 
  #mutate(number_of_cells_norm_new=normalized_number_of_cells/sum(normalized_number_of_cells)) %>% 
  ggplot(aes(x=timepoint_single_cell,y = normalized_cell_count, fill= cell_subset))+
  geom_col()+
  scale_fill_manual(values=cbPalette)+
  facet_wrap(~Participant_PPID+Group_dose_outcome)+
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/clones_that_present_in_tisses_right_normalization.png',
       width = 9,height = 7,device = 'png',dpi = 400)

#a col plot   
cd4_seq_corrected_only_trb %>% 
  mutate(cdr3aa_trb_vj=paste0(cdr3aa_trb,'_',vj_trb)) %>% 
  filter(Participant_PPID %in% shared_tissue_clns$Participant_PPID) %>%  # отсортировать только тех пациентов, у которых нашлись клоны
  filter(cdr3aa_trb_vj %in% shared_tissue_clns$cdr3aa_trb_vj ) %>% #отсортироват только те клоны, которые есть в тканях
  group_by(cell_subset,Participant_PPID,Group_dose_outcome,timepoint_single_cell) %>% 
  count() %>% 
  dplyr::rename(number_of_cells_found_in_tissues=n) %>% #может надо поделить на общее количество клеток,которые были у взяты у этих пациентов
  inner_join(number_of_cells_that_were_sorted %>% 
               dplyr::rename(total_number_of_cells=n),
             by=c('Participant_PPID','Group_dose_outcome','cell_subset','timepoint_single_cell')) %>% 
  mutate(normalized_number_of_cells=number_of_cells_found_in_tissues/total_number_of_cells) %>% 
  filter(timepoint_single_cell=='W104') %>%
  ggplot(aes(x=cell_subset,y=normalized_number_of_cells))+
  geom_col()+
  scale_fill_viridis_d()+
  facet_wrap(~Participant_PPID+Group_dose_outcome)+
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/clones_that_present_in_tisses_w104.png',
       width = 8,height = 4,device = 'png',dpi = 400)


#non-activated cells----
non_activated_cells<-read_excel('/data/vskatova/vskatova/single_cell_poised_cd4/POISED_CD4_scTCR_XH.xlsx')

non_activated_cells %<>%  
  filter(PR3=='CD154-CD137-') #to filter only non-reactive

non_activated_cells_pretty<- non_activated_cells %>% #to make the data pretty
  dplyr::rename(sample_cell_id='...1',
              Count_single_cell=Frequency,
              timepoint_single_cell=timepoint,
              Group_dose_outcome=Group,
              Cell_type=Cell_Group,
              surface_markers=PR3,
              OIT=Group2) %>% 
select(sample_cell_id,Group_dose_outcome,timepoint_single_cell,Cell_type,CTgene,
       CTnt,CTaa,CTstrict,OIT,Count_single_cell,surface_markers) %>% 
  separate(col='CTnt',into = c('cdr3nt_tra','cdr3nt_trb'),sep=('_')) %>% 
  separate(col='CTaa',into=c('cdr3aa_tra','cdr3aa_trb'),sep=('_')) %>% 
  separate(col='CTgene',into=c('vj_tra','vj_trb'),sep=('_'),remove =FALSE) %>% 
  mutate_at(c('vj_tra','vj_trb'),
            ~str_replace_all(string=., pattern='\\.',replacement = '_')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove_all(string=.,pattern='\\*\\d+')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove_all(string=.,pattern='\\_TR.C.*')) %>% 
  mutate_at(c('vj_tra','vj_trb'),~str_remove(string=.,pattern='None_')) %>% 
  mutate(vj_tra = str_remove_all(vj_tra,'/')) %>% 
  mutate(vj_trb=str_remove_all(vj_trb,'\\TRBD\\d+'))

non_activated_cells_pretty %<>% mutate_all(., ~ str_replace_na(.,replacement='NA'))  #change for NA
non_activated_cells_pretty %<>%  mutate_all(.,~ str_replace_all(.,pattern='__','_'))

#non_activated_cells_pretty %>% filter(vj_tra=='NA'& vj_trb=='NA') %>% nrow() # check for alpha and beta absence 
#non_activated_cells_pretty %>% filter(vj_tra=='NA'& vj_trb!='NA') %>% nrow() # check for alpha absence 
#non_activated_cells_pretty %>% filter(vj_trb=='NA'& vj_tra!='NA') %>% nrow() # check for beta absence 

non_activated_cells_pretty %<>% filter(!(vj_tra=='NA'& vj_trb=='NA'))# to filter ones without alpha and beta

##adding F/Q in the cdr3aa----

non_activated_cells_allergic_only<- non_activated_cells_pretty %>% 
  mutate(cdr3aa_tra=case_when(cdr3aa_tra!='NA'~ paste0('C',cdr3aa_tra,'F'),
                              TRUE~cdr3aa_tra),
         cdr3aa_trb=case_when(cdr3aa_trb!='NA'~ paste0('C',cdr3aa_trb,'F'),
                              TRUE~cdr3aa_trb)) %>%  
  filter(Group_dose_outcome!='HC') %>% 
  mutate(Participant_PPID=str_extract(sample_cell_id,pattern='P\\d+')) %>% 
  filter(Participant_PPID!='NA') #filter ones without patient_id

non_activated_cells_healthy_only<- non_activated_cells_pretty %>% 
  mutate(cdr3aa_tra=case_when(cdr3aa_tra!='NA'~ paste0('C',cdr3aa_tra,'F'),
                              TRUE~cdr3aa_tra),
         cdr3aa_trb=case_when(cdr3aa_trb!='NA'~ paste0('C',cdr3aa_trb,'F'),
                              TRUE~cdr3aa_trb)) %>%  
  filter(Group_dose_outcome=='HC') %>% 
  mutate(Participant_PPID=str_extract(sample_cell_id,pattern='C\\d+')) %>% 
  filter(Participant_PPID!='NA')

non_activated_cells_pretty_correct <- non_activated_cells_allergic_only %>% 
  bind_rows(non_activated_cells_healthy_only)

non_activated_cells_cd4_TRB<- non_activated_cells_pretty_correct %>% 
  filter(vj_trb!='NA') %>% #to filter ones without TRB sequence
  mutate(Count_single_cell=as.numeric(Count_single_cell)) %>% 
  group_by(Group_dose_outcome,OIT,timepoint_single_cell,Cell_type, #group by cdr3aa_vj TRB in 1 subset by participant because there might be   
           vj_trb,cdr3aa_trb,Participant_PPID) %>%   # duplication because of different C/D genes
  summarise(n_cdr3nt_variants_sc=n_distinct(cdr3nt_trb),
            n_alpha_cdr3aa_variants=n_distinct(cdr3aa_tra),
            n_alpha_cdr3nt_varians=n_distinct(cdr3nt_tra),
            Count_single_cell=sum(Count_single_cell))               

#identical cdr3aa in tissues
patients_who_have_bulk_tissue_seqs<-tissue_samples_bulk_trb %>% pull(Participant_PPID) %>%
  unique() %>%
  as.data.frame() %>% 
  dplyr::rename(Participant_PPID='.') #to make df with a patient who doesnt have shared clones

shared_tissue_non_reactive<-inner_join(tissue_samples_bulk_trb, #take bulk tissues
                                       non_activated_cells_cd4_TRB,# take non-reactive 
                               by=c('cdr3aa_trb','vj_trb','Participant_PPID'),
                               multiple='all')

shared_tissue_non_reactive_grouped<-shared_tissue_non_reactive %>%  # common - cdr3aa_trb, vj_trb, paticipant id 
  filter(Participant_PPID!='P114') %>% #to exclude a patient with NA
  select(cdr3aa_trb,OIT,vj_trb,Participant_PPID,Group_dose_outcome,week,timepoint_single_cell,
         Cell_type,Count_bulk,Freq_bulk,Count_single_cell,tissue) %>% 
  group_by(cdr3aa_trb,vj_trb) %>% 
  summarise(weeks_bulk=paste(unique(week),collapse =',' ),
            weeks_sc=paste(unique(timepoint_single_cell),collapse=','),
            Part_label=paste(unique(Participant_PPID),collapse = ','),#
            tissues=paste(unique(tissue),collapse =',' ),
            Count_cells=n(),
            Count_sc_2=paste(Count_single_cell,collapse =',' ),
            freq_bulk_sum=sum(Freq_bulk),
            cells_types=paste(Cell_type,collapse=','),
            n_pat_labels=n_distinct(Participant_PPID),#
            Group_dose_outcome=paste(unique(Group_dose_outcome),collapse =',' ))

# non-activated VS activated----
##week0----
w0_non_specific_shared_tissues<-shared_tissue_non_reactive %>% #w0 non activated clones
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  filter(week_bulk=='W000' & timepoint_single_cell=='BL') %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue) %>%
  summarize(n=n()) %>% 
  mutate(PE_sort='non_reactive',
         week='W000') %>% 
  inner_join(final_metadata %>% select(tissue,real_path,Participant_PPID,week), by=c('tissue','Participant_PPID','week'))
 
w0_specific_shared_tissues<-inner_join(tissue_samples_bulk_trb, #take bulk
                               cd4_seq_corrected_only_trb, #take antigen specific clones
                               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  filter(week_bulk=='W000' & timepoint_single_cell=='BL') %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='reactive',
         week='W000') %>% 
  inner_join(final_metadata %>% select(tissue,real_path,Participant_PPID,week), by=c('tissue','Participant_PPID','week'))

w0_reactive_non_reactive_shared_with_tissues<-w0_non_specific_shared_tissues %>% 
  bind_rows(w0_specific_shared_tissues)

#to add inf about clone counts in bulk and in single cell 

counts_bulk<-read_tsv(final_metadata %>% 
                        filter(week=='W000' | week=='W104' | week=="W117" ,
                               tissue!='PBMC') %>% 
                        filter(Participant_PPID %in% patients_who_have_bulk_tissue_seqs$Participant_PPID ) %>%
                        select(Participant_PPID,week,tissue,real_path) %>% 
                        pull(real_path),id = 'real_path') %>% #bulk counts
  group_by(real_path) %>% 
  summarise(n_Clones=n(),
            n_Reads=sum(cloneCount))

sc_specific_counts<-cd4_seq_corrected_only_trb %>%  #all specific cs
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',#there are different week 117 samples
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(timepoint_single_cell,OIT,cell_population,Participant_PPID) %>% #identical clones might be in different subclusters
  summarize(n_clones_sc=n()) %>% 
  mutate(PE_sort='reactive')

sc_non_specific_counts<-non_activated_cells_cd4_TRB %>% #all non specific sc 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(timepoint_single_cell,OIT,cell_population,Participant_PPID) %>% #identical clones might be in different big clusters
  summarize(n_clones_sc=n()) %>% 
  mutate(PE_sort='non_reactive')

sc_counts<-sc_specific_counts %>% bind_rows(sc_non_specific_counts) %>% 
  filter(Participant_PPID %in% patients_who_have_bulk_tissue_seqs$Participant_PPID)

## inf number of cells and clones in bulk rep and to add zero values----
w0_reactive_non_reactive_cell_numbers<- w0_reactive_non_reactive_shared_with_tissues %>% 
  select(-real_path) %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>% #to add zero values 
   ungroup() %>%
 pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','week','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>% 
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,week,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>%  #magic
  mutate(OIT=ifelse(is.na(OIT),'Placebo',OIT)) %>% 
  inner_join(sc_counts %>% 
               dplyr::rename(week=timepoint_single_cell) %>% 
               mutate(week=ifelse(week=='BL','W000',week)) %>% 
               filter(week=='W000'),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path'),
             by=c('Participant_PPID','tissue','week')) %>% 
  select(-real_path) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_Clones=as.double(n_Clones),
         normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))

#graph
w0_reactive_non_reactive_cell_numbers %>% 
  mutate(tissue=factor(x=tissue,
                           levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/reactive_vs_non_reactive_tissue_w0.png',
       width = 12,height = 6,device = 'png',dpi = 400)


##week104----
w104_non_specific<-shared_tissue_non_reactive %>% #w104 non activated clones
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  filter(week_bulk=='W104' & timepoint_single_cell=='W104') %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue) %>% 
  summarize(n=n()) %>% 
  mutate(PE_sort='non_reactive')

w104_shared_tissue_reactive<-inner_join(tissue_samples_bulk_trb, 
                                   cd4_seq_corrected_only_trb,
                                   by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  filter(week_bulk=='W104' & timepoint_single_cell=='W104') %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue) %>% 
  summarize(n=n()) %>% 
  mutate(PE_sort='reactive')

w104_total_shared_with_tissues<-w104_non_specific %>% 
  bind_rows(w104_shared_tissue_reactive) %>% 
  mutate(week='W104')

#normalization

## inf number of cells and clones in bulk rep and to add zero values----
w104_reactive_non_reactive_cell_numbers<- w104_total_shared_with_tissues %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>%  #to add zero values 
  ungroup() %>%
  pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','week','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>%
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,week,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>%  #magic
  mutate(OIT=case_when(Participant_PPID=='P101'~'Active',
                       Participant_PPID=='P105'~'Active',
                       Participant_PPID=='P108'~'Active',
                       Participant_PPID=='P113'~'Active',
                       Participant_PPID=='P114'~'Placebo',
                       Participant_PPID=='P118'~'Placebo',
                       TRUE~OIT)) %>%
  inner_join(sc_counts %>% 
               dplyr::rename(week=timepoint_single_cell) %>% 
               filter(week=='W104'),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path'),
             by=c('Participant_PPID','tissue','week')) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_Clones=as.double(n_Clones),
         normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))

w104_reactive_non_reactive_cell_numbers %>% 
  filter(OIT!='Placebo') %>% #to exclude placebo patients (2) and for one of them there is no w104 samples
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/reactive_vs_non_reactive_tissue_w104_without_placebo.png',
       width = 12,height = 6,device = 'png',dpi = 400)


##week 117----
w117_shared_non_reactive<-shared_tissue_non_reactive %>% #to select all clones from W117 in single cell and any point in BULK
  filter(timepoint_single_cell=='W117' | timepoint_single_cell=='W117B' |
           timepoint_single_cell=='W117B_') %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,timepoint_single_cell) %>% 
  summarize(n=n()) %>% 
  mutate(PE_sort='non_reactive')

w117_shared_tissue_reactive<-inner_join(tissue_samples_bulk_trb, 
                                        cd4_seq_corrected_only_trb,
                                        by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  filter(timepoint_single_cell=='W117') %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,timepoint_single_cell) %>% 
  summarize(n=n()) %>% 
  mutate(PE_sort='reactive')

w117_total_shared_with_tissues_reactive_non_reactive<-w117_shared_non_reactive %>% 
  bind_rows(w117_shared_tissue_reactive)

#normalization

## inf number of cells and clones in bulk rep and to add zero values----
w117_react_non_react_cell_numbers<- w117_total_shared_with_tissues_reactive_non_reactive %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>%  #to add zero values 
  ungroup() %>%
  pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','timepoint_single_cell','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>% 
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,timepoint_single_cell,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>% #magic
  mutate(OIT=case_when(Participant_PPID=='P114'~'Placebo',
                       Participant_PPID=='P118'~'Placebo',
                       TRUE~OIT))

w117_react_non_react_cell_counts<-w117_react_non_react_cell_numbers %>% 
  dplyr::rename(week=timepoint_single_cell) %>% 
  inner_join(sc_counts %>% 
               dplyr::rename(week=timepoint_single_cell) %>% 
               mutate(week=case_when(week=='W117B'~'W117',
                                     week=='W117B_'~'W117',
                                    TRUE~week)) %>% 
               filter(week=='W117'),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path'),
             by=c('Participant_PPID','tissue','week')) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_Clones=as.double(n_Clones),
         normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))

w117_react_non_react_cell_counts %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

  ggsave(filename = '/home/lera/poised_trb/analysis/reactive_vs_non_reactive_tissue_w117_only_IT.png',
       width = 12,height = 6,device = 'png',dpi = 400)

##to unite all the single-cell data and exclude identical cdr3aa between specific and non-specific----
#to unite peanut-reactive data
united_sc_reactive_by_patient_id<-cd4_seq_corrected_only_trb %>% 
    ungroup() %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117', #not sure
                                         timepoint_single_cell=='W117B_'~'W117',#not sure
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(Group_dose_outcome,vj_trb,cdr3aa_trb,OIT,Participant_PPID,cell_population) %>% 
    summarise(n_timepoints_sc=n_distinct(timepoint_single_cell),
              timepoints_sc=paste(unique(timepoint_single_cell),collapse=','),
              Count_single_cell_sum=sum(Count_single_cell))

#to unite non-specific data
united_sc_non_reactive_by_patients_id<-non_activated_cells_cd4_TRB %>% 
  ungroup() %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117', #not sure
                                         timepoint_single_cell=='W117B_'~'W117',#not sure
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(Group_dose_outcome,vj_trb,cdr3aa_trb,OIT,Participant_PPID,cell_population) %>% 
  summarise(n_timepoints_sc=n_distinct(timepoint_single_cell),
            timepoints_sc=paste(unique(timepoint_single_cell),collapse=','),
            Count_single_cell_sum=sum(Count_single_cell))

cells_that_are_react_andnot_reactive<-inner_join(united_sc_reactive_by_patient_id,united_sc_non_reactive_by_patients_id,
          by=c('vj_trb','cdr3aa_trb','Group_dose_outcome','OIT','Participant_PPID','cell_population'),suffix=c('.reactive','.non_reactive')) %>%
  mutate(FC_count=Count_single_cell_sum.reactive/Count_single_cell_sum.non_reactive) %>% 
    filter(FC_count<10)

united_sc_reactive_clean<-united_sc_reactive_by_patient_id %>% #to exclude from reactive
  anti_join(cells_that_are_react_andnot_reactive,
            by=c('vj_trb','cdr3aa_trb','Group_dose_outcome','OIT','Participant_PPID','cell_population'))

united_sc_non_reactive_clean<-united_sc_non_reactive_by_patients_id %>% #to exclude from non reactive
  anti_join(cells_that_are_react_andnot_reactive,
            by=c('vj_trb','cdr3aa_trb','Group_dose_outcome','OIT','Participant_PPID','cell_population'))

###to use all clones that were assigned as peanut specific without filter----

####comparison week0----
shared_w0_reactive<-united_sc_reactive_clean %>% #to calculate reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W000'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='reactive')

shared_w0_non_reactive<-united_sc_non_reactive_clean %>% #to calculate non reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W000'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='non_reactive')

shared_tissues_w0_reactive_non_reactive<-shared_w0_reactive %>% #to bind reactive and non reactive
  bind_rows(shared_w0_non_reactive)

shared_tissues_w0_react_non_react_counts<- shared_tissues_w0_reactive_non_reactive %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>%  #to add zero values 
  ungroup() %>% 
  pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','week','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>% 
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,week,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>% #magic
  mutate(OIT=case_when(Participant_PPID=='P114'~'Placebo',
                       Participant_PPID=='P105'~'Active',
                       TRUE~OIT)) %>% 
    inner_join(sc_counts %>% 
                 dplyr::rename(week=timepoint_single_cell) %>% 
                 mutate(week=ifelse(week=='BL','W000',week)) %>% 
                 filter(week=='W000'),
               by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
    inner_join(counts_bulk %>% 
                 inner_join(final_metadata %>% select(Participant_PPID,
                                                      tissue,week,real_path),by='real_path'),
               by=c('Participant_PPID','tissue','week')) %>% 
    select(-real_path) %>% 
    mutate(sc_counts=as.double(n_clones_sc),
           n_Clones=as.double(n_Clones),
           normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))
#graph
shared_tissues_w0_react_non_react_counts %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/reactive_vs_non_reactive_tissue_w0_version2.png',
       width = 12,height = 6,device = 'png',dpi = 400)

####comparison week104----

shared_w104_reactive<-united_sc_reactive_clean %>% #to calculate reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W104'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='reactive')

shared_w104_non_reactive<-united_sc_non_reactive_clean %>% #to calculate non reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W104'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='non_reactive')

shared_tissues_w104_reactive_non_reactive<-shared_w104_reactive %>% #to bind reactive and non reactive
  bind_rows(shared_w104_non_reactive)

shared_tissues_w104_react_non_react_counts<- shared_tissues_w104_reactive_non_reactive %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>%  #to add zero values 
  ungroup() %>% 
  pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','week','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>% 
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,week,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>%  #magic
mutate(OIT=case_when(Participant_PPID=='P101'~'Active',
                     Participant_PPID=='P105'~'Active',
                     Participant_PPID=='P108'~'Active',
                     Participant_PPID=='P113'~'Active',
                     Participant_PPID=='P114'~'Placebo',
                     Participant_PPID=='P118'~'Placebo',
                     TRUE~OIT)) %>%
  inner_join(sc_counts %>% 
               dplyr::rename(week=timepoint_single_cell) %>% 
               filter(week=='W104'),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path'),
             by=c('Participant_PPID','tissue','week')) %>% 
  select(-real_path) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_Clones=as.double(n_Clones),
         normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))

shared_tissues_w104_react_non_react_counts %>% 
  filter(OIT!='Placebo') %>% #to exclude placebo patients (2) and for one of them there is no w104 samples
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/home/lerucha38/POISED_single_cell_data_Xiaorui_Han/analysis/reactive_vs_non_reactive_tissue_w104_without_placebo_version2.png',
       width = 12,height = 6,device = 'png',dpi = 400)


####comparison week117----

shared_w117_reactive<-united_sc_reactive_clean %>% #to calculate reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W117'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='reactive')

shared_w117_non_reactive<-united_sc_non_reactive_clean %>% #to calculate non reactive cells
  inner_join(tissue_samples_bulk_trb %>% 
               filter(week=='W117'),
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,tissue,week)  %>% 
  summarise(n=n()) %>% 
  mutate(PE_sort='non_reactive')

shared_tissues_w117_reactive_non_reactive<-shared_w117_reactive %>% #to bind reactive and non reactive
  bind_rows(shared_w117_non_reactive)

shared_tissues_w117_react_non_react_counts<- shared_tissues_w117_reactive_non_reactive %>% 
  pivot_wider(names_from= tissue, values_from = n,values_fill = 0) %>%  #to add zero values 
  ungroup() %>% 
  pivot_longer(cols=!c('Participant_PPID','OIT','cell_population','week','PE_sort'),
               names_to = 'tissue',values_to = 'n_found_shared_clones') %>% 
  mutate(Participant_PPID=factor(Participant_PPID,levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% #magic
  complete(Participant_PPID,week,cell_population,PE_sort,tissue,fill = list(n_found_shared_clones=0)) %>%  #magic
  mutate(OIT=case_when(Participant_PPID=='P105'~'Active',
                       Participant_PPID=='P101'~'Active',
                       Participant_PPID=='P114'~'Placebo',
                       Participant_PPID=='P118'~'Placebo',
                       TRUE~OIT)) %>%
  inner_join(sc_counts %>% 
               dplyr::rename(week=timepoint_single_cell) %>% 
               filter(week=='W117'),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>%
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path'),
             by=c('Participant_PPID','tissue','week')) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_Clones=as.double(n_Clones),
         normalized_n_shr_clones=n_found_shared_clones/(sc_counts*n_Clones))

shared_tissues_w117_react_non_react_counts %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/reactive_vs_non_reactive_tissue_w117_without_placebo_version2.png',
       width = 12,height = 6,device = 'png',dpi = 400)

#activated w0 VS w104 VS w117-----

shared_tissues_w0_react_non_react_counts %>% 
  bind_rows(shared_tissues_w104_react_non_react_counts %>%
              filter(OIT!='Placebo')) %>% 
  bind_rows(shared_tissues_w117_react_non_react_counts %>% 
              filter(OIT!='Placebo')) %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones,fill=week))+
  scale_fill_viridis_d(alpha=0.4)+
  geom_boxplot(outlier.shape = NA,position = position_dodge(width = 0.8))+
  geom_point(position = position_dodge(width = 0.8))+
  stat_compare_means(label = 'p.format' )+
  facet_grid(cell_population~tissue,scales = 'free_y')

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/react_vs_non_react_onlyOIT_by_week.png',
       width = 14,height = 8,device = 'png',dpi = 400)

#to unite tissues reactive VS non-reactive----

by_patient_all_tissues_reactive<-united_sc_reactive_clean %>% 
  inner_join(tissue_samples_bulk_trb,
            by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,week)  %>% 
  summarise(n=n()) %>% 
  filter(week!='W052') %>% #because i cannot normalize this number of clones (there is no w52 single-cell data)
  mutate(PE_sort='reactive')

by_patient_all_tissues_non_reactive<-united_sc_non_reactive_clean %>% #to calculate non reactive cells
  inner_join(tissue_samples_bulk_trb ,
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  group_by(Participant_PPID,OIT,cell_population,week)  %>% 
  summarise(n=n()) %>% 
  filter(week!='W052') %>% #because i cannot normalize this number of clones
  mutate(PE_sort='non_reactive') 

tissues_react_non_react<-by_patient_all_tissues_reactive %>%
  bind_rows(by_patient_all_tissues_non_reactive) %>% 
  ungroup() %>% 
  mutate(Participant_PPID=factor(Participant_PPID,
                                 levels=sc_counts %>% pull(Participant_PPID) %>% unique())) %>% 
  complete(Participant_PPID,week,cell_population,
           PE_sort,fill = list(n=0)) %>%
  mutate(OIT=case_when(Participant_PPID=='P105'~'Active',
                       Participant_PPID=='P108'~'Active',
                       Participant_PPID=='P101'~'Active',
                       Participant_PPID=='P114'~'Placebo',
                       Participant_PPID=='P118'~'Placebo',
                       TRUE~OIT)) %>% 
  inner_join(sc_counts %>% 
               mutate(timepoint_single_cell=ifelse(timepoint_single_cell=='BL','W000',
                                                   timepoint_single_cell)) %>% 
               dplyr::rename(week=timepoint_single_cell),
             by=c('cell_population','Participant_PPID','week','PE_sort','OIT')) %>% 
  inner_join(counts_bulk %>% 
               inner_join(final_metadata %>% select(Participant_PPID,
                                                    tissue,week,real_path),by='real_path') %>% 
               group_by(Participant_PPID,week) %>% 
               summarize(n_clones_meean_per_tissue=mean(n_Clones)),
             by=c('Participant_PPID','week')) %>% 
  mutate(sc_counts=as.double(n_clones_sc),
         n_clones_meean_per_tissue=as.double(n_clones_meean_per_tissue),
         normalized_n_shr_clones=n/(sc_counts*n_clones_meean_per_tissue))
  
tissues_react_non_react %>% 
  filter(!(week=='W104' & OIT=='Placebo')) %>% #to exclude placebo patients from week 104 and w117
  filter(!(week=='W117' & OIT=='Placebo')) %>% 
  mutate(week=factor(x=week,
                       levels=c('W000','W104','W117'))) %>% 
  ggplot(aes(x=PE_sort,y=normalized_n_shr_clones))+
  scale_fill_viridis_d(alpha=0.4)+
  geom_boxplot(outlier.shape = NA,position = position_dodge(width = 0.8))+
  geom_quasirandom(position = position_dodge(width = 0.8))+
  stat_compare_means(label = 'p.format' )+
  facet_grid(cell_population~week,scales = 'free_y')

ggsave(filename = '/home/lera/single_cell_poised_cd4/analysis/react_vs_non_react_united_tissues_onlyOIT.png',
       width = 8,height = 6,device = 'png',dpi = 400)


#clone-wise comp----
#rna molecules count
sc_reactive_rna_inf<-cd4_seq_corrected_only_trb %>% #reactive
  ungroup() %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',#there are different week 117 samples
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(timepoint_single_cell,OIT,cell_population,Participant_PPID) %>%
  summarise(n_rna_molecules=sum(Count_single_cell)) %>% 
  mutate(PE_sort='reactive')

sc_non_reactive_rna_inf<-non_activated_cells_cd4_TRB %>% #non-reactive
  ungroup() %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  group_by(timepoint_single_cell,OIT,cell_population,Participant_PPID) %>%
  summarise(n_rna_molecules=sum(Count_single_cell)) %>% 
  mutate(PE_sort='non_reactive')

united_sc_rna_count<-sc_reactive_rna_inf %>% bind_rows(sc_non_reactive_rna_inf) #united

#for this clone-wise comp i need to use non-united data,but non-united data
#has clones that are in reactive and non-reactive
#so i need to exlude them again

shared_tissues_reactive<-inner_join(tissue_samples_bulk_trb, #take bulk
           cd4_seq_corrected_only_trb, #take antigen specific clones
           by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(PE_sort='reactive')

shared_tissue_non_reactive %<>% 
dplyr::rename(week_bulk=week) %>% 
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                                   Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                                   Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  mutate(PE_sort='non_reactive')

cells_that_are_react_and_not_react<-inner_join(shared_tissues_reactive,shared_tissue_non_reactive,
                                                 by=c('vj_trb','cdr3aa_trb',
                                                      'Group_dose_outcome','OIT','Participant_PPID',
                                                      'cell_population'),
                                                 suffix=c('.reactive','.non_reactive'),
                                               multiple = "all") %>%
  mutate(FC_count=Count_single_cell.reactive/Count_single_cell.non_reactive) %>% 
  filter(FC_count<10)

sc_shared_with_tissues_react_clean<-shared_tissues_reactive %>% #to exclude from reactive
  anti_join(cells_that_are_react_and_not_react,
            by=c('vj_trb','cdr3aa_trb',
                 'Group_dose_outcome','OIT','Participant_PPID','cell_population'))

sc_shared_with_tissues_non_react_clean<-shared_tissue_non_reactive %>% #to exclude from non reactive
  anti_join(cells_that_are_react_and_not_react,
            by=c('vj_trb','cdr3aa_trb','Group_dose_outcome',
                 'OIT','Participant_PPID','cell_population'))

sc_bulk_shared_normalized<-sc_shared_with_tissues_react_clean %>% 
  bind_rows(sc_shared_with_tissues_non_react_clean) %>% #to unite datasets
  select(-n_alpha_cdr3aa_variants,- n_alpha_cdr3nt_varians) %>% 
  inner_join(counts_bulk %>% select(real_path,n_Reads),by='real_path') %>% #to add inf about reads
  #filter(week_bulk=='W000') %>% #to filter week0 to try
  inner_join(united_sc_rna_count %>% ungroup(),
             by=c('PE_sort','Participant_PPID',
                  'cell_population','OIT','timepoint_single_cell')) %>% 
  dplyr::rename(n_Reads_bulk=n_Reads) %>% 
  mutate(n_Reads_bulk=as.double(n_Reads_bulk),
         Count_single_cell=as.double(Count_single_cell),
         n_rna_molecules=as.double(n_rna_molecules),
         Count_bulk=as.double(Count_bulk),
         n_normalized=Count_single_cell/(n_Reads_bulk*n_rna_molecules)) %>% 
  inner_join(final_metadata %>% 
               select(Group3,Participant_PPID),by='Participant_PPID',multiple = "all") %>% 
  unique()
  
sc_bulk_shared_normalized %>% 
  filter(week_bulk=='W104') %>% 
  filter(Group3!='Placebo') %>%
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=PE_sort,y=n_normalized))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_grid(cell_population~tissue)

ggsave(filename = '/data/vskatova/vskatova/single_cell_poised_cd4/analysis/clonewise_reac_vs_nonreact_W117.png',
       width = 12,height = 8,device = 'png',dpi = 400)

#phenotype change across all patients----
#to exclude clones that are n=in reactive and non-reactive data

non_activated_cells_cd4_TRB %<>% ungroup() %>% #to add important inf about cell type
  mutate(cell_population=case_when(Cell_type=='NR_CD154-CD137-_Teff_Me'~'T_eff',
                          Cell_type=='NR_CD154-CD137-_Na'~'T_eff',
                          Cell_type=='NR_CD154-CD137-_Treg'~'T_reg',
                          Cell_type=='NR_differentiating Treg cells'~'T_reg')) %>% 
  mutate(PE_sort='non_reactive')

cd4_seq_corrected_only_trb %<>% #reactive
  ungroup() %>% 
  mutate(cell_population=case_when(Cell_type=='PR_CD154+_Teff_Me'~'T_eff',
                                   Cell_type=='PR_CD154+Teff_Na'~'T_eff',
                                   Cell_type=='PR_CD137+CD154-_Treg'~'T_reg')) %>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',#there are different week 117 samples
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell)) %>% 
  mutate(PE_sort='reactive')

reac_non_react_clones<-inner_join(cd4_seq_corrected_only_trb,non_activated_cells_cd4_TRB,
           by=c('Participant_PPID','OIT','vj_trb','cdr3aa_trb','cell_population','Group_dose_outcome'),
           multiple='all')

sc_react_clean<-cd4_seq_corrected_only_trb %>% #to exclude from reactive
  anti_join(reac_non_react_clones,
            by=c('vj_trb','cdr3aa_trb',
                 'Group_dose_outcome','OIT','Participant_PPID','cell_population'))

sc_non_react_clean<-non_activated_cells_cd4_TRB %>% #to exclude from non reactive
  anti_join(reac_non_react_clones,
            by=c('vj_trb','cdr3aa_trb','Group_dose_outcome',
                 'OIT','Participant_PPID','cell_population'))

clones_multiple_subsets<-sc_react_clean %>% 
  select(-Cell_type,-n_cdr3nt_variants_sc) %>% 
    group_by(cdr3aa_trb,vj_trb,PE_sort,OIT,Participant_PPID,Group_dose_outcome) %>% 
    summarise(phenotypes=paste(unique(cell_subset),collapse = ','),
              numb_phenotypes=n_distinct(cell_subset),
              numb_tp_sc=n_distinct(timepoint_single_cell),
              tp_sc=paste(unique(timepoint_single_cell),collapse = ',')) %>% 
  filter(numb_tp_sc>1,
         numb_phenotypes>1)
  
 test_grep<- clones_multiple_subsets[(grepl("Treg", clones_multiple_subsets$phenotypes) |
                grepl("Th1Th17", clones_multiple_subsets$phenotypes) |
                   grepl("Th2a", clones_multiple_subsets$phenotypes)),]


  
#low-coverage trb data----

bd_116f<-fread('/home/lera/single_cell_poised_cd4/low_coverage_data/mixcr/bd_Tcells_sc_116f.clones.tsv')
bd_116f %>% filter(grepl('TRB',allVHitsWithScore)) %>% View('bd_116f_trb')
  
bd_116g<-fread('/home/lera/single_cell_poised_cd4/low_coverage_data/mixcr/bd_Tcells_sc_116g.clones.tsv')
bd_116g %>% filter(grepl('TRB',allVHitsWithScore)) %>% View('bd_116g_trb')

#healthy patients sc----
sc_patients_id<-non_activated_cells_cd4_TRB %>% pull(Participant_PPID) %>% unique()

#PBMC----

pbmc_samples<-read_tsv(file=final_metadata %>% #bulk
                         filter(tissue=='PBMC') %>% 
                         filter(Participant_PPID %in% cd4_seq_corrected$Participant_PPID) %>% 
                         pull(real_path),id="real_path") %>% 
  dplyr::rename(Count_bulk=cloneCount,
                Freq_bulk=cloneFraction,
                cdr3nt_trb=nSeqCDR3,
                cdr3aa_trb=aaSeqCDR3) %>% 
  mutate(vj_trb=paste0(bestVGene,'_',bestJGene)) %>% 
  filter(Count_bulk>1) %>%  #filter clones with count 1
  group_by(cdr3aa_trb,vj_trb,real_path) %>% # group by cdr3aa and vj to get rid of duplicated because of pcr mistakes 
  summarise(Count_bulk=sum(Count_bulk),
            Freq_bulk=sum(Freq_bulk),
            n_cdr3nt_variants_bulk=n_distinct(cdr3nt_trb)) %>% 
  inner_join(final_metadata %>% select(tissue,real_path, week,Participant_PPID),by='real_path')

shared_pbmc_activated<-inner_join(cd4_seq_corrected_only_trb,pbmc_samples,
                                  by=c('Participant_PPID','cdr3aa_trb','vj_trb'),multiple = "all")

#public data----

public_data_files<-fread('/data/vskatova/vskatova/single_cell_poised_cd4/analysis/E-MTAB-9532.sdrf.txt')

public_data_files %<>% 
  subset(.,select=which(!duplicated(colnames(.)))) %>% #to remove duplicated colnames
  dplyr::rename(cellType='Characteristics[cell type]',
                age='Characteristics[age]',
                sex='Characteristics[sex]',
                tissue='Characteristics[organism part]',
                age_group='Characteristics[developmental stage]',
                part_id='Characteristics[individual]',
                FASTQ_URI ='Comment[FASTQ_URI]',
                organism='Characteristics[organism]',
                ENA_sample_id='Comment[ENA_SAMPLE]',
                read1_file='Comment[read1 file]',
                read2_file='Comment[read2 file]',
                protocol='Factor Value[protocol]') %>% 
  filter(cellType=='T cell') %>% 
  select(cellType,age,sex,tissue,age_group,FASTQ_URI,organism,ENA_sample_id,read1_file,read2_file,protocol) %>% 
  mutate(filename_short=str_extract(read1_file,'\\Human.*_S')) %>% 
  mutate(filename_short=str_remove(filename_short,'_S'))
  
write_tsv('/data/vskatova/vskatova/single_cell_poised_cd4/analysis/metadata_puplic_sc_tcells.txt',col_names = FALSE)

report_sc<-fread('/data/vskatova/vskatova/single_cell_poised_cd4/public_sc_cd4/mixcr/report_10x_vdj_public.tsv') %>%  #mixcr report 
  mutate(filename_short=str_extract(fileName,'\\Human.*_S')) %>% 
  mutate(filename_short=str_remove(filename_short,'_S')) %>% 
  inner_join(public_data_files %>% 
               select(tissue,age,sex,protocol,filename_short),by='filename_short') %>% 
  unique() %>% 
  select(-droppedLowQual,-droppedLowQualPercents,-alignmentFailedNoBarcode,-alignmentFailedNoBarcodePercents,
         -alignmentFailedSampleNotMatched,-alignmentFailedSampleNotMatchedPercents,
         -alignmentFailedVAndJOnDifferentTargetsPercents,
         -alignmentFailedVAndJOnDifferentTargets,
         -alignmentFailedNoVHits,-alignmentFailedNoCDR3Parts,-alignmentFailedNoCDR3PartsPercents,-alignmentFailedNoVHitsPercents,
         -alignmentFailedNoJHitsPercents,-alignmentFailedNoJHits,-droppedNoClonalSeq,
         -droppedFailedMapping)

public_cd4_GI_clonsents<-read_tsv(public_cd4_GI_files,id = 'fileName') %>% #to take clonsets
  mutate(fileName=str_extract(fileName,'\\Human.*_S')) %>% 
  mutate(fileName=str_remove(fileName,'_S')) %>% 
  select(-targetQualities,
         -allDHitsWithScore,
         -allVAlignments,
         -allDAlignments,
         -allJAlignments,
         -allCAlignments,
         -minQualFR1,
         -minQualCDR1,
         -minQualCDR2,
         -minQualFR3,
         -minQualCDR3,
         -minQualFR4,
         -refPoints)
public_cd4_GI_clonsents %>% 
  mutate(allVHitsWithScore=str_remove(allVHitsWithScore,'\\*\\d+.*')) %>% 
  mutate(allJHitsWithScore=str_remove(allJHitsWithScore,'\\*\\d+.*')) %>% 
  select(fileName,tagValueCELL,readCount,uniqueMoleculeCount,allVHitsWithScore,allJHitsWithScore,allCHitsWithScore,
         nSeqCDR3,aaSeqCDR3) %>% 
  mutate(v_tra=str_extract(allVHitsWithScore,'\\TRA.*')) %>% 
  mutate(v_trb=str_extract(allVHitsWithScore,'\\TRB.*')) %>% 
  select(-allVHitsWithScore) %>% 
  pivot_wider(names_from = 'tagValueCELL' ) %>% glimpse()


final_metadata %>% 
  filter(tissue!='PBMC') %>% nrow()

