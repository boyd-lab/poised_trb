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
library(stringdist)
library(ggbeeswarm)
library(ggpubr)
library(Biostrings)
library(viridis)
library(ggcorrplot)
library(stats)
library(rworldmap)
library(ggfortify)
library(edgeR)

#metadata poised----
metadata_trb_penaut<-fread('/h/lera/poised_trb/PoisedSpecimensCellandRNAInputs.csv')#метадата от Рамоны в моей папкке

meta_trb_penaut<-metadata_trb_penaut %>%
  mutate(specimen_tissue=case_when(specimen_tissue=='Distal esophagus'~'dist_esophagus',
                                   specimen_tissue=='duodenum and stomach'~'duodenum_stomach',
                                   specimen_tissue=='Middle esophagus'~'mid_esophagus',
                                   specimen_tissue=='Proximal esophagus'~'prox_esophagus',
                                   TRUE~specimen_tissue),
         timepoint=ifelse(timepoint=='WK104','W104',timepoint)) %>% 
  mutate_at(c('participant_label','sample_label','specimen_label'),
              ~str_replace_all(string=., pattern='-',replacement = '_'))
  
new_meta<-fread('/home/lera/poised_trb/PoisedSpecimensMetadataFinal.txt') %>% 
  select(specimen_label,participant_label,specimen_tissue,
         specimen_time_point,timepoint,specimen_description,Before.initial.FC,
         Group.x,Wk104.Challenge.Result) %>% 
  dplyr::rename(week104_FC_result=Wk104.Challenge.Result,
                FC=Before.initial.FC,
                group=Group.x) %>% 
  filter(!is.na(week104_FC_result)) %>% 
  mutate(week104_FC_result=ifelse(week104_FC_result=='NO CHALLENGE','NO_CHALLENGE', week104_FC_result)) %>% 
  mutate_at(c('participant_label','specimen_label'),
            ~str_replace_all(string=., pattern='-',replacement = '_'))
  
meta_3<-fread('/home/lera/poised_trb/POISED_DBPCFC.txt') %>% 
  mutate(participant_label=
           str_replace_all(string=participant_label, pattern='-',replacement = '_')) %>% 
  inner_join(new_meta,by='participant_label')

##reparsing----
filenames<-dir('/data/vskatova/vskatova/poised_trb/poised_part_tables/')

list_names_test<-mclapply(filenames,function(filename){
  table<-fread(paste0('/data/vskatova/vskatova/poised_trb/poised_part_tables/',filename))
  pat_labels<-unique(table$participant_label)
  print(pat_labels)
  lapply(pat_labels,function(pat_label){
  
    table_1<-table %>% dplyr::filter(participant_label==pat_label)
    specimen_labels<-unique(table_1$specimen_label)
  
    lapply(specimen_labels,function(sample_lbl){
  
           table_2<-table_1 %>% dplyr::filter(specimen_label==sample_lbl)
            x<-DNAStringSet(table_2 %>% pull(trimmed_sequence))
            names(x)<- table_2 %>%  pull(trimmed_read_id)  
            raw_name<-paste0(sample_lbl,'_',pat_label) 
            full_name<-str_replace_all(raw_name,'-','_')
            file_path<-paste0('/data/vskatova/vskatova/poised_trb/poised_fasta/',full_name, '.fasta')
            print(file_path)
            Biostrings::writeXStringSet(x,filepath=file_path)
           }) }) },mc.cores = 60,mc.cleanup = TRUE)

final_metadata<-fread('/data/vskatova/vskatova/poised_trb/analysis/final_metadata.tsv')

#depth analysis by MiSeq runs---- 
#(to get stat what runs should be resequenced)
final_metadata %<>% 
  mutate(miseq_run=str_extract(sample,'\\M\\d+'))

nReads_nCLones_allFiles<-lapply(final_metadata %>% 
         mutate(path='/home/lera/poised_trb/poised_clonsets/mixcr/') %>% 
         mutate(path=paste0(path,sample)) %>% 
         pull(path),function(file_name){
           fread(file_name) %>% 
             mutate(sample_id=file_name) %>%
             mutate(sample_id=str_remove_all(sample_id,"\\/.*\\/") ) %>% 
             group_by(sample_id) %>% 
             summarise(nReads=sum(cloneCount),
                       nClones=n()) %>% 
             mutate(miseq_run=str_extract(sample_id,'\\M\\d+'))
         }) %>% bind_rows()
fwrite(nReads_nCLones_allFiles,'/home/lera/poised_trb/analysis/nReads_nClones_POISED_data.txt')

nReads_nCLones_allFiles %>% group_by(miseq_run) %>% 
  summarise(n_samples=n_distinct(sample_id),
            mean_nReads=mean(nReads),
            mean_nClones=mean(nClones)) %>% 
  fwrite('/home/lera/poised_trb/analysis/nReads_MiSeq_stat_POISED_data.txt')

#pict reads/clones count----
read_tsv(file = dir("/home/lerucha38/POISED_TRB/mixcr",pattern=".txt",
                    full.names = T),
         id="sampleId") %>%
  mutate(sampleId=str_remove_all(sampleId,"\\/.*\\/|.clonotypes.TRB.txt") ) %>% 
  group_by(sampleId) %>% 
  summarise(nClones=n(),nReads=sum(cloneCount)) %>% 
  pivot_longer(c(nClones,nReads)) %>% 
  ggplot(aes(x=sampleId,y=value))+
  geom_bar(stat="identity",fill="navyblue")+
  facet_wrap(~name,nrow=2)+
  theme(axis.text.x = element_text(angle = 90))

#tissue/pbmc overlap analysis----

##overlap by patient in different time points----

rich_patients<-final_metadata %>%
  group_by(participant_label,week104_FC_result,Group) %>% 
  count() %>% 
  filter(n>9) %>% 
  pull(participant_label)
  
rich_patients_files<-final_metadata %>% 
  filter( participant_label %in% rich_patients)  %>% 
  filter(participant_label!='BFI_0003804',
         week104_FC_result!='NO_CHALLENGE')

table_for_analysis %>% 
  ggplot(aes(x=participant_label,y=time_point,color=tissue))+
  scale_color_viridis_d(alpha=0.65,begin = 0,
                       end = 0.75)+
  geom_quasirandom(aes(size=N))+
  facet_wrap(~week104_FC_result,scales='free_x')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

#for downsample
samples_mixcr_for_analysis<- rich_patients_files %>%
 mutate(path='/home/lerucha38/POISED_TRB/fasta_poised/mixcr/',
        whole_path=paste0(path,sample)) %>% 
                 pull(whole_path)

samples_all <- read_tsv(samples_mixcr_for_analysis,id = "filename") #это целиком все файлы для оверлэпа в одной таблице

samples<-samples_all %>% 
  mutate(filename=str_remove_all(filename,"\\/.*\\/")) #преобразую файлнейм,чтобы потом из
                                                         #него брать образцы
all_closets_rich_patients <- samples %>% 
  select(filename,cloneId,cloneCount,cloneFraction,bestVGene, #работаю с таблицей для вида
         bestJGene,nSeqCDR3,aaSeqCDR3) %>% 
  dplyr::rename(count=cloneCount,
                freq=cloneFraction,
                V=bestVGene,
                J=bestJGene,
                CDR3nt=nSeqCDR3,
                CDR3aa=aaSeqCDR3) %>% 
  mutate_at(c('V','J'),
            ~str_replace_all(string=., pattern='\\*.*',replacement = '')) %>% 
  filter(filename!='M157_S026_BFI_0003982.clonotypes.TRB.txt',
         filename!='M288_S077_BFI_0003983.clonotypes.TRB.txt',
         filename!='M288_S079_BFI_0003983.clonotypes.TRB.txt',
         filename!='M168_S012_BFI_0003984.clonotypes.TRB.txt')

list_samples<-samples %>% 
  group_by(filename) %>% 
  summarise(n_clones=n(),
            n_reads=sum(cloneCount))

table_for_analysis<-list_samples %>% 
  dplyr::rename(sample=filename) %>% 
  inner_join(final_metadata %>% 
               select(tissue,time_point,Group,week104_FC_result,
                      sample,participant_label,week),by='sample') %>%
  filter(tissue!='Duodenum_Stomach')

table_for_analysis %>% #картинка где пациенты по тканям ,чтобы показать глубину секвенирования
  mutate(sample=str_extract(sample,pattern = '.*_B')) %>% 
  mutate(sample=str_remove_all(sample,pattern = '_B')) %>% 
  mutate(filename=paste0(sample,'_',week)) %>% 
  mutate(filename=as.factor(filename)) %>% 
  mutate(filename = fct_reorder(filename,week)) %>% 
  ggplot(aes(x=filename,y=n_reads,fill=participant_label))+
  scale_fill_viridis_d(alpha=0.65,begin = 0,
                       end = 0.75) +
  geom_bar(stat = 'identity')+
  theme_bw()+
  facet_wrap(~tissue,scales='free_x')+
  theme(axis.text.x = element_text(angle = 90))
  #theme(axis.text.x = element_blank())

#code for overlap
clonsents_additionalData<- # полная табличка с клонсетами и метадатой для того,чтобы брать клоны
  all_closets_rich_patients %>% 
  dplyr::rename(sample=filename) %>% 
  inner_join(table_for_analysis %>% select(tissue,sample,
                                           time_point,participant_label,Group,week104_FC_result,
                                           Group3,week), by='sample') %>% 
  filter(week!='W117B',
         week!='W104B')


lapply(participant_labels, function(part_label){
  print(part_label)
  time_points<- clonsents_additionalData %>% 
    filter(participant_label==part_label) %>% 
              pull(time_point) %>% unique()

   lapply( time_points,function(times){
    print(times)
    all_tissues<-clonsents_additionalData %>% filter(participant_label==part_label,
                                                     time_point==times) %>% 
                                   pull(tissue) %>% unique()
    print(all_tissues)
    if(length(all_tissues)>1){
      mycomparisons<-combn(all_tissues,2,c,simplify=F)
      right_comp<-mycomparisons[grep('PBMC',mycomparisons)] 
      if (length(right_comp)!=0){
        lapply(right_comp,function(right_c){
         table_clonsents_1<-clonsents_additionalData %>% filter(participant_label==part_label,
                                                             time_point==times,
                                                             tissue==right_c[1])
          table_clonsents_2<- clonsents_additionalData %>% filter(participant_label==part_label,
                                                              time_point==times,
                                                              tissue==right_c[2])
         common_clones<-
              inner_join(table_clonsents_1,table_clonsents_2, by=c('CDR3aa','V','J'))
         print(nrow(common_clones))
          common_clones_new<-
             common_clones %>% 
              mutate(overlap_number=(nrow(common_clones))/
                       (as.numeric((nrow(table_clonsents_1)))*(as.numeric((nrow(table_clonsents_2)))))) %>% 
            mutate(n_row_table1=nrow(table_clonsents_1),
                   n_row_table2=nrow(table_clonsents_2))
        }) %>% bind_rows()}
     else {
       table_clonsent_1<-clonsents_additionalData %>% filter(participant_label=='BFI_0003980')
       table_clonsent_2<-clonsents_additionalData %>% filter(participant_label=='BFI_000398')
      df<- inner_join(table_clonsent_1,table_clonsent_2,by=c('CDR3aa','V','J')) %>% 
        mutate(overlap_number=as.double(),
               n_row_table1=as.integer(),
               n_row_table2= as.integer())
      return(df)
      }}
    else {
      table_clonsent_1<-clonsents_additionalData %>% filter(participant_label=='BFI_0003980')
      table_clonsent_2<-clonsents_additionalData %>% filter(participant_label=='BFI_000398')
      DF<- inner_join(table_clonsent_1,table_clonsent_2,by=c('CDR3aa','V','J')) %>% 
        mutate(overlap_number=as.double(),
               n_row_table1=as.integer(),
               n_row_table2= as.integer())
      return(DF)
      }
      })  %>% bind_rows() 
  }) %>%  bind_rows() %>% 
    fwrite(.,paste0('/home/lerucha38/POISED_TRB/analysis/','overlaps_allTissues.txt'))
  
all_patients_overlaps<-fread('/home/lerucha38/POISED_TRB/analysis/overlaps_allTissues.txt')

deduple_overlaps<-all_patients_overlaps %>% 
  dplyr::rename(time_point=time_point.x,
                participant_label=participant_label.x,
                week104_FC_result=week104_FC_result.x,
                Group3=Group3.x,
                Group=Group.x) %>% 
  mutate(pair=paste0(tissue.x,'_',tissue.y),
         pair=case_when(pair=='dist_esophagus_PBMC'~'PBMC_dist_esophagus',
                        pair=='mid_esophagus_PBMC'~'PBMC_mid_esophagus',
                        pair=='prox_esophagus_PBMC'~'PBMC_prox_esophagus',
                        pair=='duodenum_PBMC'~'PBMC_duodenum',
                        pair=='stomach_PBMC'~'PBMC_stomach',
                        TRUE~pair))

plasma_pal <- c("red", viridis::plasma(n = 6))

deduple_overlaps %>%  
ggplot(aes(x=time_point,y=overlap_number))+
  geom_boxplot(aes(fill=Group),outlier.shape = NA,)+
  geom_point(aes(color = week104_FC_result),position=position_jitterdodge(jitter.width = 0),
             size=0.5)+
  scale_color_manual(values = plasma_pal)+
  facet_wrap(~pair,nrow=2)+
  scale_fill_viridis_d(alpha=0.45,begin = 0,
                        end = 0.75)+
  theme_bw()
  #ggtitle('patient BFI_0003980 success group A penaut 0')
  


#HLA----
hla_table<-fread('/data/vskatova/vskatova/poised_trb/analysis/HLA_typing_poised.csv',header=TRUE) %>% 
  select(-V20,-V21,-V22,-V23,-V24) %>% 
  dplyr::rename(Participant_PPID='Study ID') #прочиатем файлик

cnames<- str_replace_all(colnames(hla_table),pattern='\\/',replacement='_')#поменяем навание колонок, а то неозможно
colnames(hla_table)<-cnames

hla_table %<>% 
  filter(Participant_PPID!='') %>%  
  filter(Participant_PPID!='P002') %>% 
  inner_join(final_metadata %>% 
               select(age,Group,multi_allergy,Participant_PPID) %>% 
               distinct(.keep_all = TRUE) ,by='Participant_PPID')



##pair-wise comparison in allergic,non-allergic patients----
#нужно показать,что внутри аллергиков клоны более похожи,чем между аллергиком и неаллергиком
#32 пациента в каждой группе, аллергики на неделе 0
###week 0 ----
healthy_pat_list<-healthy_clns_dataset_from_Maxim %>% 
  pull(participant_label) %>% unique() %>% head(32)

healthy_donor_pairs<-combn(healthy_pat_list,m=2,FUN = c,simplify = F)
allergic_donor_pairs<-combn(clean_allergic_counted_w0_aa %>% pull(patients_label) %>% unique(),
                            m=2,FUN = c,simplify = F)
healthy_allergic_donor_pairs<-combn(c(healthy_pat_list,clean_allergic_counted_w0_aa %>% pull(patients_label) %>% unique()),
                                    m=2,FUN = c,simplify = F)

pair_counted<-lapply(c(healthy_donor_pairs,allergic_donor_pairs,
         healthy_allergic_donor_pairs),function(pair){
print(pair)
    common_dataset<-allergic_clns %>% bind_rows(healthy_clns_dataset_from_Maxim %>% 
                dplyr::rename(V=bestVGene,
                              J=bestJGene,
                              CDR3nt=nSeqCDR3,
                              CDR3aa=aaSeqCDR3,
                              frequency=cloneFraction,
                              count=cloneCount))
    table_1<-common_dataset %>% 
      filter(participant_label==pair[1])
    print(paste0('1st dt is',nrow(table_1)))
    
    table_2<-common_dataset %>% 
      filter(participant_label==pair[2])
    print(paste0('2nd dt is',nrow(table_2)))
    
    inner_join_table<-inner_join(table_1,table_2,by=c('V','J','CDR3aa')) %>% 
      mutate(pat_id1=pair[1],
             pat_id2=pair[2]) %>% 
      count(participant_label.x,participant_label.y)
  }) %>% bind_rows()

pair_counted_precessed<-pair_counted %>% 
  mutate(present_in_healthy= participant_label.x %in% healthy_clns_dataset_from_Maxim$participant_label,
         present_in_healthy2= participant_label.y %in% healthy_clns_dataset_from_Maxim$participant_label,
         status_pat1=ifelse(present_in_healthy==TRUE,'healthy','pa_allergic'),
         status_pat2= ifelse(present_in_healthy2==TRUE,'healthy','pa_allergic'),
         pair_group=case_when(status_pat1=='healthy' & status_pat2=='healthy'~'h-h',
                              status_pat1=='pa_allergic' & status_pat2=='pa_allergic'~'pa-pa',
                              status_pat1=='pa_allergic' & status_pat2=='healthy'~'pa-h',
                              status_pat1=='healthy' & status_pat2=='pa_allergic'~'pa-h'))
  
pair_counted_precessed<-pair_counted_precessed[!duplicated(pair_counted_precessed )] 

statTest_w0_parewised<-dunn_test(n~pair_group,
                                              data=pair_counted_precessed, 
                                              p.adjust.method="bonferroni") %>% 
  add_y_position(formula=n~pair_group,data=pair_counted_precessed)
#graph
pair_counted_precessed %>% 
  ggplot(aes(y=n,x=pair_group))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=0.9)+
  geom_text(aes(x=pair_group,label=n,y=1000),
            data=pair_counted_precessed %>% count(pair_group) ,inherit.aes=FALSE)+
  stat_compare_means(label.y = 2000)+
  stat_pvalue_manual(data=statTest_w0_parewised, 
                   label = "p.adj.signif", tip.length = 0.01,hide.ns=TRUE)
  #ylab('number of shared clonotypes among the top 30000 clns, 
       #week 104')+
  #xlab('Group')

#пытаюсь разобраться,почему в группе аллергиков нет гомогенности
#самая очевидная гипотеза- hla
###Hla processing----
hla_table_processed<-hla_table %>%
  mutate(across(where(is.character),
             ~str_replace_all(string=.,pattern ='Insufficient data',
                                                   replacement = ''))) %>% 
  mutate(across(where(is.character),
                ~str_replace_all(string=.,pattern ='-',
                                 replacement = ''))) %>% 
  mutate(across(where(is.character),
                ~str_replace_all(string=.,pattern ='Missing DRB4*',
                                 replacement = ''))) %>% 
  mutate(across(where(is.character),
                ~str_replace_all(string=.,pattern ='Absent',
                                 replacement = '')))
hla_table_processed %<>% 
mutate(across(where(is.character),
                ~str_extract(string=.,pattern ='\\w+.\\d+'))) %>% 
  mutate(hla_vec = mapply(c,A_1,A_2,B_1,B_2,C_1,C_2,DPA1_1,DPA1_2,DPB1_1,DPB1_2,
                         DQA1_1,DQA1_2,DQB1_1,DQB1_2,
                         DRB1_1,DRB1_2,DRBo_1,DRBo_2, SIMPLIFY = FALSE)) %>% #делаю вектор из аллелей
  select(-hlavec,-age,-Group,-multi_allergy)

donor_pairs<-combn(hla_table_processed$Participant_PPID,m = 2,FUN = c,simplify = F)#лист пар доноров

shared_alleles_among_allergic_donors<-lapply(donor_pairs,function(pair){
  hla_vec_1<-hla_table_processed %>% filter(Participant_PPID==pair[1]) %>% 
    pull() %>% unlist(use.names = FALSE)
 
  hla_vec_2<-hla_table_processed %>% filter(Participant_PPID==pair[2]) %>% 
    pull(hla_vec) %>% unlist(use.names = FALSE)

intersaction_HLA<- intersect(hla_vec_1,hla_vec_2)
  print(intersaction_HLA)

intersaction_HLA %<>% 
    as.data.frame() %>% 
  dplyr::rename(int_pairs='.') %>% 
  mutate(donor_1=pair[1],
       donor_2=pair[2],
      n_int=length(int_pairs))
}) %>% bind_rows()

summ_shared_alleles_among_allergic_donors<-
  shared_alleles_among_allergic_donors %>% 
  group_by(donor_1,donor_2,n_int) %>% 
  summarise(shared_alleles=paste(int_pairs,collapse = ',')) #здесь нет нулей

donor_pairs_df<-donor_pairs %>% #составляю таблицу, в которой все пары пациентов и колич пересеченных HLA
  as.data.frame() %>% 
  pivot_longer(cols=everything(),names_to = 'pairs') %>%
  select(-value) %>% 
  unique() %>% 
  mutate(pairs=str_remove(pairs,'c..'),
    donor_1=str_extract(pairs,pattern='\\w+\\d+'),
         donor_2=str_extract(pairs,pattern='....\\w+\\d+'),
    donor_2=str_remove(donor_2,'....')) %>%
  select(-pairs) %>% 
  left_join(summ_shared_alleles_among_allergic_donors,
             by=c('donor_1','donor_2')) %>% 
  mutate(n_int=ifelse(is.na(n_int),0,n_int),
         NA_cor=str_extract(shared_alleles,'NA'),
         n_int=case_when(NA_cor=='NA'~(n_int-1),
                         TRUE~n_int ),
         shared_alleles=str_remove(shared_alleles,'NA')) %>% 
  select(-NA_cor)

allergic_clns %<>% #добавить короткий ID пациента в таблицу с аллергиескими клонсетами
  inner_join(final_metadata %>% select(participant_label,Participant_PPID) %>% 
               distinct(Participant_PPID, .keep_all = TRUE),
             by='participant_label')

#надо проверить, сколько клонов и ридов в файлах у этих пациентов, у которых есть HLA
donors_with_hla<-donor_pairs_df %>% 
  dplyr::rename(Participant_PPID=donor_1) %>% 
  inner_join(final_metadata %>% 
               mutate(week=case_when(Participant_PPID=='P089'& week =='BL01'~'W000',#немного шаманства
                                     Participant_PPID=='P091'& week =='BL04'~'W000',
                                     Participant_PPID=='P100'& week =='BL01'~'W000',
                                     TRUE~week)) %>% 
               filter(tissue=='PBMC',week=='W000') %>% 
               select(tissue,week,Group,sample,Participant_PPID), by='Participant_PPID')

read_tsv(donors_with_hla %>% #посмотреть на глубину этих образцов, у которых есть HLA
           mutate(path=paste0('/home/lerucha38/POISED_TRB/fasta_poised/vdjtools_format/',sample)) %>% 
           pull(path) %>% unique(),
         id="sampleId") %>%
  mutate(sampleId=str_remove_all(sampleId,"\\/.*\\/|.txt") ) %>% 
  group_by(sampleId) %>% 
  summarise(nClones=n(),nReads=sum(count)) %>% 
  pivot_longer(c(nClones,nReads)) %>% 
  ggplot(aes(x=sampleId,y=value))+
  geom_bar(stat="identity",fill="navyblue")+
  facet_wrap(~name,nrow=2)+
  theme(axis.text.x = element_text(angle = 90))

allergic_clones_w0_hlaPresent<-lapply(donors_with_hla %>% 
                                          mutate(path=paste0('/home/lerucha38/POISED_TRB/fasta_poised/vdjtools_format/',sample)) %>% 
                                          pull(path) %>% unique(),function(file){
                                                print(file)
                                                table1<-fread(file) 
                                                print(nrow(table1))
                                                if (nrow(table1)>10000) {
                                                  table1 %<>%  mutate(sample_id=file) %>% 
                                                    group_by(CDR3nt,V,J,D,CDR3aa,sample_id) %>% 
                                                    summarise(frequency=sum(frequency),# убираю дупликаты
                                                              count=sum(count)) %>% 
                                                    ungroup() %>% 
                                                    mutate(sample_id=str_remove_all(sample_id,"\\/.*\\/"),
                                                           participant_label=str_extract(sample_id,pattern='BFI_\\d+')) %>% 
                                                    head(10000)#другой топ для получения большего количества обрз-ов
                                                }
                                                else {
                                                  return(table1 %<>% head(0))
                                                }   
                                              }) %>% bind_rows()

allergic_clones_w0_hlaPresent %<>% inner_join(final_metadata %>% #добавляю в клонсеты короткий ID
                                               select(Participant_PPID,participant_label) %>% 
                                               distinct(),
                                             by='participant_label')

donorHLA_pairs<-combn(allergic_clones_w0_hlaPresent %>% 
                        pull(Participant_PPID) %>% unique(),m = 2,FUN = c,simplify = F)
####comparison----
pair_counted_hla<-lapply(donorHLA_pairs,function(pair){
                         print(pair)
          
                         table_1<-allergic_clones_w0_hlaPresent %>% 
                           filter(Participant_PPID==pair[1])
                         print(paste0('1st dt is',nrow(table_1)))
                         
                         table_2<-allergic_clones_w0_hlaPresent %>% 
                           filter(Participant_PPID==pair[2])
                         print(paste0('2nd dt is',nrow(table_2)))
                         
                         inner_join_table<-inner_join(table_1,table_2,by=c('V','J','CDR3aa')) %>% 
                           mutate(pat_id1=pair[1],
                                  pat_id2=pair[2]) %>% 
                           count(Participant_PPID.x,Participant_PPID.y)
                       }) %>% bind_rows()

hla_clones_table<-pair_counted_hla %>% 
  dplyr::rename(donor_1=Participant_PPID.x,
                donor_2=Participant_PPID.y,
                number_of_shared_cl=n) %>% 
  inner_join(donor_pairs_df,by=c('donor_1','donor_2')) %>%
  dplyr::rename(n_shared_alleles=n_int)

#graph
hla_clones_table %>% 
  ggplot(aes(x=n_shared_alleles %>% as.factor(),
             y=number_of_shared_cl,group=n_shared_alleles))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=0.9)+
  geom_text(aes(x=n_shared_alleles %>% as.factor() ,label=n,y=400),
            data=hla_clones_table %>% count(n_shared_alleles) ,inherit.aes=FALSE)+
  stat_compare_means(label.y = 420)

###week 52----
#первый раз этот анализ был сделан на 15500 клонов,но я сейчас сделаю его на 3000 клонов, в том
#числе и на 3000 здоровых доноров
healthy_clns_dataset_from_Maxim_3000<-lapply(dir("/home/lerucha38/POISED_TRB/fasta_healthy/new_healthy/mixcr",pattern=".txt",
                                            full.names = T),function(file){
                                              print(file)
                                              table1<-fread(file) 
                                              print(nrow(table1))
                                              if (nrow(table1)>15500) {
                                                table1 %<>%  mutate(sample_id=file) %>% 
                                                  group_by(nSeqCDR3,bestVGene,bestJGene,bestDGene,aaSeqCDR3,sample_id) %>% 
                                                  summarise(cloneFraction=sum(cloneFraction),# убираю дупликаты
                                                            cloneCount=sum(cloneCount)) %>% 
                                                  ungroup() %>% 
                                                  mutate(sample_id=str_remove_all(sample_id,"\\/.*\\/"),
                                                         participant_label=str_extract(sample_id,pattern='BFI_\\d+')) %>% 
                                                  head(3000)
                                              }
                                              else {
                                                return(table1 %<>% select(-cloneId) %>% head(0))
                                              }   
                                            }) %>% bind_rows()

PA_clns_w52_top3000 %<>% inner_join(final_metadata %>% 
                                   select(Participant_PPID,participant_label) %>% unique(),
                                 by='participant_label')

PA_clns_w52_top3000 %<>% inner_join(final_metadata %>% 
                                    select(Group2,participant_label) %>% unique(),
                                  by='participant_label')

healthy_donor_pairs<-combn(healthy_clns_dataset_from_Maxim_3000 %>% 
                             pull(participant_label) %>% 
                             unique(),m=2,FUN = c,simplify = F)#пары здоровых

allergic_donor_pairs_placebo<-combn(PA_clns_w52_top3000 %>% #пары аллегиков на плацебо
                                      filter(Group2=='C') %>% 
                                    pull(participant_label) %>% unique(),
                            m=2,FUN = c,simplify = F)

allergic_donor_pairs_IT<-combn(PA_clns_w52_top3000 %>% #пары аллергиков на ИТ
                                      filter(Group2=='AB') %>% 
                                      pull(participant_label) %>% unique(),
                                    m=2,FUN = c,simplify = F)

allg_pl_IT_pairs<-combn(PA_clns_w52_top3000 %>% #пары аллегиков на плацебо+ал на ИТ
                                      pull(participant_label) %>% unique(),
                                    m=2,FUN = c,simplify = F)

healthy_allergic_placebo_pairs<-combn(c(healthy_clns_dataset_from_Maxim_3000 %>% 
                                          pull(participant_label) %>% 
                                          unique(),PA_clns_w52_top3000 %>% #пары здоровых+пл
                                                            filter(Group2=='C') %>% 
                                          pull(participant_label) %>% unique()),
                                    m=2,FUN = c,simplify = F)

healthy_allergic_IT_pairs<-combn(c(healthy_clns_dataset_from_Maxim_3000 %>% 
                                     pull(participant_label) %>% 
                                     unique(),PA_clns_w52_top3000 %>% #пары зд+ИТ
                                          filter(Group2=='AB') %>% 
                                          pull(participant_label) %>% unique()),
                                      m=2,FUN = c,simplify = F)

pair_counted_w52_top300<-lapply(c(healthy_donor_pairs,allergic_donor_pairs_placebo,
                           allergic_donor_pairs_IT,allg_pl_IT_pairs,
                           healthy_allergic_placebo_pairs,
                           healthy_allergic_IT_pairs),function(pair){
                         print(pair)
      common_dataset<-PA_clns_w52_top3000 %>% bind_rows(healthy_clns_dataset_from_Maxim_3000 %>% 
                                                                       dplyr::rename(V=bestVGene,
                                                                                     J=bestJGene,
                                                                                     CDR3nt=nSeqCDR3,
                                                                                     CDR3aa=aaSeqCDR3,
                                                                                     frequency=cloneFraction,
                                                                                     count=cloneCount))
                         table_1<-common_dataset %>% 
                           filter(participant_label==pair[1])
                         print(paste0('1st dt is ',nrow(table_1)))
                         
                         table_2<-common_dataset %>% 
                           filter(participant_label==pair[2])
                         print(paste0('2nd dt is ',nrow(table_2)))
                         
                         inner_join_table<-inner_join(table_1,table_2,by=c('V','J','CDR3aa')) %>% 
                           mutate(pat_id1=pair[1],
                                  pat_id2=pair[2]) %>% 
                           count(participant_label.x,participant_label.y)
                       }) %>% bind_rows()

pair_counted_w52_top300_pr <-  #есть нули
  pair_counted_w52_top300[!duplicated(pair_counted_w52_top300 )] 

pair_counted_w52_top300_processed<- c(healthy_donor_pairs,allergic_donor_pairs_placebo,
                               allergic_donor_pairs_IT,allg_pl_IT_pairs,
                               healthy_allergic_placebo_pairs,
                               healthy_allergic_IT_pairs) %>% as.data.frame() %>% 
  pivot_longer(cols=everything(),names_to = 'pairs',values_to = 'pat_id_1') %>%
  select(-pat_id_1) %>% 
  unique() %>% 
  mutate(pairs=str_remove(pairs,'c..'),
         donor_1=str_extract(pairs,pattern='\\w+\\d+'),
         donor_2=str_extract(pairs,pattern='\\.\\w+\\d+'), 
         donor_2=str_remove(donor_2,'\\.')) %>% 
  select(-pairs) %>%  
  unique() %>% 
  left_join(pair_counted_w52_top300_pr %>% 
              dplyr::rename(donor_1=participant_label.x,
                            donor_2= participant_label.y),
            by=c('donor_1','donor_2')) %>% 
  unique() %>% 
  mutate(n=ifelse(is.na(n),'0',n))

pair_counted_w52_top300_processed_new<-pair_counted_w52_top300_processed %>% 
  mutate(present_in_healthy= donor_1 %in% (healthy_clns_dataset_from_Maxim_3000 %>% 
           pull(participant_label) %>% 
           unique()),
         present_in_healthy2= donor_2 %in% (healthy_clns_dataset_from_Maxim_3000 %>% 
           pull(participant_label) %>% 
           unique()),
         status_pat1=ifelse(present_in_healthy==TRUE,'healthy','pa_allergic'),
         status_pat2= ifelse(present_in_healthy2==TRUE,'healthy','pa_allergic'),
         present_in_IT=donor_1 %in% (PA_clns_w52_top3000 %>% filter(Group2=='AB') %>% 
                                                   pull(participant_label)),
         present_in_IT_2=donor_2 %in% (PA_clns_w52_top3000 %>% filter(Group2=='AB') %>% 
                                                   pull(participant_label)),
         present_in_placebo=donor_1 %in% (PA_clns_w52_top3000 %>% filter(Group2=='C') %>% 
                                                   pull(participant_label)),
         present_in_placebo_2=donor_2 %in% (PA_clns_w52_top3000 %>% filter(Group2=='C') %>% 
                                                        pull(participant_label)),
         therapy_status_1=case_when(status_pat1=='pa_allergic'& present_in_IT==FALSE & 
                                      present_in_placebo==TRUE ~'placebo',
                                    status_pat1=='pa_allergic'& present_in_IT==TRUE & 
                                      present_in_placebo==FALSE ~'IT',
                                    status_pat1=='healthy'~'healthy' ),
         therapy_status_2=case_when(status_pat2=='pa_allergic' & present_in_IT_2==FALSE & 
                                      present_in_placebo_2==TRUE ~'placebo',
                                    status_pat2=='pa_allergic'& present_in_IT_2==TRUE & 
                                      present_in_placebo_2==FALSE ~'IT',
                                    status_pat2=='healthy'~'healthy'),
         pair_group=case_when(therapy_status_1=='healthy' & therapy_status_2=='healthy'~'h-h',
                              therapy_status_1=='healthy' & therapy_status_2=='placebo'~'pl-healthy',
                              therapy_status_1=='placebo' & therapy_status_2=='healthy'~'pl-healthy',
                              therapy_status_1=='healthy' & therapy_status_2=='IT'~'IT-healthy',
                              therapy_status_1=='IT' & therapy_status_2=='healthy'~'IT-healthy',
                              therapy_status_1=='placebo' & therapy_status_2=='placebo'~'pl-pl',
                              therapy_status_1=='IT' & therapy_status_2=='IT'~'IT-IT',
                              therapy_status_1=='placebo' & therapy_status_2=='IT'~'pl-IT',
                              therapy_status_1=='IT' & therapy_status_2=='placebo'~'pl-IT' )) %>% 
  mutate(n=as.integer(n))

statTest_w52_pairwised<-dunn_test(n~pair_group,
                                 data=pair_counted_w52_top300_processed_new, 
                                 p.adjust.method="bonferroni") %>% 
  add_y_position(formula=n~pair_group,data=pair_counted_w52_top300_processed_new)
#graph
pair_counted_w52_top300_processed_new %>% 
  mutate(pair_group=factor(x=pair_group,
                           levels=c('h-h','pl-pl','IT-IT','pl-IT','pl-healthy','IT-healthy'))) %>% 
  ggplot(aes(y=n,x=pair_group))+
  geom_point(size=0.8,position=position_jitter(width = 0.3,height = 0.5))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(aes(x=pair_group,label=n,y=-10),
            data=pair_counted_w52_top300_processed_new %>% count(pair_group) ,inherit.aes=FALSE)+
  #stat_compare_means()+
  stat_pvalue_manual(data=statTest_w52_pairwised, 
                     label = "p.adj.signif", tip.length = 0.01,hide.ns=FALSE)+
  #ylim(-10,50)+ #лимит,чтобы распределение посмотреть
ylab('numb of shared clonotypes, week 52')+
xlab('Group')
  
###week 104----
PA_clns_w104_top5000<-lapply(dir("/home/lerucha38/POISED_TRB/fasta_poised/vdjtools_format/downsample_0_52_104",pattern="M.*",
                                 full.names = T),function(file){
                                   print(file)
                                   table1<-fread(file) %>% 
                                     mutate(sample_id=file) %>% 
                                     group_by(CDR3nt,V,J,CDR3aa,sample_id) %>% 
                                     summarise(frequency=sum(frequency),# убираю дупликаты клонов
                                               count=sum(count)) %>% 
                                     ungroup() %>% 
                                     mutate(sample_id=str_remove_all(sample_id,"\\/.*\\/"),
                                            participant_label=str_extract(sample_id,pattern='BFI_\\d+')) %>% 
                                     inner_join(final_metadata %>% select(week,sample,Group) %>% 
                                                  dplyr::rename(sample_id=sample),by='sample_id') %>% 
                                     filter(week=='W104'|week=='W104A' ) %>% #беру только week104
                                     head(5000) #беру только 5000 клонов
                                 }) %>% rbindlist()

healthy_clns_dataset_from_Maxim_5000<-lapply(dir("/home/lerucha38/POISED_TRB/fasta_healthy/new_healthy/mixcr",pattern=".txt",
                                                 full.names = T),function(file){
                                                   print(file)
                                                   table1<-fread(file) 
                                                   print(nrow(table1))
                                                   if (nrow(table1)>15500) {
                                                     table1 %<>%  mutate(sample_id=file) %>% 
                                                       group_by(nSeqCDR3,bestVGene,bestJGene,bestDGene,aaSeqCDR3,sample_id) %>% 
                                                       summarise(cloneFraction=sum(cloneFraction),# убираю дупликаты
                                                                 cloneCount=sum(cloneCount)) %>% 
                                                       ungroup() %>% 
                                                       mutate(sample_id=str_remove_all(sample_id,"\\/.*\\/"),
                                                              participant_label=str_extract(sample_id,pattern='BFI_\\d+')) %>% 
                                                       head(5000)
                                                   }
                                                   else {
                                                     return(table1 %<>% select(-cloneId) %>% head(0))
                                                   }   
                                                 }) %>% bind_rows()

PA_clns_w104_top5000 %<>% inner_join(final_metadata %>% 
                                      select(Participant_PPID,participant_label) %>% unique(),
                                    by='participant_label')
PA_clns_w104_top5000 %<>% inner_join(final_metadata %>% 
                                      select(Group2,participant_label) %>% unique(),
                                    by='participant_label')

pair_counted_withHealthy_w104_top5000<-lapply(c(healthy_donor_pairs,allergic_donor_pairs_placebo,
                                  allergic_donor_pairs_IT,allg_pl_IT_pairs,
                                  healthy_allergic_placebo_pairs,
                                  healthy_allergic_IT_pairs),function(pair){
                                    print(pair)
      common_dataset<-PA_clns_w104_top5000 %>% bind_rows(healthy_clns_dataset_from_Maxim_5000 %>% 
                                                                                        dplyr::rename(V=bestVGene,
                                                                                                      J=bestJGene,
                                                                                                      CDR3nt=nSeqCDR3,
                                                                                                      CDR3aa=aaSeqCDR3,
                                                                                                      frequency=cloneFraction,
                                                                                                      count=cloneCount))
                                    table_1<-common_dataset %>% 
                                      filter(participant_label==pair[1])
                                    print(paste0('1st dt is ',nrow(table_1)))
                                    
                                    table_2<-common_dataset %>% 
                                      filter(participant_label==pair[2])
                                    print(paste0('2nd dt is ',nrow(table_2)))
                                    
                                    inner_join_table<-inner_join(table_1,table_2,by=c('V','J','CDR3aa')) %>% 
                                      mutate(pat_id1=pair[1],
                                             pat_id2=pair[2]) %>% 
                                      count(participant_label.x,participant_label.y)
                                  }) %>% bind_rows()



#DEEP SEQ ---- 
#the 1st run w104 and w117----
##sample sheet for MiXCR----
sampleTable<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/sample_table.tsv') 
sampleTable %>%   mutate(SAMPLE1_rev_compl=sapply(sampleTable %>% pull(SAMPLE1), 
                                                  FUN = function(x) as.character(reverseComplement(DNAString(x))))) %>% 
  select(-SAMPLE1) %>% 
  dplyr::rename(SAMPLE1=SAMPLE1_rev_compl) %>% 
  select(Sample,TagPattern,SAMPLE1,SAMPLE2) %>% 
  dplyr::rename(SAMPLE1_old=SAMPLE1, SAMPLE2_old=SAMPLE2) %>% 
  dplyr::rename(SAMPLE1=SAMPLE2_old,SAMPLE2=SAMPLE1_old) %>% 
  select(Sample,TagPattern,SAMPLE1,SAMPLE2) %>% 
  select(Sample,TagPattern,SAMPLE1) %>% 
  write_tsv('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/sample_table_1.tsv')

## rarefaction plot----
dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr', pattern = '.clns') %>% 
  as.data.frame() %>% 
  dplyr::rename(sample_downsample='.') %>% 
  mutate(sample_names=str_detect(sample_downsample,'result.*')) %>% 
  filter(sample_names==TRUE) %>% 
  write_tsv('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/rare_faction_samples.txt',col_names = FALSE)  

report_downsample_deep_deq<-read_tsv(file=dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/downsample_run1/',
                                              pattern = '.tsv',full.names = TRUE), id = 'filename')

report_downsample_deep_deq %<>% 
  mutate(week=str_extract(sample,'W\\d+')) %>% 
  select(sample,sumWeightAfter,nElementsAfter,week) %>% 
  dplyr::rename(N_sequenced_reads=sumWeightAfter,
                N_clones=nElementsAfter)
report_downsample_deep_deq %>% 
  write_tsv('/data/vskatova/vskatova/data_milab_plots/data_rarefactions_curves.tsv')

report_downsample_deep_deq %>% 
  ggplot(aes(x = N_sequenced_reads, y = N_clones, color = sample)) +
  facet_grid(~week)+
  theme_bw() + 
  scale_color_discrete(guide = 'none') +  # turn legend on or off
  geom_line() +
  geom_vline(xintercept=500000, color= "red", linetype='dashed') + 
  labs(title="Rarefaction curves") + 
  xlab("N sequenced Reads") + 
  ylab('N clones')

ggsave(filename = '/data/vskatova/vskatova/data_milab_plots/Rarefaction_curves_run1.png',
       width = 10,height = 7,device = 'png',dpi = 400)

##metadata----

first_run_files_deepSeq<-dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/',
                             pattern = '.clns.tsv') %>% 
  as.data.frame() %>% 
  dplyr::rename(sample_downsample='.') %>% 
  mutate(sample_names=str_detect(sample_downsample,'.clns.tsv')) %>% 
  filter(sample_names==TRUE) %>% 
  mutate(Participant_PPID=str_extract(sample_downsample,'P\\d+')) %>% 
  mutate(week=str_extract(sample_downsample,'W\\d+.*_')) %>%
  mutate(week=str_remove(week,'_PBMC_')) %>% 
  mutate(tissue='PBMC')

deep_seq_run1_metadata<-final_metadata %>%  #final metadata for the first seq run
  inner_join(first_run_files_deepSeq,by=c('week','tissue','Participant_PPID')) %>% 
  unique() %>% 
  select(-sample,-sample_label,-time_point,-specimen_label,-Response.at.Wk104,
         -Response.at.Wk117,-Response.at.Wk0,-'Week 104.ChallengeResult',-specimen_description,
         -"Week 104.CTD","Week 117.CTD","Week 130.CTD",
         "Week 143.CTD","Week 156.CTD",-sample_names) %>% 
  dplyr::rename(sample_name=sample_downsample,
                week117_FC_result='Week 117.ChallengeResult',
                week130_FC_result='Week 130.ChallengeResult',
                week143_FC_result='Week 143.ChallengeResult',
                week156_FC_result ='Week 156.ChallengeResult') %>% 
  select(sample_name,Participant_PPID,participant_label,week,tissue,Group,Group2,Group3,age,
         multi_allergy,Sex,count_allergens,week104_FC_result,
         week117_FC_result,week143_FC_result,week156_FC_result,
         percentCD63.high.at.BL,percentCD63.high.Wk0,percentCD63.high.Wk52,
         percentCD63.high.Wk104,percentCD63.high.Wk117) %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/',
                          sample_name))

first_run_deep_seq_metadata<-deep_seq_run1_metadata %>% #to add downsampled data path
  mutate(downsampled_500k_path=str_replace(real_path,'/POIS','/downsample/result.POIS'),
         downsampled_500k_path=str_replace(downsampled_500k_path,'.clns.tsv','.trb.tsv')) %>% 
  mutate(week=case_when((Participant_PPID=='P104'& week=='W104A') ~'W104', #to change because there is no w104B sample
                        (Participant_PPID=='P114' & week=='W104A') ~'W104',
                        (Participant_PPID=='P103' & week=='W104A') ~'W104',
                        (Participant_PPID=='P106' & week=='W104A') ~'W104',
                        TRUE~week)) %>% 
  mutate(week=case_when((Participant_PPID=='P099' & week=='W117A')~'W117', #to change because there is no w104B sample
                        (Participant_PPID=='P114'& week=='W117A') ~'W117',
                        (Participant_PPID=='P103'& week=='W117B')~'W117',
                        (Participant_PPID=='P110'& sample_name=='POIS_P110_W117_PBMC_2.clns.tsv')~'W117A',
                        (Participant_PPID=='P110'& sample_name=='POIS_P110_W117_PBMC_4.clns.tsv')~'W117B',
                        TRUE~week)) %>% unique() %>% 
  mutate(downsample_file_present=downsampled_500k_path %in% 
           dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/downsample',full.names = TRUE))

#write down the metadata
write_tsv(first_run_deep_seq_metadata,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/metadata_run1_deepSeq.tsv')

#the 2nd run BL and w36----
##sample sheet for MiXCR----
sampleTable<-read_excel('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/M355plusNova_barcode_map.xlsx')
sampleTable %<>%   
  dplyr::rename(Sample="Specimen_Parent Specimen Label",
                SAMPLE1="V barcode...10" ,
                SAMPLE2="C barcode...11") %>% 
  select(Sample,SAMPLE1,SAMPLE2) %>% 
  mutate(TagPattern='')
sampleTable %<>% 
  mutate(SAMPLE1_rev_compl=sapply(sampleTable %>% pull(SAMPLE1), 
                                  FUN = function(x) as.character(reverseComplement(DNAString(x))))) %>% 
  select(-SAMPLE1) %>% 
  dplyr::rename(SAMPLE1=SAMPLE1_rev_compl)

sampleTable %>%  
  dplyr::rename(SAMPLE1_old=SAMPLE1,
                SAMPLE2_old=SAMPLE2) %>% 
  dplyr::rename(SAMPLE1=SAMPLE2_old,
                SAMPLE2=SAMPLE1_old) %>% 
  select(Sample,TagPattern,SAMPLE1,SAMPLE2) %>% 
  write_tsv('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/sample_table_run2.tsv')

fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/sample_table_run2.tsv')

##report metrics the second run----

report_run2<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/report_mixcr_run2.tsv') 
report_run2_new<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/report_mixcr_run2_new_parameters.tsv')

##metadata----
second_run_files_deepSeq<-dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/', pattern = '.trb.tsv') %>% 
  as.data.frame() %>% 
  dplyr::rename(sample_downsample='.') %>% 
  mutate(sample_names=str_detect(sample_downsample,'.trb.tsv')) %>% 
  filter(sample_names==TRUE) %>% 
  mutate(Participant_PPID=str_extract(sample_downsample,'P\\d+')) %>% 
  mutate(week=str_extract(sample_downsample,'W\\d+')) %>% 
  mutate(week=ifelse(is.na(week),str_extract(sample_downsample,'BL\\d+'),week)) %>% 
  mutate(tissue='PBMC') %>% 
  mutate(Participant_PPID=str_replace(Participant_PPID,'P','P0')) %>% 
  mutate(Participant_PPID=case_when(Participant_PPID=='P0005'~'P005',
                                    Participant_PPID=='P0007'~'P007',
                                    Participant_PPID=='P0008'~'P008',
                                    Participant_PPID=='P0009'~'P009',
                                    Participant_PPID=='P0012'~'P012',
                                    Participant_PPID=='P0011'~'P011',
                                    Participant_PPID=='P0019'~'P019',
                                    Participant_PPID=='P0017'~'P017',
                                    Participant_PPID=='P0020'~'P020',
                                    Participant_PPID=='P0018'~'P018',
                                    Participant_PPID=='P0021'~'P021',
                                    Participant_PPID=='P0022'~'P022',
                                    Participant_PPID=='P0023'~'P023',
                                    Participant_PPID=='P0024'~'P024',
                                    Participant_PPID=='P0025'~'P025',
                                    Participant_PPID=='P0026'~'P026',
                                    Participant_PPID=='P0029'~'P029',
                                    Participant_PPID=='P0030'~'P030',
                                    Participant_PPID=='P0031'~'P031',
                                    Participant_PPID=='P0034'~'P034',
                                    Participant_PPID=='P0035'~'P035',
                                    Participant_PPID=='P0036'~'P036',
                                    Participant_PPID=='P0037'~'P037',
                                    Participant_PPID=='P0038'~'P038',
                                    Participant_PPID=='P0039'~'P039',
                                    Participant_PPID=='P0040'~'P040',
                                    Participant_PPID=='P0043'~'P043',
                                    Participant_PPID=='P0044'~'P044',
                                    Participant_PPID=='P0045'~'P045',
                                    Participant_PPID=='P0047'~'P047',
                                    Participant_PPID=='P0048'~'P048',
                                    Participant_PPID=='P0049'~'P049',
                                    Participant_PPID=='P0050'~'P050',
                                    Participant_PPID=='P0051'~'P051',
                                    Participant_PPID=='P0052'~'P052',
                                    Participant_PPID=='P0053'~'P053',
                                    Participant_PPID=='P0055'~'P055',
                                    Participant_PPID=='P0056'~'P056',
                                    Participant_PPID=='P0057'~'P057',
                                    Participant_PPID=='P0059'~'P059',
                                    Participant_PPID=='P0060'~'P060',
                                    Participant_PPID=='P0061'~'P061',
                                    Participant_PPID=='P0062'~'P062',
                                    Participant_PPID=='P0063'~'P063',
                                    Participant_PPID=='P0064'~'P064',
                                    Participant_PPID=='P0065'~'P065',
                                    Participant_PPID=='P0066'~'P066',
                                    Participant_PPID=='P0067'~'P067',
                                    Participant_PPID=='P0068'~'P068',
                                    Participant_PPID=='P0069'~'P069',
                                    Participant_PPID=='P0070'~'P070',
                                    Participant_PPID=='P0071'~'P071',
                                    Participant_PPID=='P0072'~'P072',
                                    Participant_PPID=='P0073'~'P073',
                                    Participant_PPID=='P0074'~'P074',
                                    Participant_PPID=='P0075'~'P075',
                                    Participant_PPID=='P0076'~'P076',
                                    Participant_PPID=='P0077'~'P077',
                                    Participant_PPID=='P0078'~'P078',
                                    Participant_PPID=='P0079'~'P079',
                                    Participant_PPID=='P0080'~'P080',
                                    Participant_PPID=='P0081'~'P081',
                                    Participant_PPID=='P0082'~'P082',
                                    Participant_PPID=='P0083'~'P083',
                                    Participant_PPID=='P0084'~'P084',
                                    Participant_PPID=='P0085'~'P085',
                                    Participant_PPID=='P0086'~'P086',
                                    Participant_PPID=='P0087'~'P087',
                                    Participant_PPID=='P0088'~'P088',
                                    Participant_PPID=='P0089'~'P089',
                                    Participant_PPID=='P0090'~'P090',
                                    Participant_PPID=='P0091'~'P091',
                                    Participant_PPID=='P0092'~'P092',
                                    Participant_PPID=='P0093'~'P093',
                                    Participant_PPID=='P0094'~'P094',
                                    Participant_PPID=='P0095'~'P095',
                                    Participant_PPID=='P0096'~'P096',
                                    Participant_PPID=='P0097'~'P097',
                                    Participant_PPID=='P0098'~'P098',
                                    Participant_PPID=='P0099'~'P099',
                                    Participant_PPID=='P0100'~'P100',
                                    Participant_PPID=='P0101'~'P101',
                                    Participant_PPID=='P0102'~'P102',
                                    Participant_PPID=='P0103'~'P103',
                                    Participant_PPID=='P0104'~'P104',
                                    Participant_PPID=='P0105'~'P105',
                                    Participant_PPID=='P0106'~'P106',
                                    Participant_PPID=='P0107'~'P107',
                                    Participant_PPID=='P0108'~'P108',
                                    Participant_PPID=='P0109'~'P109',
                                    Participant_PPID=='P0110'~'P110',
                                    Participant_PPID=='P0111'~'P111',
                                    Participant_PPID=='P0112'~'P112',
                                    Participant_PPID=='P0113'~'P113',
                                    Participant_PPID=='P0114'~'P114',
                                    Participant_PPID=='P0115'~'P115',
                                    Participant_PPID=='P0116'~'P116',
                                    Participant_PPID=='P0117'~'P117',
                                    Participant_PPID=='P0118'~'P118',
                                    Participant_PPID=='P0119'~'P119',
                                    Participant_PPID=='P0120'~'P120',
                                    TRUE~Participant_PPID))

second_run_files_deepSeq_new<-final_metadata %>% select(Group,Group2,Group3,age,
                                                        'Week 104.ChallengeResult','Week 117.ChallengeResult','Week 130.ChallengeResult',
                                                        'Week 143.ChallengeResult','Week 156.ChallengeResult',Participant_PPID) %>%  #final metadata for the 2nd run 
  inner_join(second_run_files_deepSeq,by=c('Participant_PPID'),multiple='all') %>% 
  unique() %>% 
  dplyr::rename(sample_name=sample_downsample,
                week104_FC_result='Week 104.ChallengeResult',
                week117_FC_result='Week 117.ChallengeResult',
                week130_FC_result='Week 130.ChallengeResult',
                week143_FC_result='Week 143.ChallengeResult',
                week156_FC_result ='Week 156.ChallengeResult') %>% 
  select(-sample_names) %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/',
                          sample_name))

#to write down the metadata for the second run
write_tsv(second_run_files_deepSeq_new,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/metadata_run2_deepSeq.tsv')

bl_w36_clones<-read_tsv(file=second_run_files_deepSeq_new %>% pull(real_path),
                        id = 'real_path') %>% 
  inner_join(second_run_files_deepSeq_new,by='real_path') %>% 
  filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID) %>% 
  dplyr::rename(v_trb=allVHitsWithScore,
                j_trb=allJHitsWithScore,
                cdr3aa_trb=aaSeqCDR3,
                cdr3nt_trb=nSeqCDR3,
                count=readCount,
                freq=readFraction) %>% 
  select(-targetQualities,-allDHitsWithScore,-allVAlignments,-allDAlignments,-allJAlignments,-allCAlignments,
         -minQualCDR3,-refPoints,-targetSequences) %>% 
  mutate(v_trb=str_remove(v_trb,'\\*\\d+.*')) %>% 
  mutate(j_trb=str_remove(j_trb,'\\*\\d+.*')) %>% 
  group_by(real_path,cdr3aa_trb,v_trb,j_trb,sample_name,Participant_PPID,
           week,Group,Group2,Group3,age,week104_FC_result,week117_FC_result,
           week130_FC_result, week143_FC_result,week156_FC_result) %>% 
  summarise(count_bulk=sum(count),
            freq_bulk=sum(freq),
            n_cdr3nt_bulk=n_distinct(cdr3nt_trb)) %>% 
  ungroup() %>% 
  mutate(vj_trb=paste0(v_trb,'_',j_trb))

#the 3rd run w52 and w0----

##barcodes----
C_barcode_file<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/barcode_list.tsv',header = TRUE) %>% 
  dplyr::rename(C_primer_barcode=id,
                C_primer=seq)

C_barcode_file_2<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/2a_c_primers.tsv',header = TRUE) %>% 
  select(C_primer_barcode,C_primer) %>% 
  mutate(C_primer=str_replace(C_primer,'CGACCTC','')) %>% 
  bind_rows(C_barcode_file) %>% 
  mutate(C_primer_barcode=ifelse(str_detect(C_primer_barcode,"10"),C_primer_barcode, str_remove(C_primer_barcode,"0")))

V_barcodes<-sampleTable %>% select("V barcode...8","V barcode...10") %>% 
  dplyr::rename(V_primer_mix="V barcode...8",
                V_barcode="V barcode...10") %>% unique()

miseq_inf<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/Miseq Information-051222.csv') 
  dplyr::rename(V_primer_mix='V primer mix',
                C_primer_barcode='C primer barcode') %>% 
  mutate(C_primer_barcode=ifelse(str_detect(C_primer_barcode,"10"),C_primer_barcode, str_remove(C_primer_barcode,"0"))) %>% 
  left_join(barcode_file,by=c('C_primer_barcode')) %>% 
  inner_join(V_barcodes,by='V_primer_mix') %>% 
  mutate(sample_label=str_replace_all(sample_label,'-','_')) %>% 
  mutate(sample_label=str_replace(sample_label,'15','015')) %>% 
  mutate(sample_label=str_replace(sample_label,'0015','015'))

miseq_run3_corrected<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/miseq_inf_new.tsv') |> 
  select(V1,participant_id,label,sample_label,V_primer_mix,
         C_primer_barcode,V_barcode) %>% 
  left_join(C_barcode_file_2,by='C_primer_barcode') %>% 
  left_join(final_metadata %>% select(sample_label,Participant_PPID,specimen_label),by='sample_label') %>% 
  unique() %>% 
  mutate(Sample=paste0(participant_id,'-',sample_label,'-PBMC'))

sampleTable_3<-miseq_run3_corrected %>% 
  select(Sample,V_barcode,C_primer) %>% 
  dplyr::rename(SAMPLE1=V_barcode,
                SAMPLE2=C_primer) %>% 
  mutate(TagPattern='') %>% 
  dplyr::rename(SAMPLE1_old=SAMPLE1,
                SAMPLE2_old=SAMPLE2) %>% 
  dplyr::rename(SAMPLE1=SAMPLE2_old,
                SAMPLE2=SAMPLE1_old) %>% 
  select(Sample,TagPattern,SAMPLE1,SAMPLE2)
  write_tsv('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/sample_table_run3.tsv')
  
#barcodes issue
barcodes<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/PSD-TCRb_S1.txt')

left_join(barcodes,sampleTable_3,by=c('SAMPLE1','SAMPLE2')) %>% 
  group_by(SAMPLE1,SAMPLE2) %>% 
  summarise(samples=paste(Sample,collapse=','),
            counts=paste(count,collapse=',')) %>% View()
##metadata-----
metadata_3d_w0_w52<-dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/',pattern ='TRB_TRB.tsv') %>% 
  as.data.frame() %>% 
  dplyr::rename(Sample_long='.') %>% 
  mutate(Sample=str_remove(Sample_long,'result.'),
         Sample=str_remove(Sample,'.clones_TRB_TRB.tsv')) %>% 
  inner_join(miseq_run3_corrected %>% select(Participant_PPID,sample_label,specimen_label,Sample),
             by='Sample') %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/',Sample_long)) %>% 
  inner_join(final_metadata %>% select(week,specimen_label,Participant_PPID,Group3) %>% unique(),
            by=c('specimen_label','Participant_PPID'))

write_tsv(metadata_3d_w0_w52,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/metadata_run3.tsv')

##report----
report_run3<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/run3_report.tsv')
read_tsv(metadata_3d_w0_w52 %>% pull(real_path),id = 'real_path') %>% 
  group_by(real_path) %>% 
  summarise(n_clones=n(),
            n_reads=sum(readCount)) %>% View()

#CD63 high meta -----

cd63_values<-read_csv('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/POISED_BAT_CD63_Wks0_117.csv') %>% 
  dplyr::rename(group_cd63='...6')

cd63_values %<>% 
  mutate(rescaled_bl=case_when(`Wk 0`<=20~'Low',
                               `Wk 0`>=50~'High',
                               TRUE~'Inter')) %>%
  mutate(rescaled_w104=case_when(`Wk 104`<=20~'Low',
                               `Wk 104`>=50~'High',
                               is.na(`Wk 104`)~'No inf',
                               TRUE~'Inter'))
#Basic repertoire features----
##downsampling---- 
nClones_nReads_run3<-read_tsv(metadata_3d_w0_w52 |> 
                                pull(real_path),id = 'real_path') |> 
  group_by(real_path) |> 
  summarise(nclones=n(),
            nReads=sum(readCount))

write_tsv(nClones_nReads_run3,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/nCLones_nReads_run3_deepSeq.tsv')

nClones_nReads_run2<-read_tsv(second_run_files_deepSeq_new |> 
                                pull(real_path),id='real_path') |> 
  group_by(real_path) |> 
  summarise(nclones=n(),
            nReads=sum(readCount))

write_tsv(nClones_nReads_run2,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/nCLones_nReads_run2_deepSeq.tsv')


nClones_nReads_run1<-read_tsv(first_run_deep_seq_metadata |> 
                                pull(real_path),id='real_path') |> 
  group_by(real_path) |> 
  summarise(nclones=n(),
            nReads=sum(readCount))

write_tsv(nClones_nReads_run1,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/nCLones_nReads_run1_deepSeq.tsv')

nClones_nReads_allRuns<-nClones_nReads_run3 |> 
  bind_rows(nClones_nReads_run2) |> 
  bind_rows(nClones_nReads_run1)

write_tsv(nClones_nReads_allRuns,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/nCLones_nReads_all_runs_deepSeq.tsv')

nClones_nReads_allRuns_1m<-nClones_nReads_allRuns |> 
  filter(nReads>1000000)

##V usage----

v_family_usage_run1<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/postanalysis_1mln_reads_21oct/postanalysis.vFamilyUsage.TRB.tsv')
#v_family_usage_run2<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample/postanalysis/i.vFamilyUsage.TRB.tsv')
#v_family_usage_run3<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/downsample/postanalysis/i.vFamilyUsage.TRB.tsv')

vusage_2weeks<-function(data,week1,week2){
  Vusage_2weeks_data<-data %>% 
    filter(week_com==week1 | week_com==week2)  %>%
    group_by(Group,Participant_PPID,Vgene) %>% 
    mutate(n_weeks_present=n()) %>%
    filter(n_weeks_present>1) %>% #filter patients with at least 2 tp
    ungroup() %>% 
    complete(Participant_PPID,Vgene,week_com,fill = list(mean_freq_by_v_gene=0)) %>% 
    mutate(OIT=case_when(Group=='C'~'Placebo',TRUE~'OIT'))
  return(Vusage_2weeks_data)
}

###week 0-w104 placebo vs OIT----
Vusage_w104_w117_initial<-v_family_usage_run1 %>%
  select(-TRBV17,-TRBV1,-TRBV3) %>% #to select low representative genes
  pivot_longer(cols = !c('sample'),
               names_to = 'Vgene',values_to = 'family_freq' ) |> 
  mutate(family_freq=ifelse(is.na(family_freq),0,family_freq)) %>% 
  mutate(sample=str_replace(sample,'TRB.result.*','clns.tsv')) %>% 
  mutate(sample=str_remove(sample,'result.')) |> 
  dplyr::rename(sample_name=sample) %>% 
  inner_join(deep_seq_run1_metadata,by='sample_name') %>%
  mutate(week_com=case_when(week=='W104A'~'W104', week=='W104B'~'W104',
                            week=='W117A'~'W117', week=='W117B'~'W117', TRUE~week)) %>% 
  group_by(Group,Participant_PPID,Vgene,week_com) %>% #to find mean freq because of replicates 104a 104b/117a 117b
  summarise(mean_freq_by_v_gene=mean(family_freq))

Vusage_104_w117<-vusage_2weeks(data=Vusage_w104_w117_initial,week1 = 'W104',week2 = 'W117') |> 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup %>% select(Group_dose_outcome,Participant_PPID) %>% unique(),by='Participant_PPID') %>% 
  filter(Group_dose_outcome=='Pea_0_fa'|Group_dose_outcome=='Pea_0_suc')

vusage_w104_w117_stat_test<-Vusage_104_w117 %>%  
  filter(week_com=='W117') |> 
  group_by(Vgene) %>% 
  wilcox_test(mean_freq_by_v_gene ~ Group_dose_outcome,comparisons = list(c('Pea_0_fa','Pea_0_suc'))) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

Vusage_104_w117 %>%   
  filter(week_com=='W117') |> 
  ggplot(aes(x=Group_dose_outcome,y=mean_freq_by_v_gene,color=Group_dose_outcome))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=2,alpha=0.5)+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  #geom_line(aes(group=Participant_PPID),alpha=0.3)+
  facet_wrap(~Vgene,scales = 'free_y',ncol = 6)+
  stat_pvalue_manual(data=vusage_w104_w117_stat_test,label = 'p.adj',
                     y.position = c(0.065, 0.05, 0.11588, 0.005, 0.025, 0.03, 0.00076, 0.07688, 
                                     0.07788, 0.06788, 0.14088, 0.05688, 0.04088, 0.03888, 0.10288, 0.11788, 0.06988, 0.06788,
                                     0.17088,0.15888, 0.25088, 0.05388))+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  theme(strip.text.x = element_text(size = 15))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/v_usage_w117_su_vs_ds.png',
       width =12,height = 7.5,device = 'png',dpi = 400)  

Vusage_initial_bl_104 %>%   #graph by OIT and Placebo baseline vs w104
  ggplot(aes(x=OIT,y=mean_freq_by_v_gene))+
  geom_point(alpha=0.5)+
  geom_boxplot()+
  facet_wrap(~Vgene+week_com,scales = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means()+
  theme_bw()

###week 117 DS vs SU----
Vusage_initial_w104_w117<-vusage_2weeks(data=Vusage_initial,week1 = 'W104',week2 = 'W117') %>% 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup %>% select(Group_dose_outcome,Participant_PPID) %>% unique(),by='Participant_PPID') %>% 
  filter(Group_dose_outcome=='Pea_0_fa'|Group_dose_outcome=='Pea_0_suc')

#graph w 117 SU and DS
Vusage_initial_w104_w117 %>% 
  filter(Vgene!='TRBV2' & Vgene!='TRBV20' &  Vgene!='TRBV18' &  Vgene!='TRBV27'  &  Vgene!='TRBV10'
         & Vgene!='TRBV12'  & Vgene!='TRBV13'  & Vgene!='TRBV14'  & Vgene!='TRBV19'
         & Vgene!='TRBV24'  & Vgene!='TRBV25'  & Vgene!='TRBV28'  & Vgene!='TRBV6') %>% 
  ggplot(aes(x=week_com,y=mean_freq_by_v_gene))+
  geom_point()+
  geom_line(aes(group=Participant_PPID),alpha=0.4)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Vgene+Group_dose_outcome,scales = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means()+
  labs(title='V gene usage between SU and DS week 104- 117 on downsampled data',y='V gene family freq ', x='SU vs DS ')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/v_usage_su_ds_w104_w117_paired.png',
       width =12,height = 10,device = 'png',dpi = 400)  


##CDR3 metrics----
cdr3metrics<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/postanalysis/i.cdr3metrics.TRB.tsv')

cdr3metrics_right<- cdr3metrics %>% 
  mutate(sample=str_replace(sample,'.clns','.clns.tsv')) %>% 
  mutate(sample=str_remove(sample,'result.')) %>% 
  dplyr::rename(sample_name=sample) %>% 
  pivot_longer(cols = !c('sample_name'),
               names_to = 'cdr3_metrics',values_to = 'value' )

cdr3metrics_w104_w117<- cdr3metrics_right %>% 
  inner_join(deep_seq_run1_metadata,by='sample_name') %>% 
  unique() %>% 
  mutate(week_com=case_when(week=='W104A'~'W104',
                            week=='W104B'~'W104',
                            week=='W117A'~'W117',
                            week=='W117B'~'W117',
                            TRUE~week)) %>% 
  group_by(week117_FC_result,Group,Participant_PPID,cdr3_metrics,week_com) %>% 
  summarise(mean_value=mean(value)) %>% 
  ungroup()

#graph
cdr3metrics_w104_w117 %>% 
  filter(week_com=='W104')  %>% 
  mutate(OIT=case_when(Group=='C'~'Placebo',
                       TRUE~'OIT')) %>% 
  filter(cdr3_metrics=='Added nucleotides'|
           cdr3_metrics=='Charge of CDR3'|
           cdr3_metrics=='Length of CDR3, aa'|
           cdr3_metrics=='Length of CDR3, nt'|
           cdr3_metrics=='N2Hydrophobicity of CDR3') %>% 
  ggplot(aes(x=OIT,y=mean_value,group=OIT))+
  geom_point()+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~cdr3_metrics,scales = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means()+
  labs(title='CDR3 metrics between placebo and OIT ',
       y='value ',
       x='Status')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/cdr3_metrics.png',
       width =7,height = 5,device = 'png',dpi = 400)  

##diversity metrics----
#nclones_nReads_func(list_samples = dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/postanalysis_1mln_reads_21oct/',
    #pattern='1mil_reads_TRB.tsv',full.names = TRUE)) |> View()

diversity_metrics<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/postanalysis_1mln_reads_21oct/postanalysis.diversity.TRB.tsv')

diversity_metrics_right<- diversity_metrics %>% 
  mutate(sample=str_replace(sample,'TRB.result.*','clns.tsv')) %>% 
  mutate(sample=str_remove(sample,'result.')) |> 
  dplyr::rename(sample_name=sample) %>% 
  pivot_longer(cols = !c('sample_name'),
               names_to = 'diversity_metrics',values_to = 'value' )

diversity_metrics_w104_w117<- diversity_metrics_right %>% 
  inner_join(deep_seq_run1_metadata,by='sample_name') %>% 
  unique() %>% 
  mutate(week_com=case_when(week=='W104A'~'W104',
                            week=='W104B'~'W104',
                            week=='W117A'~'W117',
                            week=='W117B'~'W117',
                            TRUE~week)) %>% 
  group_by(week117_FC_result,Group,Participant_PPID,diversity_metrics,week_com) %>% 
  summarise(mean_value=mean(value)) %>% 
  ungroup()

#graph BL/w0 placebo and OIT

more_than800k_reads_run2<-nclones_nReads_func(list_samples = dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample_1mln_reads/',
  pattern='1mil_reads_TRB.tsv',full.names = TRUE)) |>
  filter(nReads>800000)

more_than900k_reads_run3<-nclones_nReads_func(list_samples = dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/downsample_1mln_reads',
                                                                 pattern='1mil_reads_TRB.tsv',full.names = TRUE)) |> 
  filter(nReads>900000)

diversity_run3<-fread('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/downsample_1mln_reads/postanalysis.diversity.TRB.tsv') |> 
  mutate(sample=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run3/downsample_1mln_reads/',sample)) |> 
  mutate(sample=str_replace(sample,'.clns','_TRB.tsv')) |> 
  filter(sample %in% more_than900k_reads_run3$real_path)
  
diversity_w0<-diversity_run3 |> 
  mutate(Sample=str_remove(sample,'.TRB.down_1mil_reads_TRB.tsv'),
         Sample=str_remove(Sample,'\\/.*\\/result.')) |>
  inner_join(metadata_3d_w0_w52 |> select(Sample,week,Participant_PPID),by='Sample') |> 
  filter(week!='W052')
  
  
#graph w0-w104 placebo vs OIT
w0_w104_diversity<-diversity_w0 |> 
  select(-sample) |> 
  pivot_longer(cols =!c('Sample','week','Participant_PPID') ,
               names_to = 'diversity_metrics',values_to = 'mean_value') |> 
  dplyr::rename(week_com=week) |> 
  bind_rows(diversity_metrics_w104_w117) |> 
  select(Participant_PPID,diversity_metrics,week_com,mean_value) |> 
  inner_join(final_metadata |> select(Group,Participant_PPID) |> unique(),by='Participant_PPID') |> 
  mutate(OIT=case_when(Group=='C'~'Placebo',
                         TRUE~'OIT')) %>% 
  filter(week_com=='W000'| week_com=='W104') |> 
  group_by(Participant_PPID,OIT,diversity_metrics) |> 
  mutate(n_tp=n()) |> 
  filter(n_tp==2) |> 
  filter(diversity_metrics=='Chao1 estimate' |  diversity_metrics=='Shannon-Wiener diversity' |
           diversity_metrics=='Observed diversity')

w0_w104_diversity |> 
  ggplot(aes(x=week_com,y=mean_value))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID,color=OIT),alpha=0.3)+
  geom_point(alpha=0.3,aes(color=OIT))+
  scale_color_manual(values=c('#E69F00',"#0072B2"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means(paired = TRUE)+
  facet_wrap(~diversity_metrics+OIT,scales = 'free_y',nrow = 3)+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 15))+
  labs(title='Div. metrics between week 0 and week 104 in pl. and OIT',
       y='value ')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/diversity_metrics_w0_w104_pl_oit.png',
       width =6,height = 8,device = 'png',dpi = 400) 

#graph w104 placebo and OIT
diversity_metrics_w104_w117 %>% 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> 
               select(Participant_PPID,Group_dose_outcome)|> unique(),
             by='Participant_PPID') |> 
  filter(week_com=='W104')  %>% 
  mutate(OIT=case_when(Group=='C'~'Placebo',
                       TRUE~'OIT')) %>% 
  filter(diversity_metrics=='Chao1 estimate' | 
           diversity_metrics=='Shannon-Wiener diversity' |
           diversity_metrics=='Observed diversity') %>% 
  mutate(diversity_metrics=case_when(diversity_metrics=='Shannon-Wiener diversity'~'ShW diversity',TRUE~diversity_metrics)) |> 
  ggplot(aes(x=OIT,y=mean_value,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=2,alpha=0.8)+
  scale_color_manual(values=c('#E69F00',"#0072B2"))+
  facet_wrap(~diversity_metrics,scales = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means()+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 15))+
  labs(title='Diversity metrics between placebo and OIT at week104 ',y='value ')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/diversity_metrics_placebo_vs_OIT_w104.png',
       width =8,height = 4,device = 'png',dpi = 400) 

#w 117 su vs ds

diversity_metrics_w104_w117 %>% 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Participant_PPID,Group_dose_outcome)|> unique(),by='Participant_PPID') |> 
  filter(week_com=='W117')  %>% 
  filter(Group_dose_outcome=='Pea_0_fa' | Group_dose_outcome=='Pea_0_suc' ) |>
  filter(diversity_metrics=='Chao1 estimate' | 
           diversity_metrics=='Normalized Shannon-Wiener index' | 
           diversity_metrics=='Shannon-Wiener diversity' |
          diversity_metrics=='Observed diversity' ) |> 
  mutate(diversity_metrics=case_when(diversity_metrics=='Shannon-Wiener diversity'~'ShW diversity',TRUE~diversity_metrics)) |> 
  ggplot(aes(x=Group_dose_outcome,y=mean_value,color=Group_dose_outcome))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=2,alpha=0.9)+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  facet_wrap(~diversity_metrics,scales = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  stat_compare_means()+
  labs(title='Diversity metrics between SU and DS at week 117 ', y='Value ')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 15))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/diversity_metrics_w117_su_vs_ds.png',
       width =8,height = 4,device = 'png',dpi = 400) 


#summed freq of PR clones in bulk pbmc----
cd4_seq_corrected_only_trb %>% 
  mutate(cell_subset=ifelse(cell_subset=='PR_TNa_2','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_3','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_TNa_4','PR_TNa_1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_Tfh13','PR_Th2conv',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Treg_act2','PR_Treg_act1',cell_subset),
         cell_subset=ifelse(cell_subset=='PR_Teff_Me_act2','PR_Teff_Me_act1',cell_subset))

process_mixcr_clones<-function(list_samples){
  data_clones_mixcr<-read_tsv(list_samples,id = 'real_path' ) %>% 
    dplyr::rename(v_trb=allVHitsWithScore, j_trb=allJHitsWithScore,
                  cdr3aa_trb=aaSeqCDR3,cdr3nt_trb=nSeqCDR3,
                  count=readCount, freq=readFraction) %>% 
    select(real_path,v_trb,j_trb,cdr3aa_trb,cdr3nt_trb,count,freq) %>% 
    mutate(v_trb=str_remove(v_trb,'\\*\\d+.*')) %>% 
    mutate(j_trb=str_remove(j_trb,'\\*\\d+.*')) %>% 
    mutate(vj_trb=paste0(v_trb,'_',j_trb))
  return(data_clones_mixcr)
}

grouping_mixcr_clones<-function(file_mixcr_clones){
  file_mixcr_clones %>% 
  group_by(real_path,cdr3aa_trb,v_trb,j_trb) %>% 
    summarise(count_bulk=sum(count),
              freq_bulk=sum(freq),
              n_cdr3nt_bulk=n_distinct(cdr3nt_trb)) %>% ungroup() %>% 
    mutate(vj_trb=paste0(v_trb,'_',j_trb))
}

make_mean_clonesFreq_forReplicates<-function(data_full_repertoire,metadata_bulk_data,week_common,week_2,week_3){
  data_full_rep_metadata<-data_full_repertoire %>%
   inner_join(metadata_bulk_data %>% select(Participant_PPID,real_path,week) %>% 
              unique(),by='real_path')
  
  print(nrow(data_full_rep_metadata))
  data_replicated_mean<-data_full_rep_metadata %>% 
    select(cdr3aa_trb,vj_trb,Participant_PPID,week,freq_bulk) %>% 
    filter(week==week_2| week==week_3) %>% 
    pivot_wider(names_from =week, values_from = freq_bulk,
                id_cols = c('cdr3aa_trb','vj_trb','Participant_PPID')) %>%
    mutate({{week_2}}:=ifelse(is.na(get(!!week_2)),0,get(!!week_2)),
           {{week_3}}:=ifelse(is.na(get(!!week_3)),0,get(!!week_3))) %>% 
    mutate(freq_bulk=rowMeans(select(.,c((!!week_2),(!!week_3)))))
  
  data_corrected_freq<- data_full_rep_metadata %>% 
    select(cdr3aa_trb,vj_trb,Participant_PPID,week,freq_bulk) %>% 
    filter(week==week_common) %>% 
    bind_rows(data_replicated_mean) %>% 
    mutate(week=week_common)
  return(data_corrected_freq) 
}

find_sc_clones_in_bulk<-function(data_bulk,data_sc,week_bulk,week_sc,metadata_bulk_data){
  sc_clones_in_bulk<-data_bulk %>% 
    inner_join(metadata_bulk_data %>% select(Participant_PPID,real_path,week) %>% unique(),by='real_path') %>% 
    inner_join(data_sc |>  
                 #filter(timepoint_single_cell!=week_sc_added) |> 
                 filter(timepoint_single_cell==week_sc),
               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
    select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
    unique() %>% 
    group_by(Participant_PPID,OIT,week) %>% 
    summarise(sum_freq_in_bulk=sum(freq_bulk)) %>% 
    ungroup() %>% 
    right_join(metadata_bulk_data %>% 
                 filter(week==week_bulk) %>% 
                 #filter(downsample_file_present==TRUE) %>% 
                 filter(Participant_PPID %in% (data_sc %>% 
                                                 filter(timepoint_single_cell==week_sc) %>% pull(Participant_PPID))) %>% 
                 select(Participant_PPID,Group3),
               by='Participant_PPID') %>% 
    mutate(OIT=case_when(Group3=='Placebo'~'Placebo',
                         Group3=='Peanut 0'~'Active',
                         Group3=='Peanut 300'~'Active',
                         TRUE~OIT)) %>% 
    mutate(sum_freq_in_bulk=ifelse(is.na(sum_freq_in_bulk),0,sum_freq_in_bulk)) %>% unique() %>% 
    mutate(week=ifelse(is.na(week),week_bulk,week)) %>% 
    mutate(week_common=week_bulk)
  sc_clones_in_bulk_mean<-sc_clones_in_bulk %>% 
    group_by(week_common,Participant_PPID,OIT,Group3) %>% 
    mutate(mean_freq_in_bulk=mean(sum_freq_in_bulk))
  return(sc_clones_in_bulk_mean)
}

##BL PR baseline vs week 0-----
###total----
w0_not_downsampled<-process_mixcr_clones(list_samples = metadata_3d_w0_w52 %>% filter(week=='W000') %>% 
                                                  pull(real_path))

w0_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones = w0_not_downsampled)

prBL_in_w0_not_downsampled<-find_sc_clones_in_bulk(data_bulk= w0_not_downsampled,  data_sc = cd4_seq_corrected_only_trb, 
                                                            week_sc = 'BL',week_bulk = 'W000',
                                                            metadata_bulk_data = metadata_3d_w0_w52)

pr_bl_bulk_bl_w0_not_downsampled<-prBL_in_w0_not_downsampled %>% 
  bind_rows(prCL_BL_in_bulk_BL_not_downsample) %>% 
  select(-week,-mean_freq_in_bulk) %>% unique() %>% 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  mutate(week_common=ifelse(week_common=='BL01','baseline','week 0'))

united_w0_and_bl %>% 
  mutate(week_common=ifelse(week_common=='BL01','baseline','week 0')) |> 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  mutate(week_common=factor(x=week_common,levels=c('baseline','week 0'))) %>% 
  mutate(oit_tp=paste0(OIT,'_',week_common)) %>% 
    ggplot(aes(x=week_common,y=mean_freq_in_bulk_normalized))+
    geom_boxplot(outlier.shape = NA,aes(color=oit_tp))+
    geom_point(alpha=0.7,aes(color=oit_tp))+
    geom_line(aes(group=Participant_PPID),alpha=0.8,size=0.3, colour="grey")+
    scale_color_manual(values=c("#D55E00", "#E69F00","#D55E00","#0072B2"))+
    stat_compare_means(paired = TRUE,method ='wilcox.test')+
    facet_wrap(~OIT,scales = 'free_y',)+
    labs( y='Normalized summed frequency (by patient ID)', x='Week', fill='Week')+
    theme_classic()+
    theme(strip.background = element_blank())+
    xlab(element_blank())+
    theme(legend.position="none")+ 
    theme(strip.text.x = element_text(size = 15))
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/BL_PRclones_in_w0_bl_bulk_norm.png',
       width = 6,height = 4,device = 'png',dpi = 400) 

###by subset----
pr_cl_bL_w0_bulk_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= w0_not_downsampled %>% 
                                                            inner_join(metadata_3d_w0_w52 %>%
                                                                         select(real_path,Participant_PPID,week) %>% unique(),by='real_path'), data_sc = cd4_seq_corrected_only_trb, 
                                                       week_sc = 'BL',week_bulk = 'W000')

pr_bl_bl_w0_bulk_by_subset<-pr_cl_bL_bulk_by_subset %>% 
  #dplyr::rename(week=week_common) |> 
  bind_rows(pr_cl_bL_w0_bulk_by_subset) %>% 
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>%
  group_by(Participant_PPID,OIT,cell_subset) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  dplyr::rename(week=week_common) |> 
  left_join(nClones_ncReads_realpath |> 
              mutate(week=case_when((Participant_PPID=='P034' & week=='BL02')~'BL01',
                                    (Participant_PPID=='P074' & week=='BL02')~'BL01',
                                    TRUE~week)),by=c('Participant_PPID','week')) |> 
  mutate(oit_tp=paste0(OIT,'_',week)) |> 
  mutate(sum_freq_in_bulk_by_subset_norm=sum_freq_in_bulk_by_subset/nclones)

stat_pr_bl_bl_w0_bulk_by_subset<-pr_bl_bl_w0_bulk_by_subset %>% 
  group_by(cell_subset,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,comparisons =c('BL01','W000')) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() #all pvalues are not significant

pr_bl_bl_w0_bulk_by_subset %>% 
ggplot(aes(x=week,y=sum_freq_in_bulk_by_subset_norm))+
  geom_boxplot(outlier.shape = NA,aes(color=oit_tp))+
  geom_point(alpha=0.7,aes(color=oit_tp))+
  geom_line(aes(group=Participant_PPID),alpha=0.8,size=0.3, colour="grey")+
  scale_color_manual(values=c("#D55E00", "#E69F00","#D55E00","#0072B2"))+
  facet_wrap(~cell_subset+OIT,scales = 'free',nrow =  3 )+
  scale_y_log10(expand = expansion(mult = c(0.1, 0.1)), oob = scales::squish_infinite)+
  labs( y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/Subsets_BL_PR_inBL_W0_bulk_norm.png',
       width = 16,height = 9,device = 'png',dpi = 200)

##PR BL in W0,W52,W104 bulk data----
###total----
#make w104 NOT downsampled data
w104_bulk_not_downsampled<-process_mixcr_clones(list_samples = first_run_deep_seq_metadata %>% filter(week=='W104'| week=='W104B'| week=='W104A') %>% 
                       pull(real_path))

w104_bulk_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones = w104_bulk_not_downsampled)

prCL_BL_in_bulk_w104_not_downsample<-find_sc_clones_in_bulk(data_bulk= w104_bulk_not_downsampled,  data_sc = cd4_seq_corrected_only_trb, 
                                        week_sc = 'BL',week_bulk = 'W104',
                                        metadata_bulk_data = first_run_deep_seq_metadata %>% 
                                          mutate(week_bulk_samples=case_when(week=='W104'~'W104',
                                                                             week=='W104A'~'W104',
                                                                             week=='W104B'~'W104',TRUE~week)))


cd4_seq_corrected_only_trb<-сd4_seq_corrected_only_trb

#make baseline NOT downsampled data
bl01_repertoire_not_downsampled<-process_mixcr_clones(list_samples=second_run_files_deepSeq_new %>%  filter(week=='BL01' | week=='BL02' | week=='BL03' ) %>% 
                                      filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID) %>% pull(real_path))

bl01_repertoire_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones =bl01_repertoire_not_downsampled )

prCL_BL_in_bulk_BL_not_downsample<-find_sc_clones_in_bulk(data_bulk= bl01_repertoire_not_downsampled,
                                                          data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'BL01',
                                      metadata_bulk_data = second_run_files_deepSeq_new %>% 
                                        mutate(week_bulk_samples=case_when(week=='BL01'~'BL01',
                                                                           week=='BL02'~'BL01',
                                                                           week=='BL03'~'BL01',TRUE~week)))
#to calculate for w36
prCL_BL_in_bulk_w36_not_downsample<-find_sc_clones_in_bulk(data_bulk= w36_bulk_not_downsampled,
                                                           data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'W36',
                                                           metadata_bulk_data=second_run_files_deepSeq_new %>% 
                                                             mutate(week=case_when(week=='W34'~'W36',
                                                                                   week=='W38'~'W36',
                                                                                   TRUE~week)))

#to calculate for w52
prCL_BL_in_bulk_w52_not_downsample<-find_sc_clones_in_bulk(data_bulk= w52_bulk_not_downsampled,
                                                          data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'W052',
                                                          metadata_bulk_data =metadata_3d_w0_w52)


#normalization
nClones_ncReads_realpath<-nClones_nReads_allRuns |> 
  left_join(metadata_3d_w0_w52 |> 
              bind_rows(second_run_files_deepSeq_new) |> 
              bind_rows(first_run_deep_seq_metadata),by='real_path') |> 
  select(real_path,nclones,nReads,Participant_PPID,week)


#to unite bl and w0
united_w0_and_bl<-prBL_in_w0_not_downsampled %>% 
  left_join(nClones_ncReads_realpath,by=c('Participant_PPID','week')) |> #to add inf about normalization
  bind_rows(prCL_BL_in_bulk_BL_not_downsample |> 
              left_join(nClones_ncReads_realpath,by=c('Participant_PPID','week'))) %>% 
  mutate(week_common2='W000') |> 
  group_by(Participant_PPID,OIT,Group3,week_common2) |> 
  mutate(mean_freq_in_bulk_2=mean(mean_freq_in_bulk),
            nclones_mean=mean(nclones)) |> 
  mutate(mean_freq_in_bulk_normalized=mean_freq_in_bulk_2/nclones_mean)
  

pr_bl_bulk_w0_w52_w104_not_downsampled<-united_w0_and_bl %>% 
  dplyr::rename(mean_freq_in_bulk=mean_freq_in_bulk_2) |> 
  bind_rows(prCL_BL_in_bulk_w52_not_downsample |> 
              left_join(nClones_ncReads_realpath,by=c('Participant_PPID','week')) |> 
              mutate(mean_freq_in_bulk_normalized=mean_freq_in_bulk/nclones)) |> 
  bind_rows(prCL_BL_in_bulk_w104_not_downsample |> 
              left_join(nClones_ncReads_realpath,by=c('Participant_PPID','week')) |> 
              group_by(Participant_PPID,OIT,week_common) |> 
              mutate(mean_nclones=mean(nclones),
                        mean_freq_in_bulk_normalized=mean_freq_in_bulk/mean_nclones)) |> 
  mutate(week_common=ifelse(is.na(week_common),week_common2,week_common)) |> 
  select(Participant_PPID,OIT,mean_freq_in_bulk,week_common) |> 
 # select(Participant_PPID,OIT,mean_freq_in_bulk_normalized,week_common) |>  # for normalization
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==3) 


stat_pr_bl_w0_w52_w104_bulk<-pr_bl_bulk_w0_w52_w104_not_downsampled %>% ungroup() |> 
  group_by(OIT) %>% 
  wilcox_test(mean_freq_in_bulk ~ week_common, paired = T,comparisons = list(c('W052','W000'),
                                                                             c('W104','W052'),
                                                                             c('W000','W104'))) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

pr_bl_bulk_w0_w52_w104_not_downsampled %>% ungroup() |> 
  mutate(week_common=factor(x=week_common,labels = c('W000','W052','W104'))) %>% 
  ggplot(aes(x=week_common,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,aes(color=OIT))+
  scale_color_manual(values=c('#E69F00','#0072B2'))+
  stat_pvalue_manual(data =stat_pr_bl_w0_w52_w104_bulk,label = 'p.adj',
                     y.position =c(0.025,0.028) ,hide.ns = TRUE,)+
  facet_wrap(~OIT,scales = 'free')+
  labs(title='Summed freq of BL PR clones in w0, w52 and w104 bulk ',
       y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PRclones_Sumfreq_Bl_in_w0_w52_w104_notDownsampled.png',
       width = 6,height = 4,device = 'png',dpi = 400) 

###T eff vs Tregs PR bl in bl, w0, w52, w104----
find_sc_clones_in_bulk_by_large_subset<-function(data_bulk,data_sc,week_bulk,week_sc){
  
  sc_clones_in_bulk<-data_bulk %>% 
    right_join(data_sc %>% filter(timepoint_single_cell== week_sc),
               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
    unique() %>% 
    select(cell_subset,Cell_type,cdr3aa_trb,vj_trb,Participant_PPID,freq_bulk) %>% 
    unique() %>% 
    mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) %>% 
    group_by(Participant_PPID,Cell_type) %>% 
    summarise(sum_freq_in_bulk_by_subset=sum(freq_bulk)) %>% 
    ungroup() %>% 
    mutate(sum_freq_in_bulk_by_subset=ifelse(is.na(sum_freq_in_bulk_by_subset),0,sum_freq_in_bulk_by_subset)) %>% unique() %>% 
    mutate(week_common=week_bulk)
  
  sc_clones_in_bulk_filtered<-sc_clones_in_bulk %>% 
    filter(Participant_PPID %in% (data_bulk %>% pull(Participant_PPID) %>% unique()))
  
  #  sc_clones_in_bulk_mean<-sc_clones_in_bulk %>% 
  #  group_by(week_common,Participant_PPID,OIT,Group3,cell_subset) %>% 
  #  mutate(mean_freq_in_bulk_by_subset=mean(sum_freq_in_bulk_by_subset))
  return(sc_clones_in_bulk_filtered)
}

#bl
pr_bl_in_bulk_BL_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= bl01_repertoire_not_downsampled |> 
                                                                                 inner_join(second_run_files_deepSeq_new %>% 
                                                                                              select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                                                                               data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'BL01')

#w0
pr_bl_in_bulk_w0_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w0_not_downsampled |> 
                                                                                 inner_join(metadata_3d_w0_w52 %>% 
                                                                                              select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                                                                               data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'W000')

#w52
pr_bl_in_bulk_w52_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w52_bulk_not_downsampled |> 
                                                                                  inner_join(metadata_3d_w0_w52 %>% 
                                                                                               select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                                                                                data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'W052')

#w104
pr_bl_in_bulk_w104_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled ,
                                                                                 data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'W104')
#unite bl and week 0
pr_bl_w0_bl_not_down_large_sub<-pr_bl_in_bulk_BL_not_down_large_subset %>% 
  bind_rows(pr_bl_in_bulk_w0_not_down_large_subset) %>%
  mutate(week_common2="W000") |> 
  group_by(Participant_PPID,Cell_type,week_common2) |> 
  summarise(mean_bl_w0_freq_bulk=mean(sum_freq_in_bulk_by_subset)) |> 
  dplyr::rename(week_common=week_common2,
                sum_freq_in_bulk_by_subset=mean_bl_w0_freq_bulk)

#for normalization
nclones_nReads_func<-function(list_samples){
  data_clones_mixcr<-read_tsv(list_samples,id = 'real_path' )
  data_stat<-data_clones_mixcr |> group_by(real_path) |> 
    summarise(nClones=n(),
              nReads=sum(readCount))
  return(data_stat)
}

pr_bl_w0_w52_w104_not_downsampled_large_subsets<-pr_bl_w0_bl_not_down_large_sub %>% 
  bind_rows(pr_bl_in_bulk_w52_not_down_large_subset) %>% 
  bind_rows(pr_bl_in_bulk_w104_not_down_large_subset) |> 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Participant_PPID,OIT) |> unique(),
             by='Participant_PPID') |> 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  group_by(Participant_PPID,OIT,Cell_type) %>% 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==3)

stat_pr_0_bl_w0_w52_w104_not_downsampled_large_subsets<-pr_bl_w0_w52_w104_not_downsampled_large_subsets %>% 
  mutate(week_common=ifelse(week_common=='W000','baseline/W000',week_common)) |> 
  group_by(OIT,Cell_type) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,comparisons = list(c('W052','baseline/W000'),
                                                                                      c('W104','W052'),
                                                                                      c('baseline/W000','W104'))) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

pr_bl_w0_w52_w104_not_downsampled_large_subsets %>% 
  mutate(week_common=ifelse(week_common=='W000','baseline/W000',week_common)) |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7)+
  scale_color_manual(values=c('#E69F00',"#0072B2"))+
  stat_pvalue_manual(data=stat_pr_0_bl_w0_w52_w104_not_downsampled_large_subsets,
                     label ='p.adj.signif',hide.ns = TRUE, y.position=c(0.022))+
  facet_wrap(~OIT+Cell_type,scales = 'free')+
  labs( y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PRclones_Sumfreq_by_large_sub_bl_in_w0_w52_w104_not_downsampled.png',
       width = 7.5,height = 5,device = 'png',dpi = 400) 

###Th1/Th2a/Th2conv----
find_sc_clones_in_bulk_by_subset<-function(data_bulk,data_sc,week_bulk,week_sc_add,week_sc){
  
  sc_clones_in_bulk<-data_bulk %>% 
    right_join(data_sc %>% 
              #  filter(timepoint_single_cell!=week_sc_add) |> 
                 filter(timepoint_single_cell== week_sc),
               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
    unique() %>% 
    select(cell_subset,cdr3aa_trb,vj_trb,Participant_PPID,freq_bulk) %>% 
    unique() %>% 
    mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) %>% 
    group_by(Participant_PPID,cell_subset) %>% 
    summarise(sum_freq_in_bulk_by_subset=sum(freq_bulk)) %>% 
    ungroup() %>% 
    mutate(sum_freq_in_bulk_by_subset=ifelse(is.na(sum_freq_in_bulk_by_subset),0,sum_freq_in_bulk_by_subset)) %>% unique() %>% 
    mutate(week_common=week_bulk)
  
  sc_clones_in_bulk_filtered<-sc_clones_in_bulk %>% 
    filter(Participant_PPID %in% (data_bulk %>% pull(Participant_PPID) %>% unique()))
  
  #  sc_clones_in_bulk_mean<-sc_clones_in_bulk %>% 
  #  group_by(week_common,Participant_PPID,OIT,Group3,cell_subset) %>% 
  #  mutate(mean_freq_in_bulk_by_subset=mean(sum_freq_in_bulk_by_subset))
  return(sc_clones_in_bulk_filtered)
}


pr_cl_bL_bulk_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= bl01_repertoire_not_downsampled %>% 
                                                            inner_join(second_run_files_deepSeq_new %>% 
                                                                         select(real_path,Participant_PPID,week) %>% unique(),by='real_path'), 
                                                          data_sc = cd4_seq_corrected_only_trb, 
                                                          week_sc = 'BL',week_bulk = 'BL01')

pr_cl_bL_bulk_w052_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= w52_bulk_not_downsampled %>% 
                                                                 inner_join(metadata_3d_w0_w52 %>% select(real_path,week,Participant_PPID),
                                                                            by='real_path'),
                                                               data_sc = cd4_seq_corrected_only_trb, 
                                                               week_sc = 'BL',week_bulk = 'W052')

#w104
#preparing metadata
first_run_deep_seq_metadata %<>% 
  mutate(week_bulk_samples=case_when(week=='W104'~'W104', week=='W104A'~'W104',week=='W104B'~'W104',TRUE~week))
first_run_deep_seq_metadata %<>% 
  mutate(OIT=case_when(Group3=='Placebo'~'Placebo', Group3=='Peanut 0'~'Active',Group3=='Peanut 300'~'Active'))

#because w117 and w104 have replicates 

#make new freq based on replicates
w104_bulk_corrected_freq_not_downsampled<-make_mean_clonesFreq_forReplicates(data_full_repertoire=w104_bulk_not_downsampled,
                                                                             week_common = 'W104',week_2 = 'W104A',week_3 = 'W104B',
                                                                             metadata_bulk_data = first_run_deep_seq_metadata)

#clones calculations  
pr_cl_BL_bulk_w104<-find_sc_clones_in_bulk_by_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled, 
                                                     data_sc = cd4_seq_corrected_only_trb,
                                                     week_sc = 'BL',week_bulk = 'W104')
#Bl and w0 uniting
pr_bl_united_Bl_w0<-pr_cl_bL_w0_bulk_by_subset %>% #by subset by pat_id bl
  bind_rows(pr_cl_bL_bulk_by_subset) |> 
  mutate(week_common2='W000') |> 
  group_by(Participant_PPID,cell_subset,week_common2) |> 
  summarise(bl_w0_mean_greq_bulk=mean(sum_freq_in_bulk_by_subset)) |>
  dplyr::rename(week_common=week_common2,
                sum_freq_in_bulk_by_subset=bl_w0_mean_greq_bulk)

pr_bl_united_w0_w52_w104_bulk_subsets<-pr_bl_united_Bl_w0 |> 
  bind_rows(pr_cl_bL_bulk_w052_by_subset) %>% 
  bind_rows(pr_cl_BL_bulk_w104) %>% 
  filter(cell_subset=='PR_Th2conv'|
           cell_subset=='PR_Teff_Me_Th1_CTL'|
           cell_subset=='PR_Teff_Me_Th2a') |> 
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>% 
  group_by(Participant_PPID,OIT,cell_subset) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==3) |> 
  filter(Participant_PPID %in% (c((w104_bulk_corrected_freq_not_downsampled |> pull(Participant_PPID) |> unique()),
                                  (cd4_seq_corrected_only_trb |> filter(timepoint_single_cell=='BL') |> pull(Participant_PPID) |> unique())))) |> 
  filter(Participant_PPID %in% (c((metadata_3d_w0_w52 |> filter(week=='W000') |> pull(Participant_PPID) |> unique()),
                                  (cd4_seq_corrected_only_trb |> filter(timepoint_single_cell=='BL') |> pull(Participant_PPID) |> unique())))) |> 
  filter(Participant_PPID %in% (c((metadata_3d_w0_w52 |> filter(week=='W052') |> pull(Participant_PPID) |> unique()),
                                  (cd4_seq_corrected_only_trb |> filter(timepoint_single_cell=='BL') |> pull(Participant_PPID) |> unique())))) |> 
  filter(Participant_PPID %in% (c((second_run_files_deepSeq_new |> filter(week=='BL01'|week=='BL02'| week=='BL03') |> pull(Participant_PPID) |> unique()),
                                  (cd4_seq_corrected_only_trb |> filter(timepoint_single_cell=='BL') |> pull(Participant_PPID) |> unique())))) 

stat_pr_bl_w0_w52_w104_bulk_subsets<-pr_bl_united_w0_w52_w104_bulk_subsets %>% 
  mutate(week_common=ifelse(week_common=='W000','basleine/W000',week_common)) |> 
  group_by(cell_subset,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,
              comparisons =list(c('basleine/W000','W052'),c('basleine/W000','W104'),c('W052','W104'))) %>% 
  add_significance(p.col = 'p') %>% 
  add_y_position()

pr_bl_united_w0_w52_w104_bulk_subsets %>% 
  mutate(week_common=ifelse(week_common=='W000','basleine/W000',week_common)) |> 
  mutate(week_common=factor(x=week_common,labels = c('basleine/W000','W052','W104'))) %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7,size=2)+
  facet_wrap(~cell_subset+OIT,scales = 'free_y',nrow =  3 )+
  scale_color_manual(values=c('#E69F00',"#0072B2"))+
  stat_pvalue_manual(data=stat_pr_bl_w0_w52_w104_bulk_subsets |> 
                       select(-p.adj.signif,-p.adj),
                     label = 'p.signif',hide.ns = TRUE,y.position = c(0.005,0.002))+
  #scale_y_log10(expand = expansion(mult = c(0.1, 0.1)), oob = scales::squish_infinite)+
  labs(y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/bl_PR_in_w0_w52_W104_bulk_Th2_Th2a_th1.png',
       width = 7,height = 6.5,device = 'png',dpi = 400)


#connections with cd63 
#baseline
pr_cl_bL_bulk_by_subset_not_downsampled %>% 
  inner_join(cd63_values,by='Participant_PPID') %>% 
  ggplot(aes(x=group_cd63,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(aes(color=group_cd63),outlier.shape = NA)+
  geom_point(alpha=0.5,aes(color=group_cd63),size=3)+
  facet_wrap(~cell_subset,scales = 'free_y',nrow =  2 )+
  scale_fill_manual(values=c("#105BCC","#659406","#C26A27"))+
  scale_color_manual(values=c("#105BCC","#659406","#C26A27"))+
  theme_bw()+
  stat_compare_means()+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(title='Freq of baseline PR clones in BL bulk among different BAT groups')+
  xlab('%CD63high group')+
  theme(legend.position="none")

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/cd63_subsets_bl_bulk.png',
       width = 12,height = 6,device = 'png',dpi = 400)



##PR w104 in W0,W52,w104 bulk data----
###total----
##added new variable 'week_sc_added = 'BL' that needs to delet clones that were at the baseline. To took just new ones
prCL_w104_in_bulk_BL_not_downsample<-find_sc_clones_in_bulk(data_bulk= bl01_repertoire_not_downsampled,
                                                          data_sc = cd4_seq_corrected_only_trb, week_sc = 'W104',week_bulk = 'BL01',
                                                          metadata_bulk_data = second_run_files_deepSeq_new %>% 
                                                            mutate(week_bulk_samples=case_when(week=='BL01'~'BL01',
                                                                                               week=='BL02'~'BL01',
                                                                                               week=='BL03'~'BL01',TRUE~week)))

prCL_w104_in_bulk_w0_not_downsample<-find_sc_clones_in_bulk(data_bulk=w0_not_downsampled,week_sc_added = 'BL',
                                                            data_sc = cd4_seq_corrected_only_trb, week_sc = 'W104',week_bulk = 'W000',
                                                            metadata_bulk_data =metadata_3d_w0_w52)

prCL_w104_in_bulk_w52_not_downsample<-find_sc_clones_in_bulk(data_bulk=w52_bulk_not_downsampled,week_sc_added = 'BL',
                                                            data_sc = cd4_seq_corrected_only_trb, week_sc = 'W104',week_bulk = 'W052',
                                                            metadata_bulk_data =metadata_3d_w0_w52)


prCL_w104_in_bulk_w104_not_downsample<-find_sc_clones_in_bulk(data_bulk= w104_bulk_not_downsampled,  data_sc = cd4_seq_corrected_only_trb, 
                                                            week_sc = 'W104',week_bulk = 'W104',
                                                            metadata_bulk_data = first_run_deep_seq_metadata %>% 
                                                              mutate(week_bulk_samples=case_when(week=='W104'~'W104',
                                                                                                 week=='W104A'~'W104',
                                                                                                week=='W104B'~'W104',TRUE~week)))
#to unite w0 and bl
pr104_in_united_w0_and_bl<-prCL_w104_in_bulk_BL_not_downsample %>% 
  bind_rows(prCL_w104_in_bulk_w0_not_downsample) %>% 
  mutate(week_common2='W000') |> 
  group_by(Participant_PPID,OIT,Group3,week_common2) |> 
  summarise(mean_freq_in_bulk_2=mean(mean_freq_in_bulk))

prw104_w0_w52_w104_not_downsampled<-pr104_in_united_w0_and_bl %>% 
  dplyr::rename(mean_freq_in_bulk=mean_freq_in_bulk_2,
                week_common=week_common2) |> 
  bind_rows(prCL_w104_in_bulk_w52_not_downsample) %>% 
  bind_rows(prCL_w104_in_bulk_w104_not_downsample) %>% 
  select(-week,-sum_freq_in_bulk) %>% unique() %>% 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==3) 

stat_pr_w104_w0_w52_w104_not_downsampled<-prw104_w0_w52_w104_not_downsampled %>% 
  group_by(OIT) %>% 
  wilcox_test(mean_freq_in_bulk ~ week_common, paired = T,comparisons = list(c('W052','W000'),
                                                                             c('W104','W052'),
                                                                             c('W000','W104'))) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

prw104_w0_w52_w104_not_downsampled %>% 
  mutate(week_common=factor(x=week_common,labels = c('W000','W052','W104'))) %>% 
  ggplot(aes(x=week_common,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,aes(color=OIT))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  stat_pvalue_manual(data =stat_pr_w104_w0_w52_w104_not_downsampled,label = 'p.adj',hide.ns = TRUE)+
  facet_wrap(~OIT,scales = 'free')+
  labs(y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_PRclones_Sumfreq_in_w0_w52_w104_not_downsampled.png',
       width = 7,height = 4,device = 'png',dpi = 400)

###Th1/Th2a/Th2conv----

pr_cl_w104_bulk_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= bl01_repertoire_not_downsampled %>% 
                                                            inner_join(second_run_files_deepSeq_new %>% 
                                                                         select(real_path,Participant_PPID,week) %>% unique(),by='real_path'), data_sc = cd4_seq_corrected_only_trb, 
                                                          week_sc = 'W104',week_bulk = 'BL01',week_sc_add = 'BL')
#w0
pr_cl_w104_bulk_w0_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk= w0_not_downsampled %>% 
                                                                       inner_join(metadata_3d_w0_w52 %>% 
                                                                                    select(Participant_PPID,week,real_path) %>% unique(),by='real_path'), 
                                                                     data_sc = cd4_seq_corrected_only_trb,
                                                                     week_sc = 'W104',week_bulk = 'W000',week_sc_add = 'BL')

#w52
pr_cl_w104_bulk_w052_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= w52_bulk_not_downsampled %>% 
                                                                 inner_join(metadata_3d_w0_w52 %>% select(real_path,week,Participant_PPID),
                                                                            by='real_path'),
                                                               data_sc = cd4_seq_corrected_only_trb, 
                                                               week_sc = 'W104',week_bulk = 'W052',week_sc_add = 'BL')


#w104  
pr_cl_w104_bulk_w104_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled, 
                                                                       data_sc = cd4_seq_corrected_only_trb,
                                                                       week_sc = 'W104',week_bulk = 'W104')


#united bl and w0
pr104_united_w0_bl<-pr_cl_w104_bulk_by_subset |> 
  bind_rows(pr_cl_w104_bulk_w0_not_downsampled) |> 
  mutate(week_common2='W000') |> 
  group_by(Participant_PPID,cell_subset,week_common2) |> 
  summarise(mean_freq_bl_w0=mean(sum_freq_in_bulk_by_subset)) |> 
  dplyr::rename(sum_freq_in_bulk_by_subset=mean_freq_bl_w0,
                week_common=week_common2)

#w0,w52 and w104 together
w104_pr_clones_w0_w52_w104_bulk<-pr104_united_w0_bl %>% 
  bind_rows(pr_cl_w104_bulk_w052_by_subset) %>%
  bind_rows(pr_cl_w104_bulk_w104_not_downsampled) %>% 
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>% 
  group_by(Participant_PPID,OIT,cell_subset) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==3) |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|cell_subset=='PR_Teff_Me_Th1_CTL')

#statistics
stat_w104_pr_clones_w0_w52_w104_bulk<-w104_pr_clones_w0_w52_w104_bulk %>% 
  group_by(cell_subset,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,
              comparisons =list(c('W000','W052'),c('W000','W104'),c('W052','W104'))) %>% 
  add_significance(p.col = 'p') %>% 
  add_y_position()

w104_pr_clones_w0_w52_w104_bulk %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7,size=2)+
  facet_wrap(~cell_subset+OIT,scales = 'free_y',nrow =  3 )+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  stat_pvalue_manual(data=stat_w104_pr_clones_w0_w52_w104_bulk,
                     label = 'p.signif',hide.ns = TRUE)+
  labs( y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/Subsets_w104_PR_in_w0_w52_w104_bulk.png',
       width =7 ,height = 6.5,device = 'png',dpi = 400)

###w104 PR  SU VS DS----

#w117 repertoire preparation 
w117_not_downsampled<-process_mixcr_clones(list_samples=first_run_deep_seq_metadata %>%  filter(week=='W117' | week=='W117A' | week=='W117B' ) %>% 
                                             filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID) %>% pull(real_path))

w117_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones = w117_not_downsampled)

w117_corrected_freq_not_downsampled<-make_mean_clonesFreq_forReplicates(data_full_repertoire=w117_not_downsampled,
                                                                        week_common = 'W117',week_2 = 'W117A',week_3 = 'W117B',
                                                                        metadata_bulk_data = first_run_deep_seq_metadata)
#w104 clones findings
pr_cl_w104_bulk_w104_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled, 
                                                                       data_sc = cd4_seq_corrected_only_trb,week_sc_add = 'BL',
                                                                       week_sc = 'W104',week_bulk = 'W104')
#w117 clones findings
pr_cl_w104_bulk_w117_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk= w117_corrected_freq_not_downsampled, 
                                                                       data_sc = cd4_seq_corrected_only_trb,week_sc_add = 'BL',
                                                                       week_sc = 'W104',week_bulk = 'W117')

#w104 and w117 together
w104_pr_clones_w104_and_w117_bulk<-pr_cl_w104_bulk_w104_not_downsampled %>% 
  bind_rows(pr_cl_w104_bulk_w117_not_downsampled) %>%  
  left_join(first_run_deep_seq_metadata %>% ungroup() %>% select(Participant_PPID,Group3) %>% unique(),by='Participant_PPID') %>% 
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% select(Group_dose_outcome,Participant_PPID) %>% unique(),by='Participant_PPID') %>% 
  group_by(Participant_PPID,Group3,Group_dose_outcome,cell_subset) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|cell_subset=='PR_Teff_Me_Th1_CTL')

#statistics
stat_w104_pr_clones_w104_and_w117_bulk<-w104_pr_clones_w104_and_w117_bulk %>% 
  filter(Group_dose_outcome=='Pea_0_fa'| Group_dose_outcome=='Pea_0_suc') %>% 
  group_by(cell_subset,Group_dose_outcome) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,comparisons =c('W104','W117')) %>% 
  add_significance(p.col = 'p') |> 
  add_y_position()

w104_pr_clones_w104_and_w117_bulk %>% 
  filter(cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  filter(Group_dose_outcome=='Pea_0_fa'|Group_dose_outcome=='Pea_0_suc') %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=Group_dose_outcome))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=Participant_PPID),alpha=0.3)+
  geom_point(alpha=0.7,size=2)+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  stat_pvalue_manual(data=stat_w104_pr_clones_w104_and_w117_bulk,
                      label = 'p.signif',hide.ns = TRUE,y.position = c(0.002))+
  facet_wrap(~cell_subset+Group_dose_outcome,
             scales = 'free_y',nrow =  1 )+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  labs(y='Summed freq (by Patient Id)')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th1_w104_PR_in_W104_W117_bulk_SU_vs_DS.png',
       width = 5,height = 4,device = 'png',dpi = 400)

#the same plot but facets will be weeks

#statistics
stat_w104_pr_clones_w104_and_w117_bulk_su_vs_ds<-w104_pr_clones_w104_and_w117_bulk %>% 
  filter(Group_dose_outcome=='Pea_0_fa'| Group_dose_outcome=='Pea_0_suc') %>% 
  group_by(cell_subset,week_common) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ Group_dose_outcome,
              comparisons =c('Pea_0_fa','Pea_0_suc')) %>% 
  add_significance(p.col = 'p') |> 
  add_y_position()


w104_pr_clones_w104_and_w117_bulk %>% 
  filter(cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  filter(Group_dose_outcome=='Pea_0_fa'|Group_dose_outcome=='Pea_0_suc') %>% 
  ggplot(aes(x=Group_dose_outcome,y=sum_freq_in_bulk_by_subset,color=Group_dose_outcome))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.7,size=2)+
  stat_compare_means()+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  #stat_pvalue_manual(data=stat_w104_pr_clones_w104_and_w117_bulk,
                   #  label = 'p.signif',hide.ns = TRUE,y.position = c(0.002))+
  facet_wrap(~cell_subset+week_common,
             scales = 'free_y',nrow =  1 )+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  labs(y='Summed freq (by Patient Id)')+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th1_w104_PR_in_W104_W117_bulk_SU_vs_DS_different_grouping.png',
       width = 5,height = 3.5,device = 'png',dpi = 400)

## bl pr in bl bulk, w104 in w104, w117 in w117-----

###by large subsets----

#to transform single-cell data
сd4_seq_corrected_only_trb %<>% 
  mutate(timepoint_single_cell=case_when(timepoint_single_cell=='W117B'~'W117',
                                         timepoint_single_cell=='W117B_'~'W117',
                                         TRUE~timepoint_single_cell))
#w104
pr_W104_in_bulk_w104_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled,
                                                                                 data_sc = cd4_seq_corrected_only_trb, week_sc = 'W104',
                                                                                 week_bulk = 'W104')

#week 117
pr_W117_in_bulk_w117_not_down_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w117_corrected_freq_not_downsampled,
                                                                                   data_sc = cd4_seq_corrected_only_trb, week_sc = 'W117',week_bulk = 'W117')

#bl-w104
bl_w104_large_subsets<-pr_bl_in_bulk_BL_not_down_large_subset |> 
  bind_rows(pr_W104_in_bulk_w104_not_down_large_subset) |> 
  group_by(Participant_PPID,Cell_type) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==2) |>
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> 
               select(Participant_PPID,Group_dose_outcome) |> unique(),by='Participant_PPID') |> 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  inner_join(cd4_seq_corrected |> select(OIT,Participant_PPID) |> unique(),by='Participant_PPID')

#statistics
stat_w104_pr_clones_w104_bulk<-bl_w104_large_subsets %>% 
  group_by(Cell_type,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, 
              paired = T,comparisons =list(c('BL01','W104'))) %>% 
  add_significance(p.col = 'p') |> 
  add_y_position() |> 
  select(-y.position) |> 
  unique()

bl_w104_large_subsets |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,size=2,aes(color=OIT))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  #stat_pvalue_manual(data=stat_w104_pr_clones_w104_bulk |> 
                     #  select(-'p.adj.signif',-'p.adj'),
                    # label = 'p.signif',hide.ns = TRUE,y.position = c(0.008))+
  facet_wrap(~OIT+Cell_type,scales = 'free_y',nrow =  3 )+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  labs(y='Summed freq (by Patient Id)')+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/large_subs_bl_w104_OIT_pl.png',
       width = 8,height = 6,device = 'png',dpi = 400)

#bl-w104-w117
bl_w104_w117_large_subsets<-pr_bl_in_bulk_BL_not_down_large_subset |> 
  bind_rows(pr_W104_in_bulk_w104_not_down_large_subset) |> 
  bind_rows(pr_W117_in_bulk_w117_not_down_large_subset) |> ungroup() |> 
  group_by(Participant_PPID,Cell_type) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==3) |>
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Participant_PPID,Group_dose_outcome) |> unique(),by='Participant_PPID') |> 
  filter(Group_dose_outcome=='Pea_0_fa'|Group_dose_outcome=='Pea_0_suc') |> 
  filter(Cell_type!='PR_CD154+Teff_Na')

#statistics
stat_w104_pr_clones_w104_and_w117_bulk<-bl_w104_w117_large_subsets %>% 
  group_by(Cell_type,Group_dose_outcome) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, 
              paired = T,comparisons =list(c('W104','W117'),c('BL01','W104'),c('BL01','W104'))) %>% 
  add_significance(p.col = 'p') |> 
  add_y_position() |> 
  select(-y.position) |> 
  unique()
  
bl_w104_w117_large_subsets |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=Group_dose_outcome))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,size=2,aes(color=Group_dose_outcome))+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  stat_pvalue_manual(data=stat_w104_pr_clones_w104_and_w117_bulk |> 
                       select(-'p.adj.signif',-'p.adj'),
                     label = 'p.signif',hide.ns = TRUE,y.position = c(0.008))+
  facet_wrap(~Group_dose_outcome+Cell_type,scales = 'free_y',nrow =  3 )+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  labs(y='Summed freq (by Patient Id)')+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/large_subs_bl_w104_w117.png',
       width = 8,height = 6,device = 'png',dpi = 400)

### Th1/ Th2conv/ Th2a----
#w117
#w104  
pr_cl_w117_bulk_w117_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk= w117_corrected_freq_not_downsampled, 
                                                                       data_sc = cd4_seq_corrected_only_trb,
                                                                       week_sc = 'W117',week_bulk = 'W117')

bl_w104_th1_th2a_th2conv<-pr_cl_bL_bulk_by_subset |> 
  bind_rows(pr_cl_w104_bulk_w104_not_downsampled) |> 
  filter(cell_subset=='PR_Th2conv' | cell_subset=='PR_Teff_Me_Th1_CTL' | cell_subset=='PR_Teff_Me_Th2a') |>
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
              select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>% 
  group_by(Participant_PPID,OIT,cell_subset) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==2)

stat_pr_bl1_w104_bulk_subsets<-bl_w104_th1_th2a_th2conv %>% 
  group_by(cell_subset,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,
              comparisons =list(c('BL01','W104'))) %>% 
  add_significance(p.col = 'p') %>% 
  add_y_position()

bl_w104_th1_th2a_th2conv %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.7,size=2)+
  geom_line(aes(group=Participant_PPID))+
  stat_pvalue_manual(data=stat_pr_bl1_w104_bulk_subsets,
                     label = 'p.signif',hide.ns = TRUE,y.position = c(0.0052))+
  facet_wrap(~cell_subset+OIT,scales = 'free_y',nrow =  3 )+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  labs(y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))  

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th1_th2a_th2c_bl_w104.png',
       width = 7,height = 6,device = 'png',dpi = 400)

bl_w104_w117_th1_th2a_th2conv %>% 
  ggplot(aes(x=Group_dose_outcome,y=sum_freq_in_bulk_by_subset,color=Group_dose_outcome))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7,size=2)+
  stat_compare_means()+
  facet_wrap(~cell_subset+week_common,scales = 'free_y',nrow =  3 )+
  scale_color_manual(values=c("#9966CC","#009E73"))+
  labs(y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th1_th2a_Th2conv_bl-bl_104-w104_w117-w117_PR_bulk_SU_vs_DS.png',
       width = 7,height = 6,device = 'png',dpi = 400)

#DOWNSAMPLE 500k reads-----
###BL bulk samples downasamling----
second_run_files_deepSeq_new %<>% #to add downsampled data path
  mutate(downsampled_500k_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample/',sample_name),
         downsample_file_present=downsampled_500k_path %in% 
         dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample',full.names = TRUE))

bl01_repertoire_downsampled<-process_mixcr_clones(second_run_files_deepSeq_new %>% filter(week=='BL01' | week=='BL02'|  week=='BL03') %>% 
                                                    pull(downsampled_500k_path))
bl01_repertoire_downsampled<-grouping_mixcr_clones(file_mixcr_clones = bl01_repertoire_downsampled)

pr_cl_bL_bulk<-find_sc_clones_in_bulk(data_bulk= bl01_repertoire_downsampled ,data_sc = cd4_seq_corrected_only_trb, week_sc = 'BL',week_bulk = 'BL01',
                       metadata_bulk_data = second_run_files_deepSeq_new %>% 
                         mutate(week_bulk_samples=case_when(week=='BL01'~'BL01',
                                                            week=='BL02'~'BL01',
                                                            week=='BL03'~'BL01',TRUE~week)) %>% dplyr::rename(not_downsampled_path=real_path,
                                                                                                         real_path=downsampled_500k_path))

###104 bulk samples downsampling----

w104_repertoire_downsampled<-process_mixcr_clones(first_run_deep_seq_metadata %>% filter(week=='W104' | week=='W104B'|  week=='W104A') %>% 
                                                    pull(downsampled_500k_path))
w104_repertoire_downsampled<-grouping_mixcr_clones(file_mixcr_clones = w104_repertoire_downsampled)


pr_cl_w104_bulk<-find_sc_clones_in_bulk(data_bulk= w104_repertoire_downsampled,  data_sc = cd4_seq_corrected_only_trb, 
                                        week_sc = 'W104',week_bulk = 'W104',
                                      metadata_bulk_data = first_run_deep_seq_metadata %>% 
                                        mutate(week_bulk_samples=case_when(week=='W104'~'W104',
                                                                           week=='W104A'~'W104',
                                                                           week=='W104B'~'W104',TRUE~week)) %>% 
                                       dplyr::rename(not_downsampled_path=real_path,
                                                      real_path=downsampled_500k_path))
pr_bulk_w104<-pr_cl_bL_bulk %>% 
  bind_rows(pr_cl_w104_bulk) %>% 
  select(-week,-sum_freq_in_bulk) %>% unique() %>% 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  mutate(timepoint=case_when(week_common=='BL01'~'baseline',
                             week_common=='BL02'~'baseline', 
                             week_common=='W104'~'week 104'))

pr_bulk_w104 %>% 
  ggplot(aes(x=timepoint,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(fill=timepoint))+
  geom_point(alpha=0.7)+
  geom_line(aes(group=Participant_PPID),alpha=0.5)+
  #geom_text(data=freq_pr_bulk_bl_2 %>%
             # group_by(OIT,timepoint) %>% count(),aes(x=timepoint,label=n,y=-0.01))+
  scale_fill_manual(values=c("#336666", "#FF9C99"))+
  stat_compare_means(paired = TRUE,method ='wilcox.test')+
  facet_grid(~OIT,scales = 'free_y')+
  labs(title='Summed frequency of PR clones 
       in the following repertoires on the downsampled data', y='Summed frequency',
       x='Week', fill='Week')+
  theme_bw()
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PRclones_Sumfreq_Bl_inBL_w104_inW104.png',
       width = 6,height = 4,device = 'png',dpi = 400)

##BL pr clones in the BL and w104 samples----
####total----
pr_cl_BL_bulk_w104<-find_sc_clones_in_bulk(data_bulk= w104_repertoire_downsampled, 
                                             data_sc = cd4_seq_corrected_only_trb,
                                             week_sc = 'BL',week_bulk = 'W104',
                                             metadata_bulk_data = first_run_deep_seq_metadata %>% 
                                               mutate(week_bulk_samples=case_when(week=='W104'~'W104',
                                                                                  week=='W104A'~'W104',
                                                                                  week=='W104B'~'W104',TRUE~week)) %>% 
                                               dplyr::rename(not_downsampled_path=real_path,
                                                             real_path=downsampled_500k_path))


pr_BL_bulk_BL_w104<-pr_cl_bL_bulk %>% 
  bind_rows(pr_cl_BL_bulk_w104) %>% 
  select(-week,-sum_freq_in_bulk) %>% unique() %>% 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  mutate(timepoint=case_when(week_common=='BL01'~'baseline',
                             week_common=='BL02'~'baseline', 
                             week_common=='W104'~'week 104'))

pr_BL_bulk_BL_w104 %>% 
  ggplot(aes(x=timepoint,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(fill=timepoint))+
  geom_point(alpha=0.7,size=0.3)+
  geom_line(aes(group=Participant_PPID),alpha=0.5)+
  scale_fill_manual(values=c("#336666", "#FF9C99"))+
  stat_compare_means(paired = TRUE,method ='wilcox.test')+
  facet_grid(~OIT,scales = 'free_y')+
  geom_text(data=pr_BL_bulk_BL_w104 %>% 
              group_by(OIT,timepoint) %>% 
              summarise(median_freq=median(mean_freq_in_bulk)) %>% unique() %>% mutate(median_freq=paste0('Me ',median_freq)),
            aes(x=timepoint,y=0.025,label=median_freq),size=3)+
  labs(title='Summed frequency of baseline PR clones 
       in the following repertoires on the downsampled data', y='Summed frequency',
       x='Week', fill='Week')+
  theme_bw()

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PRclones_Sumfreq_baseline_inBL_and_w104_bulk.png',
       width = 6.5,height = 4,device = 'png',dpi = 400)

##w104 clones in BL and w104 bulk----
#baseline
pr_cl_W104_bulk_baseline<-find_sc_clones_in_bulk(data_bulk= bl01_repertoire_downsampled ,data_sc = cd4_seq_corrected_only_trb,
                                      week_sc = 'W104',week_bulk = 'BL01',
                                      metadata_bulk_data = second_run_files_deepSeq_new %>% 
                                        mutate(week_bulk_samples=case_when(week=='BL01'~'BL01',week=='BL02'~'BL01',
                                                                           week=='BL03'~'BL01',TRUE~week)) %>% 
                                        dplyr::rename(not_downsampled_path=real_path,
                                                      real_path=downsampled_500k_path))
#w104
pr_cl_w104_bulk_w104<-find_sc_clones_in_bulk(data_bulk= w104_repertoire_downsampled, 
                                        data_sc = cd4_seq_corrected_only_trb,
                                        week_sc = 'W104',week_bulk = 'W104',
                                        metadata_bulk_data = first_run_deep_seq_metadata %>% 
                                          mutate(week_bulk_samples=case_when(week=='W104'~'W104', week=='W104A'~'W104',
                                                                             week=='W104B'~'W104',TRUE~week)) %>% 
                                          dplyr::rename(not_downsampled_path=real_path,
                                                        real_path=downsampled_500k_path))
pr_w104_in_bulk_baseline_w104<-pr_cl_W104_bulk_baseline %>% 
  bind_rows(pr_cl_w104_bulk_w104) %>% 
  select(-week,-sum_freq_in_bulk) %>% unique() %>% 
  group_by(Participant_PPID,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  mutate(timepoint=case_when(week_common=='BL01'~'baseline', week_common=='BL02'~'baseline',  week_common=='W104'~'week 104'))

pr_w104_in_bulk_baseline_w104 %>% 
  ggplot(aes(x=timepoint,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(fill=timepoint),alpha=0.8)+
  geom_point(alpha=0.7,size=2)+
  geom_line(aes(group=Participant_PPID),alpha=0.5)+
  scale_fill_manual(values=c("#336666", "#FF9C99"))+
  stat_compare_means(paired = TRUE,method ='wilcox.test')+
  facet_grid(~OIT,scales = 'free_y')+
  labs(title='Summed frequency of week 104 PR clones 
       in the following repertoires on the downsampled data', y='Summed frequency',
       x='Week', fill='Week')+
  theme_bw()

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PRclones_Sumfreq_w104_inBL_and_inW104.png',
       width = 6.5,height = 4,device = 'png',dpi = 400)




#cd63high and freq of PR clones correlation----
###Bl/week0 together----
united_w0_and_bl %>% 
  dplyr::rename(mean_freq_in_bulk=mean_freq_in_bulk_2,
                week_common=week_common2) |> 
  inner_join(cd63_values,by='Participant_PPID') %>%
  ggplot(aes(x=rescaled_bl,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(color=rescaled_bl))+
  geom_quasirandom(size=2,alpha=0.7,aes(color=rescaled_bl))+
  scale_color_manual(values=c("#105BCC","#659406","#C26A27"))+
  ylab('Summed frequency (by patient ID)')+
  xlab('%CD63high groups')+
  theme_classic()+
  ylim(0,0.022)+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/BAT_rescaled_prBL_blbulk.png',
       width = 6,height = 4,device = 'png',dpi = 400) 

###baseline clones by subsets----

pr_bl_united_Bl_w0 |> 
  inner_join(cd63_values,by='Participant_PPID') %>% 
  filter(sum_freq_in_bulk_by_subset<1.150444e-01) |> 
  ggplot(aes(x=rescaled_bl,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=rescaled_bl))+
  geom_quasirandom(size=2,alpha=0.7,aes(color=rescaled_bl))+
  scale_color_manual(values=c("#105BCC","#659406","#C26A27"))+
  facet_wrap(~cell_subset,scales = 'free')+
  stat_compare_means()+
  ylab('Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/BAT_rescaled_prBL_bl_bulk_subsets.png',
       width = 10,height = 6,device = 'png',dpi = 400) 
  
###scatter plot + correlation----
BAT_baseline_bl_frequncy <- ggscatter(prCL_BL_in_bulk_BL_not_downsample %>% 
                  inner_join(cd63_values,by='Participant_PPID'), x = "Wk 0", y = "sum_freq_in_bulk",
                  color = "rescaled_bl",palette = "jco", add = "reg.line",  
                conf.int = TRUE)

bat_baseline_bl_frequency<-BAT_baseline_bl_frequncy + 
  stat_cor(aes(color=rescaled_bl),method = "pearson", 
           label.x = 3)

prCL_BL_in_bulk_BL_not_downsample %>% 
  inner_join(cd63_values,by='Participant_PPID') %>% 
  select(Participant_PPID,OIT,'Wk 0',sum_freq_in_bulk,rescaled_bl) %>% 
  write_tsv('/data/vskatova/vskatova/data_milab_plots/Grouped_scatter_bulk_PR_bl_BAT.tsv')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/Grouped_BAT_results_w0_bulk_PRclones.png',
       width = 6,height = 4.5,device = 'png',dpi = 400)

prCL_w104_in_bulk_w104_not_downsample %>% 
  inner_join(cd63_values,by='Participant_PPID') %>% View()
  filter(rescaled_w104!='No inf') %>% 
  ggplot(aes(x=OIT,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(fill=OIT))+
  scale_fill_manual(values=c("#105BCC","#659406","#C26A27"))+
  geom_quasirandom(size=2,alpha=0.7)+
  stat_compare_means()+
  theme_bw()+
  facet_wrap(~rescaled_w104,scales ='free')+
  labs(title='Summed frequency of w104 PR clones in the w104 bulk')+
  ylab('Summed frequency')+
  xlab(' %CD63high groups at w104')+
  theme(legend.position="none")
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/rescaledBAT_results_w104_oit_placebo.png',
       width = 7,height = 4,device = 'png',dpi = 400)

#CTD cumulative dose correlation ----
cumulative_tol_dose<-fread('/data/vskatova/vskatova/single_cell_poised_cd4/POISED_CTD_XH.csv') #load the dataset
cumulative_tol_dose %<>% 
  mutate(week=str_extract(ID,'_.*'),
         week=str_remove(week,'_'),
         Participant_PPID=str_extract(ID,'P\\d+'))
##baseline and final CTD----
###total----
#baseline pr clones with final CTD (total frequency without subsets division)
ctd_pr_bl<-united_w0_and_bl %>% 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Group_dose_outcome,Participant_PPID) |> unique(),
             by='Participant_PPID') |> 
  filter(Group_dose_outcome=='Pea_0_fa'| Group_dose_outcome=='Pea_0_suc') |> 
  inner_join(cumulative_tol_dose %>% filter(week=='W104'),by='Participant_PPID')

ctd_pr_bl |> 
ggscatter( x = "CTD_final", y = "mean_freq_in_bulk_2",color = 'Group_dose_outcome',
           palette = c("#9966CC","#009E73"), add = "reg.line", 
           add.params = list(color = 'black', 
                             fill = "lightgray"),  conf.int = TRUE)+
  stat_cor( method = "pearson")+
  xlab('final CTD')+
  ylab('Summed frequency (by patient id)')+
  labs(title='Correlation between summed freq of Bl PR clones \
        in w0 bulk TRB PBMC and final CTD')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/final_CTD_freq_PRclones_baseline_ds_su.png',
       width = 6.5,height = 5,device = 'png',dpi = 400)

###Teff vs Tregs----
#baseline pr clones with final CTD (T effector vs Tregs)
pr_bl_w0_bl_not_down_large_sub |> 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Group_dose_outcome,Participant_PPID) |> unique(),
             by='Participant_PPID') |> 
  filter(Group_dose_outcome=='Pea_0_fa'| Group_dose_outcome=='Pea_0_suc') |> 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  inner_join(cumulative_tol_dose %>% filter(week=='W104'),by='Participant_PPID') |> 
  ggscatter( x = "CTD_final", y = "sum_freq_in_bulk_by_subset",color = 'Group_dose_outcome',
             palette = c("#9966CC","#009E73"), add = "reg.line", 
             add.params = list(color = 'black', 
                               fill = "lightgray"),  conf.int = TRUE)+
  stat_cor( method = "pearson")+
  facet_wrap(~Cell_type,scales = 'free')+
  xlab('final CTD')+
  ylab('Summed frequency (by patient id)')
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/final_CTD_freq_PRclones_baseline_ds_su_Teff_Tregs.png',
       width = 7,height = 5,device = 'png',dpi = 400)

###Th2a-Th2conv-Th1CTL----
#baseline pr clones with final CTD (Th1-Th2a-Th2conv)

pr_bl_united_w0_w52_w104_bulk_subsets |> 
  filter(week_common=='W000') |> 
  inner_join(cumulative_tol_dose %>% filter(week=='W104'),by='Participant_PPID') |> 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Group_dose_outcome,Participant_PPID) |> unique(),
             by='Participant_PPID') |> 
  filter(Group_dose_outcome=='Pea_0_suc' | Group_dose_outcome=='Pea_0_fa') |> 
  ggscatter( x = "CTD_final", y = "sum_freq_in_bulk_by_subset",color = 'Group_dose_outcome',
             palette = c("#9966CC","#009E73"), add = "reg.line", 
             add.params = list(color = 'black', 
                               fill = "lightgray"),  conf.int = TRUE)+
  stat_cor( method = "pearson")+
  facet_wrap(~cell_subset,scales = 'free')+
  xlab('final CTD')+
  ylab('Summed frequency (by patient id)')
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/final_CTD_freq_PRclones_baseline_ds_su_Th1_Th2a_Th2conv.png',
       width = 9,height = 5,device = 'png',dpi = 400) 


#baseline pr clones with final CTD (Th1_CTL/Th2conv and other combinations)

Th2conv_Th1CTL_rate<-pr_bl_united_w0_w52_w104_bulk_subsets |> 
  filter(week_common=='W000') |> 
  inner_join(cumulative_tol_dose %>% filter(week=='W104'),by='Participant_PPID') |> 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Group_dose_outcome,Participant_PPID) |> unique(),
             by='Participant_PPID') |> 
  filter(Group_dose_outcome=='Pea_0_suc' | Group_dose_outcome=='Pea_0_fa') |> 
  select(Participant_PPID,cell_subset,week_common,sum_freq_in_bulk_by_subset,Group_dose_outcome,CTD_final) |> 
  pivot_wider(names_from =cell_subset,values_from = sum_freq_in_bulk_by_subset ) |> 
  filter(PR_Th2conv!=is.na(PR_Th2conv)) |>
  mutate(rate_Th2conv_Th1=PR_Th2conv/PR_Teff_Me_Th1_CTL) |> 
  mutate(rate_Th2a_Th1=PR_Teff_Me_Th2a/PR_Teff_Me_Th1_CTL)

Th2conv_Th1CTL_rate |> 
  ggscatter( x = "CTD_final", y = "rate_Th2a_Th1",color = 'Group_dose_outcome',
           palette = c("#9966CC","#009E73"), add = "reg.line", 
           add.params = list(color = 'black', 
                             fill = "lightgray"),  conf.int = TRUE)+
  stat_cor( method = "pearson")+
  #facet_wrap(~cell_subset,scales = 'free')+
  xlab('final CTD')+
  ylab('Rate Th2a/Th1_CTL')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/final_CTD_freq_PRclones_baseline_ds_su_Th1_Th2a_Th2conv.png',
       width = 9,height = 5,device = 'png',dpi = 400)

#baseline pr clones with baseline CTD
united_w0_and_bl %>% 
  inner_join(cumulative_tol_dose %>% filter(week=='W104'),by='Participant_PPID') |> 
  ggscatter( x = "CTD", y = "mean_freq_in_bulk_2",color = 'OIT',
             palette = c("#9966CC","#009E73"), add = "reg.line", 
             add.params = list(color = 'black', 
                               fill = "lightgray"),  conf.int = TRUE)+
  stat_cor( method = "pearson")+
  xlab('baseline CTD')+
  ylab('Summed frequency (by patient id)')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/baseline_CTD_freq_PRclones_baseline.png',
       width = 6.5,height = 5,device = 'png',dpi = 400)
  
prCL_w104_in_bulk_w104_not_downsample %>% 
  inner_join(cd63_values,by='Participant_PPID') %>% View()
filter(rescaled_w104!='No inf') %>% 
  ggplot(aes(x=OIT,y=mean_freq_in_bulk))+
  geom_boxplot(outlier.shape = NA,aes(fill=OIT))+
  scale_fill_manual(values=c("#105BCC","#659406","#C26A27"))+
  geom_quasirandom(size=2,alpha=0.7)+
  stat_compare_means()+
  theme_bw()+
  facet_wrap(~rescaled_w104,scales ='free')+
  labs(title='Summed frequency of w104 PR clones in the w104 bulk')+
  ylab('Summed frequency')+
  xlab(' %CD63high groups at w104')+
  theme(legend.position="none")



#Enrichment of PR in PBMC/tissue----

final_metadata %<>% #prepare metadata file for tissue analysis
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr_V4/',sample)) %>% 
  mutate(real_path=str_replace(real_path,'.txt','clones_trb_prod.tsv'))

##tissues between each other week 0----    
#for normalization
tissues_nCLones_nReads<-read_tsv(final_metadata %>% filter(tissue!='PBMC') %>% #upload tissue repertoires week 0 
           filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID ) %>% 
            pull(real_path),id = 'real_path') |> 
  group_by(real_path) |> 
  summarise(nClones=n(),
            nReads=sum(readCount)) |> 
  inner_join(final_metadata %>% 
               select(real_path,Participant_PPID,week,tissue) %>% unique(),by='real_path')
  
week0_tissues<-process_mixcr_clones(list_samples = final_metadata %>% filter(tissue!='PBMC') %>% #upload tissue repertoires week 0 
                       filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID ) %>% 
                       filter(week=='W000') %>% pull(real_path))
week0_tissues<-grouping_mixcr_clones(week0_tissues)
week0_tissues %<>%  inner_join(final_metadata %>% 
                                select(real_path,Participant_PPID,week,tissue) %>% unique(),by='real_path')

pr_cl_bL_bulk_by_subset_not_downsampled<-find_sc_clones_in_bulk_by_subset(data_bulk =bl01_repertoire_not_downsampled %>% 
                                   inner_join(second_run_files_deepSeq_new %>% 
                                                select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                                 data_sc = cd4_seq_corrected_only_trb,week_bulk = 'BL01',week_sc = 'BL')

#for tissues i will write separate code because in the function there are no grouping variable tissue
###by large subset----

find_sc_cl_in_tissues_large_subset<-function(data_bulk,tissues_list,data_sc,week_sc,week_bulk){
  sc_clones_in_all_tissues<-lapply(tissues_list,function(tissue_chosen){
   
     data_bulk_tissue<-data_bulk |> 
      filter(tissue==tissue_chosen)
  
  sc_clones_by_tissue<-find_sc_clones_in_bulk_by_large_subset(data_bulk= data_bulk_tissue,
                                                          data_sc = data_sc, 
                                                          week_sc =week_sc,
                                                          week_bulk =week_bulk) |> 
      mutate(tissue=tissue_chosen)
  }) |> bind_rows()
  return(sc_clones_in_all_tissues)
}
####bl----
prBl_in_w0_tissues_by_large_subset<-find_sc_cl_in_tissues_large_subset(data_bulk = week0_tissues,
                                                                       data_sc = cd4_seq_corrected_only_trb,
                                                              week_bulk = 'W000',week_sc = 'BL',
                                                          tissues_list = week0_tissues |> pull(tissue) |> unique()) |> 
  inner_join(tissues_nCLones_nReads |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(sum_freq_in_bulk_by_subset_normalized=sum_freq_in_bulk_by_subset/nClones)


prBl_in_w0_tissues_by_large_subset %>% 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=Cell_type,y=sum_freq_in_bulk_by_subset_normalized))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_quasirandom(size=2,alpha=0.5,aes(color=tissue))+
  facet_wrap(~tissue)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  stat_compare_means()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))


ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/bl_pr_tissues_bl_large_subsets_normalized.png',
       width = 12,height = 5,device = 'png',dpi = 400)

####w104----

prW104_in_w104_tissues_by_large_subset<-find_sc_cl_in_tissues_large_subset(data_bulk = week104_tissues,
                                                                       data_sc = cd4_seq_corrected_only_trb,
                                                                       week_bulk = 'W104',week_sc = 'W104',
                                                                       tissues_list = week0_tissues |> pull(tissue) |> unique()) |> 
  inner_join(tissues_nCLones_nReads |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(sum_freq_in_bulk_by_subset_normalized=sum_freq_in_bulk_by_subset/nClones)

prW104_in_w104_tissues_by_large_subset %>% 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  filter(Participant_PPID %in% (cd4_seq_corrected_only_trb |> ungroup() |>
           filter(OIT=='Active') |> pull(Participant_PPID) |> unique())) |> 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_subset_normalized))+
  geom_boxplot(outlier.shape = NA,aes(color=tissue))+
  geom_quasirandom(size=2,aes(color=tissue))+
  facet_wrap(~Cell_type)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized frequency of PR clones')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_break(c(3e-06,6e-06 ), scales=0.25)

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_pr_tissues_w104_large_subsets_normalized.png',
       width = 12,height = 5,device = 'png',dpi = 400)

####bl-w104 comparison-----
w0_w104_large_subsets<-prBl_in_w0_tissues_by_large_subset |> 
  bind_rows(prW104_in_w104_tissues_by_large_subset %>% 
              filter(Cell_type!='PR_CD154+Teff_Na') |> 
              filter(Participant_PPID %in% (cd4_seq_corrected_only_trb |> ungroup() |>
                                              filter(OIT=='Active') |> pull(Participant_PPID) |> unique()))) |> 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>% 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Participant_PPID,OIT),by='Participant_PPID') |> 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  filter(OIT=='Active') %>% unique() |> 
  #group_by(Participant_PPID,Cell_type,tissue,OIT) |> 
  group_by(Participant_PPID,tissue,OIT) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==2) |> ungroup() |> 
 # ungroup() |> 
#  group_by(Participant_PPID,tissue,week_common,OIT) |> 
#  summarise(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset_normalized)) |> View()
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset_normalized))+
  geom_boxplot(outlier.shape = NA,aes(color=tissue))+
  geom_line(aes(group=Participant_PPID,color=tissue))+
  geom_point(size=2,aes(color=tissue))+
 # stat_compare_means(paired = TRUE)+
  #facet_wrap(~tissue,nrow = 1)+
  facet_wrap(~Cell_type+tissue,scales = 'free_y',nrow = 2)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized frequency of PR clones by tissue')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))
  #scale_y_break(c(3e-06,1e-05 ), scales=0.25)

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w0_w104_large_subsets_normalized.png',
       width = 14,height = 7,device = 'png',dpi = 400)

###Figure S14----

plot_S14<-plot_grid(w0_W104_tissues_normalized_no_break,w0_w104_large_subsets, 
                   labels = c('A', 'B'),nrow = 2)

ggsave(plot_S14,filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/S14_figure_GI.png',
       width =14 ,height = 15,device = 'png',dpi = 400)

###by cell subset----
find_sc_cl_in_tissues<-function(data_bulk,tissues_list,data_sc,week_sc,week_bulk){
  sc_clones_in_all_tissues<-lapply(tissues_list,function(tissue_chosen){
    data_bulk_tissue<-data_bulk |> 
      filter(tissue==tissue_chosen)
  sc_clones_by_tissue<-find_sc_clones_in_bulk_by_subset(data_bulk= data_bulk_tissue,
                                     data_sc = data_sc, 
                                     week_sc =week_sc,week_bulk =week_bulk) |> 
    mutate(tissue=tissue_chosen)
  }) |> bind_rows()
  return(sc_clones_in_all_tissues)
}

prBl_in_w0_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week0_tissues,data_sc = cd4_seq_corrected_only_trb,
                                                    week_bulk = 'W000',week_sc = 'BL',
                      tissues_list = week0_tissues |> pull(tissue) |> unique()) |> 
  inner_join(tissues_nCLones_nReads |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue)) |> 
  mutate(sum_freq_in_bulk_by_subset_normalized=sum_freq_in_bulk_by_subset/nClones)


prBl_in_w104_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week104_tissues, 
                                                      data_sc = cd4_seq_corrected_only_trb,
                                                    week_bulk = 'W104', week_sc = 'BL',
                                                    tissues_list = week104_tissues |> pull(tissue) |> unique()) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue))
  
#calculate the statistics  
stat_pr_bl_tissues_subsets_w0<-dunn_test(data=prBl_in_w0_tissues_by_subset %>% 
                                           ungroup() %>% 
                                   group_by(cell_subset),formula=sum_freq_in_bulk_by_subset_normalized ~ tissue) %>% 
  add_y_position(formula=sum_freq_in_bulk_by_subset_normalized ~ tissue) |> 
  add_significance()


#pr_bl_tissues_w0 %>% 
prBl_in_w0_tissues_by_subset %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>% 
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_subset_normalized))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_quasirandom(size=2,alpha=0.5,aes(color=tissue))+
  stat_pvalue_manual(data=stat_pr_bl_tissues_subsets_w0,hide.ns = TRUE,step.increase = 0.001,
                    label = 'p.adj.signif',y.position = c(0.00000015))+
  facet_wrap(~cell_subset,scales = 'free_y')+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/bl_pr_tissues_bl_normalized.png',
       width = 13.5,height = 6,device = 'png',dpi = 400)

###by cell subset bar plot----
cbPalette_alpha <- c( "#009999", "#CCFFFF", "#FFCC99", "#66CC99", "#993300","#FFFFCC", "#CC79A7",  "#5F31CC", "#AD3757" )
                               
pr_bl_tissues_w0_clones<-inner_join(week0_tissues, cd4_seq_corrected_only_trb %>% filter(timepoint_single_cell=='BL'),
                             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  select(cdr3aa_trb,vj_trb,freq_bulk,Participant_PPID,week,tissue,cell_subset) %>% unique() %>% 
  group_by(Participant_PPID,week,tissue) %>% 
  mutate(n_clones_in_tissue=n()) %>% ungroup() %>% 
  group_by(Participant_PPID,week,cell_subset,tissue) %>% 
  summarise(n_clones_normalized=n()/dplyr::first(n_clones_in_tissue))

write_tsv(pr_bl_tissues_w0_clones,'/data/vskatova/vskatova/data_milab_plots/barplot_pr_clones_inTissues.tsv')  

pr_bl_tissues_w0_clones %>%   
  ggplot(aes(x=tissue,y=n_clones_normalized,fill=cell_subset))+
  geom_col()+
  scale_fill_manual(values=cbPalette_alpha)+
  facet_wrap(~Participant_PPID)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title='Normalized counts of baseline PR clones found in week 0 tissue samples ')

ggsave(filename = '/data/vskatova/vskatova/data_milab_plots/barplots_PRclones_inTissues.png',
       width = 7,height = 5,device = 'png',dpi = 400)
  
###by tissue
test<-pr_bl_tissues_w0 %>% 
  group_by(tissue,Participant_PPID,week,Group_dose_outcome,OIT) %>% 
  mutate(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset))

test %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>% 
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_tissue))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(fill=tissue))+
  geom_quasirandom(size=2,alpha=0.5)+
  stat_compare_means()+
  scale_fill_manual(values = c("#009999", "#CCFFFF", "#FFCC99", "#AD3757", "#5F31CC","#FFFFCC"))+
  theme_bw()+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
  

##tissues between each other week 104----
week104_tissues<-process_mixcr_clones(list_samples = final_metadata %>% filter(tissue!='PBMC') %>% #upload tissue repertoires week 0 
                                      filter(Participant_PPID %in% cd4_seq_corrected_only_trb$Participant_PPID ) %>% 
                                      filter(week=='W104') %>% pull(real_path))

week104_tissues<-grouping_mixcr_clones(week104_tissues)
week104_tissues %<>%  inner_join(final_metadata %>% 
                                 select(real_path,Participant_PPID,week,tissue) %>% unique(),by='real_path')

#for tissues i will write separate code because in the function there are no grouping variable tissue
pr_w104_tissues_w104<-inner_join(week104_tissues, cd4_seq_corrected_only_trb %>% 
                                   filter(timepoint_single_cell=='W104'),
                             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  select(cdr3aa_trb,vj_trb,freq_bulk,Participant_PPID,week,tissue,cell_subset) %>% unique() %>% 
  group_by(Participant_PPID,week,tissue,cell_subset) %>% 
  mutate(sum_freq_in_bulk_by_subset=sum(freq_bulk)) %>% ungroup() %>% 
  complete(Participant_PPID,tissue,cell_subset,fill = list(sum_freq_in_bulk_by_subset=0)) %>% 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               filter(timepoint_single_cell=='W104') %>% 
               select(cell_subset,timepoint_single_cell,Participant_PPID,Group_dose_outcome,OIT) %>% unique(),by=c('Participant_PPID','cell_subset')) %>% 
  filter(Participant_PPID %in% week104_tissues$Participant_PPID) %>% 
  mutate(week_common='W104') %>% 
  inner_join(tissues_nCLones_nReads |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue)) |> 
  mutate(sum_freq_in_bulk_by_subset_normalized=sum_freq_in_bulk_by_subset/nClones)

###by tissue-----
pr_w104_tissues_w104_by_tissue<-pr_w104_tissues_w104 %>%  
  group_by(Participant_PPID,tissue,week_common,Group_dose_outcome,OIT) %>% 
  summarise(sum_freq_in_bulk_by_tissue_normalized=sum(sum_freq_in_bulk_by_subset_normalized)) |> 
  
  
pr_w104_tissues_w104_by_tissue %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>%
  filter(OIT=='Active') %>% #because we want to check just OIT group
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_tissue_normalized,color=tissue))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_quasirandom(size=2,alpha=0.5,aes(color=tissue))+
  stat_compare_means()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_pr_tissues_w104_by_tissue.png',
       width = 8,height = 4,device = 'png',dpi = 400)

###by cell subset----

#calculate the statistics  
stat_pr_104_tissues_subsets_w104<-dunn_test(data=pr_w104_tissues_w104 %>% 
                                              select(Participant_PPID,tissue,cell_subset,
                                                     sum_freq_in_bulk_by_subset_normalized,OIT) |> unique() |> 
                                              filter(OIT=='Active') %>% 
                                              ungroup() %>% 
                                           group_by(cell_subset),
                                           formula=sum_freq_in_bulk_by_subset_normalized ~ tissue) %>% 
  add_y_position(formula=sum_freq_in_bulk_by_subset_normalized ~ tissue)

pr_w104_tissues_w104 %>%  
  select(Participant_PPID,tissue,cell_subset,sum_freq_in_bulk_by_subset_normalized,OIT) |> unique() |> 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>%
  filter(OIT=='Active') %>% 
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_subset_normalized))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_quasirandom(size=2,alpha=0.5,aes(color=tissue))+
  stat_pvalue_manual(data=stat_pr_104_tissues_subsets_w104,label ='p.adj.signif',
                     hide.ns = TRUE,y.position = c(0.00000065,0.00000075,0.00000085))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(~cell_subset,scales = 'free_y')+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))  

ggsave(filename ='/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_pr_tissues_w104_by_cell_subset_norm.png',
       width = 13.5,height = 6,device = 'png',dpi = 400)

##tissues between each other bl and week 104 together----

pr_bl_tissues_w0 %>% 
  bind_rows(pr_w104_tissues_w104) %>% 
  group_by(Participant_PPID,tissue,cell_subset,Group_dose_outcome,OIT) %>% 
  mutate(n_tp=n()) %>% 
  filter(n_tp>1) %>% 
  group_by(Participant_PPID,tissue,Group_dose_outcome,OIT,week) %>% 
  summarise(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset)) %>% 
  filter(OIT=='Active') %>% 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue)) %>% 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>%
  filter(OIT=='Active') %>% 
  ggplot(aes(x=week,y=sum_freq_in_bulk_by_tissue,fill=week))+
  geom_boxplot(outlier.shape = NA,alpha=0.7)+
  geom_point(size=2,alpha=0.5)+
  geom_line(aes(group=Participant_PPID)) +
  facet_wrap(~tissue,scales = 'free_y')+
  stat_compare_means(paired = TRUE)+
  # stat_pvalue_manual(data=stat_pr_bl_tissues_w0, y.position = 'y.position', label = 'p.adj')+
  # geom_text(data=pr_bl_tissues_w0_by_tissue %>% group_by(tissue) %>% summarise(n=n()),
  #aes(x=tissue,label=n,y=-0.001))+
  scale_fill_manual(values = c("#009999", "#CCFFFF", "#FFCC99", "#AD3757", "#5F31CC","#FFFFCC"))+
  theme_bw()+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
  
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w0_w104_by_tissue.png',
       width = 7.5,height = 4.5,device = 'png',dpi = 400)

##frequency pr clones in PBMC vs tissues baseline----
pr_bl_pbmc_tissues_w0_by_tissue<-prBl_in_w0_tissues_by_subset |> 
  group_by(Participant_PPID,tissue,week_common) |> 
  summarise(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset)) |> 
  inner_join(tissues_nCLones_nReads |> 
               mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                                       tissue=='mid_esophagus'~'mid_esoph',
                                       tissue=='dist_esophagus'~'dist_esoph',
                                       TRUE~tissue)) |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(sum_freq_in_bulk_by_tissue_normalized=sum_freq_in_bulk_by_tissue/nClones) |> 
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> select(Participant_PPID,OIT) |> unique(),
             by='Participant_PPID')

pr_bl_pbmc_tissues_w0_by_tissue_pbmc<-pr_bl_pbmc_tissues_w0_by_tissue |> 
  bind_rows(united_w0_and_bl |> 
              ungroup() |> 
              mutate(tissue='pbmc') |> 
              select(tissue,Participant_PPID,OIT,week_common2,mean_freq_in_bulk_normalized) |> 
              dplyr::rename(sum_freq_in_bulk_by_tissue_normalized=mean_freq_in_bulk_normalized,
                            week_common=week_common2)) |> 
  unique() |> 
  group_by(week_common,Participant_PPID,OIT) |> 
  mutate(n_tissues=n()) |> 
  filter(n_tissues>1)

stat_pr_bl_tissues_vs_pbmc_w0<-dunn_test(data=pr_bl_pbmc_tissues_w0_by_tissue %>% ungroup(),
                                 formula=sum_freq_in_bulk_by_subset_normalized ~ tissue) %>% 
  add_y_position(formula=sum_freq_in_bulk_by_subset_normalized ~ tissue)

pr_bl_pbmc_tissues_w0_by_tissue_pbmc %>%
  mutate(tissue=factor(x=tissue,levels=c('pbmc','prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>%
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_tissue_normalized))+
  geom_boxplot(outlier.shape = NA,aes(color=tissue))+
  geom_quasirandom(size=2,aes(color=tissue))+
  scale_color_manual(values=c("#A6761D","#1B9E77","#D95F02","#7570B3", "#E7298A", "#66A61E","#E6AB02"))+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_break(c(2e-06,6e-06 ), scales=0.25)

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/bl_pr_w0_tissues_bl01_pbmc_normalized.png',
       width = 7,height = 5,device = 'png',dpi = 400)

##frequency pr clones in PBMC vs tissues w104----


nReads_nClones_meaned_1st_run<-nClones_nReads_run1 |>   #nClones_nReads_run1 somwwhere from the upper part
  inner_join(first_run_deep_seq_metadata,by='real_path') |> 
  group_by(Participant_PPID,Group) |> 
  summarise(nClones=mean(nclones))

pr_w104_pbmc_tissues_w104_by_tissue<-pr_w104_tissues_w104 %>% 
  group_by(Participant_PPID,tissue,week_common) %>% 
  summarise(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset)) %>% 
#  group_by(Participant_PPID,week_common) %>% mutate(n=n()) %>% 
 # filter(n==6) %>% 
  inner_join(tissues_nCLones_nReads |> 
               mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                                       tissue=='mid_esophagus'~'mid_esoph',
                                       tissue=='dist_esophagus'~'dist_esoph',
                                       TRUE~tissue)) |> 
               dplyr::rename(week_common=week),by=c('Participant_PPID','week_common','tissue')) |> 
  mutate(sum_freq_in_bulk_by_tissue_normalized=sum_freq_in_bulk_by_tissue/nClones) |> 
  bind_rows(pr_cl_w104_bulk_w104_not_downsampled %>% mutate(tissue='pbmc') %>% 
              mutate(week_common='W104') |> 
              group_by(Participant_PPID,week_common,tissue) |> 
              summarise(sum_freq_in_bulk_by_tissue=sum(sum_freq_in_bulk_by_subset)) |> 
              inner_join(nReads_nClones_meaned_1st_run,by='Participant_PPID') |> 
              mutate(sum_freq_in_bulk_by_tissue_normalized=sum_freq_in_bulk_by_tissue/nClones)) |> 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>% 
  filter(OIT=='Active') |> 
  group_by(Participant_PPID,OIT) |> 
  mutate(n_tissues=n()) |> 
  filter(n_tissues>1)

pr_w104_pbmc_tissues_w104_by_tissue %>%
  mutate(tissue=factor(x=tissue,levels=c('pbmc','prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>%
  ggplot(aes(x=tissue,y=sum_freq_in_bulk_by_tissue_normalized))+
  geom_boxplot(outlier.shape = NA,alpha=0.7,aes(color=tissue))+
  geom_quasirandom(size=2,aes(color=tissue))+
  #scale_fill_brewer(palette = "Dark2")+
  scale_color_manual(values=c("#A6761D","#1B9E77","#D95F02","#7570B3", "#E7298A", "#66A61E"))+
  labs(y='Normalized summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))  

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_pr_w104_tissues_104_pbmc_normalized.png',
       width = 6,height = 4,device = 'png',dpi = 400)

###W0-w104 comp----

w0_W104_tissues_normalized_no_break<-pr_bl_pbmc_tissues_w0_by_tissue_pbmc |> 
  bind_rows(pr_w104_pbmc_tissues_w104_by_tissue) |> 
  filter(tissue!='pbmc') |> 
  mutate(tissue=factor(x=tissue,
                      levels=c('prox_esoph','mid_esoph','dist_esoph','stomach','duodenum'))) %>% 
  filter(OIT=='Active') %>% unique() |> 
  group_by(Participant_PPID,tissue,OIT) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==2) |> 
  ungroup() |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_tissue_normalized))+
  geom_boxplot(outlier.shape = NA,aes(color=tissue))+
  geom_line(aes(group=Participant_PPID,color=tissue))+
  geom_point(size=2,aes(color=tissue))+
  facet_wrap(~tissue,nrow = 1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(y='Normalized frequency of PR clones by tissue')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")
  #scale_y_break(c(5e-06,1e-05 ), scales=0.25)

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w0_W104_tissues_normalized_no_break.png',
       width = 8,height = 5,device = 'png',dpi = 400)

##frequency of bl pr clones in w0 and w104 tissues-----

prBl_in_w0_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week0_tissues,data_sc = cd4_seq_corrected_only_trb,
                                                    week_bulk = 'W000',week_sc = 'BL',
                                                    tissues_list = week0_tissues |> pull(tissue) |> unique()) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue))

prBl_in_w104_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week104_tissues, 
                                                      data_sc = cd4_seq_corrected_only_trb,
                                                      week_bulk = 'W104', week_sc = 'BL',
                                                      tissues_list = week104_tissues |> pull(tissue) |> unique()) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue))

prBL_w0_w104_tissues<-prBl_in_w0_tissues_by_subset |> 
  bind_rows(prBl_in_w104_tissues_by_subset) |> 
  group_by(Participant_PPID,cell_subset,tissue) |> 
  mutate(n_timepoints=n_distinct(week_common)) |> 
  filter(n_timepoints>1) |> 
  inner_join(cd4_seq_corrected_only_trb |> 
               ungroup() |> select(Participant_PPID,OIT) |> unique(),by='Participant_PPID')

stat_pr_bl_w0_w52_w104_bulk_by_subset<-prBL_w0_w104_tissues %>% 
  filter(OIT!='Placebo') |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|
           cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  ungroup() |> 
  group_by(OIT,tissue,cell_subset) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,comparisons = list(c('W104','W000'))) %>% 
  #adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

prBL_w0_w104_tissues %>% 
  filter(OIT!='Placebo') |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|
           cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  mutate(week_common=factor(x=week_common,labels = c('W000','W104'))) %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,aes(color=OIT))+
  scale_color_manual(values=c("#E69F00"))+
  facet_wrap(~tissue+cell_subset,scales = 'free',ncol = 3)+
  labs(y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/bl_pr_w0_w104_in_tissues.png',
       width = 8,height = 8.5,device = 'png',dpi = 400)

##frequency of w104 pr clones in w0 and w104----
pr_w104_in_w0_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week0_tissues,
                                                       data_sc = cd4_seq_corrected_only_trb,
                                                    week_bulk = 'W000',week_sc = 'W104',
                                                    tissues_list = week0_tissues |> pull(tissue) |> unique()) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue))

pr_w104_in_w104_tissues_by_subset<-find_sc_cl_in_tissues(data_bulk = week104_tissues, 
                                                      data_sc = cd4_seq_corrected_only_trb,
                                                      week_bulk = 'W104', week_sc = 'W104',
                                                      tissues_list = week104_tissues |> pull(tissue) |> unique()) |> 
  mutate(tissue=case_when(tissue=='prox_esophagus'~'prox_esoph',
                          tissue=='mid_esophagus'~'mid_esoph',
                          tissue=='dist_esophagus'~'dist_esoph',
                          TRUE~tissue))

pr_w104_w0_w104_tissues<-pr_w104_in_w0_tissues_by_subset |> 
  bind_rows(pr_w104_in_w104_tissues_by_subset) |> 
  group_by(Participant_PPID,cell_subset,tissue) |> 
  mutate(n_timepoints=n_distinct(week_common)) |> 
  filter(n_timepoints>1) |> 
  inner_join(cd4_seq_corrected_only_trb |> 
               ungroup() |> select(Participant_PPID,OIT) |> unique(),by='Participant_PPID')

stat_pr_w104_w0_w52_w104_bulk_by_subset<-pr_w104_w0_w104_tissues %>% 
  filter(OIT!='Placebo') |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|
           cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  ungroup() |> 
  group_by(OIT,tissue,cell_subset) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,comparisons = list(c('W104','W000'))) %>% 
  #adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

pr_w104_w0_w104_tissues %>% 
  filter(OIT!='Placebo') |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|
           cell_subset=='PR_Teff_Me_Th1_CTL') |> 
  mutate(week_common=factor(x=week_common,labels = c('W000','W104'))) %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,aes(color=OIT))+
  scale_color_manual(values=c("#E69F00"))+
  facet_wrap(~tissue+cell_subset,scales = 'free',ncol = 3)+
  labs(y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/w104_pr_w0_w104_in_tissues.png',
       width = 8,height = 8.5,device = 'png',dpi = 400)

##frequency BL clones in BL, w104 clones in W104----
###total----
pr_bl_pbmc_tissues_w0_by_tissue |> 
  bind_rows(pr_w104_pbmc_tissues_w104_by_tissue) |> 
  group_by(Participant_PPID,tissue,OIT) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp>1) |> 
  filter(OIT=='Active') |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_tissue_normalized))+
  geom_boxplot(aes(color=tissue))+
  geom_point(aes(color=tissue))+
  geom_line(aes(group=Participant_PPID,color=tissue))+
  facet_grid(~tissue,scales='free')+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  stat_compare_means(paired = TRUE)+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/corr_PR_cells_w0_w104_in_tissues.png',
       width = 9,height = 4.5,device = 'png',dpi = 400)

#CLUSTERS----

#downsample 700k
report_downsample_700k_run2<-read_tsv(file=dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample_700kReads/',
                                              pattern = '.tsv',full.names = TRUE), id = 'filename')
bl_bulk_downsampled_700k_run2_sampleList<-report_downsample_700k_run2 %>% 
  mutate(week=str_extract(sample,'W\\d+')) %>% 
  select(sample,sumWeightAfter,nElementsAfter,week) %>% 
  dplyr::rename(N_sequenced_reads=sumWeightAfter,
                N_clones=nElementsAfter) %>% 
  filter(N_sequenced_reads==700000) %>% 
  mutate(week=ifelse(is.na(week),str_extract(sample,'BL\\d+'),week)) %>% 
  mutate(sample=str_replace(sample,'.clns','.trb.tsv')) %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample_700kReads/',sample)) %>% 
  filter(week=='BL01'| week=='BL02') %>% 
  pull(real_path)

bl_700k_downsampled<-process_mixcr_clones(list_samples = bl_bulk_downsampled_700k_run2_sampleList)

#Sys.setenv(PATH="/software/anaconda3/bin") 
#setwd("/data/vskatova/vskatova/ALICE-master")
#source("ALICE.R")

find_alice_clones<-function(processed_mixcr_clonsets,week_bulk,list_samples){
  
  lapply(list_samples,function(sample){
    clonsets<-processed_mixcr_clonsets %>% 
      filter(real_path==sample) %>% 
      dplyr::rename(CDR3.amino.acid.sequence=cdr3aa_trb,
                    CDR3.nucleotide.sequence=cdr3nt_trb,
                    Read.count=count,
                    bestVGene=v_trb,
                    bestJGene=j_trb) %>% 
      filter(Read.count>1) %>% 
      as.data.table() %>% 
      list()
    table_2<-ALICE_pipeline_OLGA (DTlist=clonsets, cores=20 )
    table_2_2<-table_2[[1]]
    table_2_2 %<>% 
      mutate(week_bulk)
    
    path<-str_replace(sample,'poised_trb_deepSeq/mixcr_run1/downsample_700kReads/','poised_trb_deepSeq/analysis/alice_hits_700k/')
    print(path)
    txt_pbmc_path<-str_replace(path,'.trb.tsv','_alice_hits.txt')
    fwrite(table_2_2,txt_pbmc_path)
  })
  }
##baseline----
done_files_bl<-dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/alice_hits_700k/') %>% 
  as.data.frame() %>% 
  dplyr::rename(sample='.') %>% 
  mutate(sample=str_replace(sample,'_alice_hits.txt','.trb.tsv')) %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample_700kReads/',sample)) %>% 
  pull(real_path)

bl_files<-bl_700k_downsampled %>%  pull(real_path) %>% unique() #to make list of all files

alice_baseline_clones<-find_alice_clones(processed_mixcr_clonsets = bl_700k_downsampled,week_bulk = 'BL', #run the command
                                         list_samples = bl_files[!bl_files %in% done_files_bl])

##W104----
report_downsample_700k_run1<-read_tsv(file=dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/downsample_700kReads/',
                                               pattern = 'report.tsv',full.names = TRUE), id = 'filename')
report_downsample_700k_run1_sampleList<-report_downsample_700k_run1 %>% 
  mutate(week=str_extract(sample,'W\\d+')) %>% 
  select(sample,sumWeightAfter,nElementsAfter,week) %>% 
  dplyr::rename(N_sequenced_reads=sumWeightAfter,
                N_clones=nElementsAfter) %>% 
  filter(N_sequenced_reads==700000) %>% 
  mutate(sample=str_replace(sample,'.clns','.trb.tsv')) %>% 
  mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run1/downsample_700kReads/',sample)) %>% 
  filter(week=='W104') %>% 
  pull(real_path)

w104_700k_downsample<-process_mixcr_clones(list_samples = report_downsample_700k_run1_sampleList)
w104_files<-w104_700k_downsample %>%  pull(real_path) %>% unique() #to make list of all files
alice_w104_clones<-find_alice_clones(processed_mixcr_clonsets = w104_700k_downsample,#run the command
                                         week_bulk = 'W104', list_samples = w104_files)


##baseline alice clusters processing----
BL_Alice_clusters<-read_csv(file=dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/alice_hits_700k',
                                               pattern = 'BL',full.names = TRUE), id = 'filename') %>% 
  dplyr::rename(cdr3aa_trb=CDR3.amino.acid.sequence,cdr3nt_trb=CDR3.nucleotide.sequence) %>% 
  select(real_path,filename,cdr3aa_trb,cdr3nt_trb,vj_trb,week_bulk) %>% 
  mutate(Participant_PPID=str_extract(filename,'P\\d+')) %>% 
  group_by(real_path,filename,cdr3aa_trb,vj_trb,week_bulk,Participant_PPID) %>% 
  summarise(n_cdr3nt_alice=n_distinct(cdr3nt_trb)) %>% 
  mutate(cdr3aa_vj_trb=paste0(cdr3aa_trb,'_',vj_trb))

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

ALICE_graph_Baseline<- make.sequence.graph( BL_Alice_clusters %>%
                                              filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% 
                                              pull(cdr3aa_vj_trb), max_errs = 1)

ALICE_graph_Baseline<-set.vertex.attribute(ALICE_graph_Baseline, 'participant_label', 
                                            V(ALICE_graph_Baseline),BL_Alice_clusters %>% 
                                             filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% pull(Participant_PPID))

ALICE_graph_Baseline<-set.vertex.attribute(ALICE_graph_Baseline, 'n_CDR3nt', 
                                            V(ALICE_graph_Baseline), BL_Alice_clusters %>% 
                                             filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% pull(n_cdr3nt_alice))

write_graph(ALICE_graph_Baseline,
            '/data/vskatova/vskatova/data_milab_plots/baseline_clones_alice.gml',format = 'gml')

##logoseq----
# I need to set cluster id
adj_matrix <- stringdistmatrix(BL_Alice_clusters %>%
                                 filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% 
                                 pull(cdr3aa_vj_trb), 
                               BL_Alice_clusters %>%
                                 filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% 
                                 pull(cdr3aa_vj_trb), method = "hamming")
adj_matrix <- 1*(adj_matrix <= 1)
rownames(adj_matrix) <- BL_Alice_clusters %>%
  filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% 
  pull(cdr3aa_vj_trb)
colnames(adj_matrix) <- BL_Alice_clusters %>%
  filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7') %>% 
  pull(cdr3aa_vj_trb)
diag(adj_matrix) <- 0
g1 <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "undirected",
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
components <- components(g1)
BL_Alice_clusters_mini<- BL_Alice_clusters %>%
  filter(vj_trb=='TRBV5-1_TRBJ2-1'| vj_trb=='TRBV7-8_TRBJ2-7')

BL_Alice_clusters_mini$cluster_id<- components$membership

BL_Alice_clusters_mini_id<-BL_Alice_clusters_mini %>% 
  group_by(cluster_id) %>% 
  mutate(n_clones_in_cluster=n_distinct(cdr3aa_vj_trb)) %>% 
  filter(n_clones_in_cluster>=10) %>% 
  inner_join(vdjdb_database,by=c('cdr3aa_trb','vj_trb'))

write_tsv(BL_Alice_clusters_mini_id,'/data/vskatova/vskatova/data_milab_plots/logoseq.tsv')

seq_list <- list(cluster3 = BL_Alice_clusters_mini_id %>% filter(cluster_id=='3') %>% pull(cdr3aa_trb), 
                 cluster4 = BL_Alice_clusters_mini_id %>% filter(cluster_id=='4') %>% pull(cdr3aa_trb))

ggseqlogo(seq_list, ncol=2,method = 'prob')+
  theme_classic()

ggsave(filename = '/data/vskatova/vskatova/data_milab_plots/logoplot.png',
       width = 7,height = 2,device = 'png',dpi = 400)


bl_700k_downsampled<-grouping_mixcr_clones(bl_700k_downsampled)

second_run_files_deepSeq_new %<>% 
  mutate(downsample_700k_path=paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/mixcr_run2/downsample_700kReads/',sample_name))

PRclones_in_bulk_bl<-right_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
                                  filter(timepoint_single_cell=='BL'), bl_700k_downsampled %>% 
                                  inner_join(second_run_files_deepSeq_new %>% 
                                               dplyr::rename(not_downsampled_path=real_path,
                                                             real_path=downsample_700k_path) %>% 
                                               select(real_path,Participant_PPID) %>% unique(), by='real_path'),
                                  by=c('Participant_PPID','cdr3aa_trb','vj_trb'))

#baseline PR in bulk trb baseline
recentage_of_PRclones_in_bulk_bl<-PRclones_in_bulk_bl %>% 
  mutate(found_in_bulk=ifelse(is.na(Count_single_cell),'no','yes')) %>% 
  select(Group_dose_outcome,Participant_PPID,cdr3aa_trb,vj_trb,OIT,real_path,found_in_bulk) %>% unique() %>% 
  filter(found_in_bulk=='yes') %>% 
  group_by(Group_dose_outcome,OIT,Participant_PPID,real_path) %>% 
  summarise(nPR_found_inBulk=n()) %>% 
  inner_join(bl_700k_downsampled %>% group_by(real_path) %>% summarise(n_clones_inBulkSample=n()),by='real_path') %>% 
  mutate(Normalized_n_prCL=nPR_found_inBulk/dplyr::first(n_clones_inBulkSample))

#baseline PR in alice baseline
Bl_prClones_in_aliceClusters<-inner_join(BL_Alice_clusters,cd4_seq_corrected_only_trb %>%  ungroup() %>% 
             filter(timepoint_single_cell=='BL'),by=c('cdr3aa_trb','vj_trb','Participant_PPID'))

##w104 clusters processing----
W104_Alice_clusters<-read_csv(file=dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/alice_hits_700k',
                                     pattern = 'W104',full.names = TRUE), id = 'filename') %>% 
  dplyr::rename(cdr3aa_trb=CDR3.amino.acid.sequence,cdr3nt_trb=CDR3.nucleotide.sequence) %>% 
  select(real_path,filename,cdr3aa_trb,cdr3nt_trb,vj_trb) %>% 
  mutate(Participant_PPID=str_extract(filename,'P\\d+')) %>% 
  group_by(real_path,filename,cdr3aa_trb,vj_trb,Participant_PPID) %>% 
  summarise(n_cdr3nt_alice=n_distinct(cdr3nt_trb))

w104_prClones_in_aliceClusters<-inner_join(W104_Alice_clusters,cd4_seq_corrected_only_trb %>%  ungroup() %>% 
                                        filter(timepoint_single_cell=='W104'),by=c('cdr3aa_trb','vj_trb','Participant_PPID'))
#VDJdb----
#upload vdjdb database
vdjdb_database<-fread('/data/vskatova/vskatova/SearchTable-2023-08-03 13_20_36.002.tsv') %>% 
  mutate_at(c('V','J'), ~str_replace_all(string=., pattern='\\*\\d+',replacement = '')) %>% 
  mutate(vj_trb=paste0(V,'_',J)) %>% 
  dplyr::rename(cdr3aa_trb=CDR3)

cd4_vdjdb<-vdjdb_database %>% 
  filter(`MHC class`=='MHCII')

inner_join(cd4_seq_corrected_only_trb %>% separate(col=vj_trb,into = c('V','J'),sep = '_'),
           cd4_vdjdb,by=c('cdr3aa_trb','V','J')) %>% View()

week0_alice_clusters_annotation<-inner_join(BL_Alice_clusters,vdjdb_database,by=c('cdr3aa_trb','vj_trb')) %>% 
  filter(Score>=2) %>% 
  dplyr::rename(MHC_class='MHC class') %>% 
  select(-complex.id,-Reference,-Meta,-Epitope,-`Epitope gene`,-Method,-CDR3fix) %>% unique() %>% 
  group_by(cdr3aa_trb,vj_trb,MHC_class) %>% 
    summarise(n_pat=n_distinct(Participant_PPID),
              epitopes=paste(unique(`Epitope species`),collapse=','))

week104_alice_clusters_annotation<-inner_join(W104_Alice_clusters,vdjdb_database,by=c('cdr3aa_trb','vj_trb')) %>% 
  filter(Score>=2) %>% 
  dplyr::rename(MHC_class='MHC class') %>% 
  select(-complex.id,-Reference,-Meta,-Epitope,-`Epitope gene`,-Method,-CDR3fix) %>% unique() %>% 
  group_by(cdr3aa_trb,vj_trb,MHC_class) %>% 
  summarise(n_pat=n_distinct(Participant_PPID),
            score=dplyr::first(Score),
            epitopes=paste(unique(`Epitope species`),collapse=',')) %>%
  mutate(epitopes=case_when(epitopes=='InfluenzaA,synthetic,EBV'~'EBV,InfluenzaA',
                            epitopes=='InfluenzaA,synthetic'~'InfluenzaA',TRUE~epitopes))

w104_annotated_clusters<-right_join(week104_alice_clusters_annotation,W104_Alice_clusters,
           by=c('cdr3aa_trb','vj_trb')) %>% 
  group_by(cdr3aa_trb,vj_trb,MHC_class,epitopes) %>% 
  summarise(n_pat=n_distinct(Participant_PPID),
            epitopes=paste(unique(epitopes),collapse=',')) %>% 
  mutate(week='W104')

bl_annotated_clusters<-right_join(week0_alice_clusters_annotation,BL_Alice_clusters,
                                    by=c('cdr3aa_trb','vj_trb')) %>% 
  group_by(cdr3aa_trb,vj_trb,MHC_class,epitopes) %>% 
  summarise(n_pat=n_distinct(Participant_PPID),
            epitopes=paste(unique(epitopes),collapse=',')) %>% 
  mutate(week='W000')

bl_annotated_clusters %>% 
  bind_rows(w104_annotated_clusters) %>% 
  ggplot(aes(x=week))+
  geom_bar(aes(fill = epitopes))+
  scale_fill_manual(values=c("#009999","#FFCC99","#5F31CC",'#CCFFCC','#CC0996','#999999'))+
  theme_bw()
 # scale_y_log10()
  
  


#PCA of trajectories----
##only bulk---




##sc+bulk----
#w0
bl01_repertoire_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones= bl01_repertoire_not_downsampled)

bl01_bulk<-bl01_repertoire_not_downsampled %>% #to find PR clones that are in bulk w104 from any sc timepoint
  inner_join(second_run_files_deepSeq_new %>% select(Participant_PPID,real_path,week) %>% unique(),by='real_path') %>% 
  inner_join(cd4_seq_corrected_only_trb,by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
  unique()

#w36
w36_bulk_not_downsampled<-process_mixcr_clones(list_samples = second_run_files_deepSeq_new %>% filter(week=='W36'| week=='W38'| week=='W34') %>% 
                                                  pull(real_path))
w36_bulk_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones= w36_bulk_not_downsampled)

w36_bulk<-w36_bulk_not_downsampled %>% 
  inner_join(second_run_files_deepSeq_new %>% select(Participant_PPID,real_path,week) %>% unique(),by='real_path') %>% 
  inner_join(cd4_seq_corrected_only_trb,by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
  unique()

#w52
w52_bulk_not_downsampled<-process_mixcr_clones(list_samples= metadata_3d_w0_w52 %>% 
                                                 inner_join(final_metadata,by=c('Participant_PPID','specimen_label')) %>% 
                                                 filter(week=='W052') %>% 
                                            pull(real_path.x))

w52_bulk_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones= w52_bulk_not_downsampled)

w52<-w52_bulk_not_downsampled %>% 
  inner_join(metadata_3d_w0_w52,by='real_path') %>% 
  inner_join(final_metadata %>% 
               select(week,Participant_PPID,specimen_label,Group),
             by=c('Participant_PPID','specimen_label')) %>% 
  inner_join(cd4_seq_corrected_only_trb,by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
  unique() %>% 
  filter(week=='W052')

#w104
w104_bulk<-w104_bulk_corrected_freq_not_downsampled %>% #to find PR clones that are in bulk w104
  inner_join(cd4_seq_corrected_only_trb,by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
  unique()
 
bl_w52_w104_whole_rep_patIds<-bl01_repertoire_not_downsampled %>% #to define patients with 3 time points
  inner_join(second_run_files_deepSeq_new,by='real_path') %>% 
  mutate(week='W000') %>% 
  bind_rows(w52_bulk_not_downsampled %>% 
              inner_join(metadata_3d_w0_w52,by='real_path') %>% 
              inner_join(final_metadata %>% 
                           select(week,Participant_PPID,specimen_label,Group),
                         by=c('Participant_PPID','specimen_label')) %>% filter(week=='W052')) %>% 
  bind_rows(w104_bulk_corrected_freq_not_downsampled) %>% 
  group_by(Participant_PPID) %>% 
  summarise(n_week=n_distinct(week),
            week=paste(unique(week),collapse = ',')) %>% 
  filter(n_week==3) %>% 
  pull(Participant_PPID)

bl_w52_w104_freq<-bl01_bulk %>% 
  bind_rows(w52) %>% 
  bind_rows(w104_bulk) %>% 
  filter(Participant_PPID %in% bl_w52_w104_whole_rep_patIds ) %>% 
  pivot_wider(id_cols = c('cdr3aa_trb','vj_trb','Participant_PPID'),
              names_from = 'week',values_from = 'freq_bulk',values_fill = 0) %>% 
  mutate(cdr3aa_vj_trb=paste0(cdr3aa_trb,'_',vj_trb)) %>% 
  select(cdr3aa_vj_trb,BL01,W052,W104,Participant_PPID) %>% 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID')
 

##graph frequency basic----
bl_w52_w104_freq %>% 
  pivot_longer(cols = c('BL01','W052','W104'),
               names_to ='week',values_to = 'freq_bulk') %>% 
  mutate(week=factor(x=week,levels=c('BL01','W052','W104'))) %>% 
    ggplot(aes(x=week,group= cdr3aa_vj_trb, y=freq_bulk)) +
    geom_line(alpha=0.3)+
    geom_point(alpha=0.5,size=2)+
    ylim(0,0.001)

ggsave(filename = '/data/vskatova/vskatova/data_milab_plots/bl_w36_w104.png',
       width = 5,height = 3,device = 'png',dpi = 400)

##PCA code----
#lets make separately for Placebo and OIT
#OIT first

fviz_nbclust(only_log10_freq_OIT, kmeans, method = "wss") #elbow plot
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/elbowplot_w0_w52_w104.png',
       width = 5,height = 3,device = 'png',dpi = 400)

k2<-kmeans(x=only_log10_freq_OIT,centers = 2,nstart = 25)
k3<-kmeans(x=only_log10_freq_OIT,centers = 3,nstart = 25)
k4<-kmeans(x=only_log10_freq_OIT,centers = 4,nstart = 25)
k5<-kmeans(x=only_log10_freq_OIT,centers = 5,nstart = 25)
k6<-kmeans(x=only_log10_freq_OIT,centers = 6,nstart = 25)

pl2<-fviz_cluster(k2, geom = "point", data = only_log10_freq_OIT)+
  ggtitle("k = 2")
pl3<-fviz_cluster(k3, geom = "point", data = only_log10_freq_OIT)+
  ggtitle("k = 3")
pl4<-fviz_cluster(k4, geom = "point", data = only_log10_freq_OIT)+
  ggtitle("k = 4")
pl5<-fviz_cluster(k5, geom = "point", data = only_log10_freq_OIT)+
  ggtitle("k = 5")
pl6<-fviz_cluster(k6, geom = "point", data = only_log10_freq_OIT)+
  ggtitle("k = 6")
grid.arrange(pl2, pl3, pl4, pl5,pl6, nrow = 2)

freq_Placebo_log10<-bl_w52_w104_freq %>% #to calculate log10 frequency
  mutate_at(c('BL01','W052','W104'),~ ifelse(.x == 0, min(.x[.x>0])/2, .x )) %>% 
  mutate(BL01=log10(BL01),W052=log10(W052),W104=log10(W104)) %>% 
  filter(OIT=='Placebo')


only_log10_freq_OIT<-freq_OIT_log10 %>%  select(BL01,W052,W104) # to select just values for the test
only_log10_freq_Placebo<-freq_Placebo_log10 %>%  select(BL01,W052,W104) # to select just values for the test

#autoplot(kmeans(only_log10_freq_Placebo,3), data = only_log10_freq_Placebo) # for UMAP
 # ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/3clusters_pca_w0_w52_w104_placebo.png',
    #  3 width = 4,height = 3,device = 'png',dpi = 400)

freq_OIT_log10 %>%  #graph with line and connected dots
  mutate(cluster_id=kmeans(only_log10_freq,3 )[['cluster']]) %>%
  mutate(cdr3aa_vj_trb_ppid=paste0(Participant_PPID,'_',cdr3aa_vj_trb)) %>% 
  pivot_longer(cols = c('BL01','W052','W104'),
                 names_to ='week',values_to = 'freq_bulk') %>% 
  mutate(week=factor(x=week,levels=c('BL01','W052','W104'))) %>% 
  ggplot(aes(x=week,group= cdr3aa_vj_trb_ppid, y=freq_bulk)) +
  geom_line(alpha=0.3)+
  geom_point(alpha=0.3,size=2)+
  facet_wrap(~cluster_id)+
  labs(y='log10_freq_in_bulk_repertoire', x='week')
  
ggsave(filename = '/data/vskatova/vskatova/single_cell_poised_cd4/analysis/BL_w36_w104_freq_6clusters.png',
         width = 7,height = 5,device = 'png',dpi = 400)

log10_bar_error<- #freq_Placebo_log10 %>% 
  freq_OIT_log10 %>%
   mutate(cluster_id=kmeans(only_log10_freq_OIT,centers=3,iter.max = 100,nstart = 100,)[['cluster']]) %>%
   mutate(cdr3aa_vj_trb_ppid=paste0(Participant_PPID,'_',cdr3aa_vj_trb)) %>% 
   pivot_longer(cols = c('BL01','W052','W104'), names_to ='week',values_to = 'freq_bulk')

write_tsv(log10_bar_error,'/data/vskatova/vskatova/data_milab_plots/log10_frequency_in_clusters_bl_w52_w104.tsv')

log10_bar_error<-fread('/data/vskatova/vskatova/data_milab_plots/log10_frequency_in_clusters_bl_w52_w104.tsv')

data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)) }
  data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
  }  

df2 <- data_summary(log10_bar_error, varname="freq_bulk", 
                    groupnames=c("cluster_id", "week"))
df2$week=as.factor(df2$week)
df2$cluster_id=as.factor(df2$cluster_id)

df2 %>% 
  mutate(week=factor(x=week,levels=c('BL01','W052','W104'))) %>%
  ggplot(aes(x=week, y=freq_bulk, group=cluster_id, color=cluster_id)) + 
  geom_line(alpha=0.8) +
  geom_point(size=3)+
  facet_grid(~cluster_id)+
  geom_errorbar(aes(ymin=freq_bulk-sd, ymax=freq_bulk+sd),position=position_dodge(0.05))+
  scale_color_manual(values=c("#198020","#5F31CC","#C26A27","#068A94","#A324B2","#659406"))+
  labs(title="Average PCA cluster frequency over time", y="log10 avg clone frequency", x = "week")+
  theme_bw()

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/pca_3clust_placebo_sd.png',
       width = 6,height = 4,device = 'png',dpi = 400)

write_tsv(log10_bar_error,'/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/6clusters_pca.tsv')

##digging into interesting clusters----

subsets_frequncy_subsets<-log10_bar_error %>% 
  separate(col=cdr3aa_vj_trb,into = c('cdr3aa_trb','v','j'),sep = '_') %>% 
  mutate(vj_trb=paste0(v,'_',j)) %>% 
  select(-v,-j,-cdr3aa_vj_trb_ppid) %>% 
  inner_join(cd4_seq_corrected_only_trb,
             by=c('cdr3aa_trb','vj_trb','Participant_PPID','OIT')) %>% 
  select(cdr3aa_trb,vj_trb,Participant_PPID,cluster_id,OIT,cell_subset,Cell_type) %>% 
  unique() %>% 
  group_by(cluster_id) %>% 
  dplyr::mutate(n_clones_in_subset=n())

#barplot with proportions of cell subsets
subsets_frequncy_subsets_barplot<-subsets_frequncy_subsets %>% 
  group_by(cluster_id,cell_subset) %>% 
  dplyr::summarise(Normalized_n_clonest=n()/dplyr::first(n_clones_in_subset)) %>% 
  ungroup() 
  
subsets_frequncy_subsets_barplot %>%   
  mutate(cluster_id=as.factor(cluster_id)) %>% 
  ggplot(aes(x=cluster_id,y=Normalized_n_clonest,fill=cell_subset))+
  geom_col()+
  scale_fill_manual(values=cbPalette_alpha)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title='Proportion of cell subsets in different trajectory clusters ')+
  ylab('normalized n clones')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/cell_subsets_in_3_PCA_clusters_OIT.png',
       width = 6,height = 4,device = 'png',dpi = 400)

subsets_frequncy_subsets_by_ppid<- subsets_frequncy_subsets %>% 
  group_by(cluster_id,Participant_PPID,cell_subset) %>% 
  dplyr::summarise(Normalized_n_clonest_by_patient=n()/dplyr::first(n_clones_in_subset))

subsets_frequncy_subsets_by_ppid %>% 
  ungroup() %>% 
  mutate(cluster_id=as.factor(cluster_id)) %>%
  ggplot(aes(x=cluster_id,y=Normalized_n_clonest_by_patient,color=cluster_id))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=c('#009999','#993300','#66CC99'))+
  geom_point(size=2.5,alpha=1)+
  theme_bw()+
  facet_wrap(~cell_subset,nrow = 2,scale='free_y')+
  stat_compare_means()+
  labs(title='Proportion of cell subsets in different trajectory clusters',
       y='Normalized n of clones by patient',
       x='PCA trajectory cluster id')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  theme(legend.position="none")

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/PCA_clusters_n_clones_byPpid.png',
       width = 12,height = 6,device = 'png',dpi = 400)

#SU vs DS pca----
w117_bulk_not_downsampled<-process_mixcr_clones(list_samples = first_run_deep_seq_metadata %>% filter(week=='W117'| week=='W117A'| week=='W117B') %>% 
                                                 pull(real_path))

w117_bulk_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones= w117_bulk_not_downsampled)

w117_bulk_not_downsampled<-make_mean_clonesFreq_forReplicates(data_full_repertoire=w117_bulk_not_downsampled,
                                                                             week_common = 'W117',week_2 = 'W117A',week_3 = 'W117B',
                                                                             metadata_bulk_data = first_run_deep_seq_metadata)
w117_bulk<-w117_bulk_not_downsampled %>% 
  inner_join(cd4_seq_corrected_only_trb,by=c('cdr3aa_trb','vj_trb','Participant_PPID'),multiple = "all") %>% 
  select(-cell_subset,-Cell_type,-n_cdr3nt_variants_sc,-Count_single_cell,-timepoint_single_cell) %>% 
  unique()

patients_4_tp<-bl01_repertoire_not_downsampled %>% #to define patients with 4 time points (bl-w52-w104-w117)
  inner_join(second_run_files_deepSeq_new,by='real_path') %>% 
  mutate(week='W000') %>% 
  bind_rows(w52_bulk_not_downsampled %>% 
              inner_join(metadata_3d_w0_w52,by='real_path') %>% 
              inner_join(final_metadata %>% 
                           select(week,Participant_PPID,specimen_label,Group),
                         by=c('Participant_PPID','specimen_label')) %>% filter(week=='W052')) %>% 
  bind_rows(w104_bulk_corrected_freq_not_downsampled) %>% 
  bind_rows(w117_bulk_not_downsampled) %>% 
  group_by(Participant_PPID) %>% 
  summarise(n_week=n_distinct(week),
            week=paste(unique(week),collapse = ',')) %>% 
  filter(n_week==4) %>% #to filter ones with bl w52 w104 w117
  pull(Participant_PPID)

bl_w52_w104_w117freq<-bl01_bulk %>% 
  bind_rows(w52) %>% 
  bind_rows(w104_bulk) %>% 
  bind_rows(w117_bulk) %>% 
  filter(Participant_PPID %in% patients_4_tp ) %>% 
  pivot_wider(id_cols = c('cdr3aa_trb','vj_trb','Participant_PPID'),
              names_from = 'week',values_from = 'freq_bulk',values_fill = 0) %>% 
  mutate(cdr3aa_vj_trb=paste0(cdr3aa_trb,'_',vj_trb)) %>% 
  select(cdr3aa_vj_trb,BL01,W052,W104,W117,Participant_PPID) %>% 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               select(Participant_PPID,OIT,Group_dose_outcome) %>% unique(),by='Participant_PPID')

##pca code----
#lets make separately for peanut 300 and peanut 0
#lets find out what number of clusters is optimal

freq_peanut300_log10<-bl_w52_w104_w117freq %>% #to calculate log10 frequency
  mutate_at(c('BL01','W052','W104','W117'),~ ifelse(.x == 0, min(.x[.x>0])/2, .x )) %>% 
  mutate(BL01=log10(BL01),W052=log10(W052),W104=log10(W104),W117=log10(W117)) %>% 
  filter(Group_dose_outcome=='Pea_300_fa'| Group_dose_outcome=='Pea_300_suc' )

fviz_nbclust(only_freq_peanut300_log10, kmeans, method = "wss")

only_freq_peanut300_log10<-freq_peanut300_log10 %>%  select(BL01,W052,W104,W117) # to select just values for the test

#only_log10_freq_Placebo<-freq_Placebo_log10 %>%  select(BL01,W052,W104) # to select just values for the test



autoplot(kmeans(only_freq_peanut300_log10,3), data = only_freq_peanut300_log10) # for UMAP

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/3clusters_pca_w0_w52_w104_placebo.png',
       width = 4,height = 3,device = 'png',dpi = 400)

freq_peanut300_log10 %>%  #graph with line and connected dots
  mutate(cluster_id=kmeans(only_freq_peanut300_log10,6 )[['cluster']]) %>%
  mutate(cdr3aa_vj_trb_ppid=paste0(Participant_PPID,'_',cdr3aa_vj_trb)) %>% 
  pivot_longer(cols = c('BL01','W052','W104','W117'),
               names_to ='week',values_to = 'freq_bulk') %>% 
  mutate(week=factor(x=week,levels=c('BL01','W052','W104','W117'))) %>% 
  ggplot(aes(x=week,group= cdr3aa_vj_trb_ppid, y=freq_bulk)) +
  geom_line(alpha=0.3)+
  geom_point(alpha=0.3,size=2)+
  facet_wrap(~cluster_id)+
  labs(y='log10_freq_in_bulk_repertoire', x='week')

ggsave(filename = '/data/vskatova/vskatova/single_cell_poised_cd4/analysis/BL_w36_w104_freq_6clusters.png',
       width = 7,height = 5,device = 'png',dpi = 400)

freq_peanut300_log10_bar_error<- freq_peanut300_log10 %>% #freq_Placebo_log10 %>% 
  mutate(cluster_id=kmeans(only_freq_peanut300_log10,centers=6,iter.max = 100,nstart = 100,)[['cluster']]) %>%
  mutate(cdr3aa_vj_trb_ppid=paste0(Participant_PPID,'_',cdr3aa_vj_trb)) %>% 
  pivot_longer(cols = c('BL01','W052','W104','W117'), names_to ='week',values_to = 'freq_bulk')

write_tsv(freq_peanut300_log10_bar_error,'/data/vskatova/vskatova/data_milab_plots/log10_frequency_in_clusters_bl_w52_w104_w117_pean300.tsv')

log10_bar_error<-fread('/data/vskatova/vskatova/data_milab_plots/log10_frequency_in_clusters_bl_w52_w104_w117_pean300.tsv')

data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)) }
  data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}  

df2 <- data_summary(log10_bar_error, varname="freq_bulk", 
                    groupnames=c("cluster_id", "week"))
df2$week=as.factor(df2$week)
df2$cluster_id=as.factor(df2$cluster_id)

df2 %>% 
  mutate(week=factor(x=week,levels=c('BL01','W052','W104','W117'))) %>%
  ggplot(aes(x=week, y=freq_bulk, group=cluster_id, color=cluster_id)) + 
  geom_line(alpha=0.8) +
  geom_point(size=3)+
  facet_grid(~cluster_id)+
  geom_errorbar(aes(ymin=freq_bulk-sd, ymax=freq_bulk+sd),position=position_dodge(0.05))+
  scale_color_manual(values=c("#198020","#5F31CC","#C26A27","#068A94","#A324B2","#659406"))+
  labs(title="Average PCA cluster frequency over time", y="log10 avg clone frequency", x = "week")+
  theme_bw()

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/pca_3clust_placebo_sd.png',
       width = 6,height = 4,device = 'png',dpi = 400)

subsets_frequncy_subsets<-log10_bar_error %>% 
  separate(col=cdr3aa_vj_trb,into = c('cdr3aa_trb','v','j'),sep = '_') %>% 
  mutate(vj_trb=paste0(v,'_',j)) %>% 
  select(-v,-j,-cdr3aa_vj_trb_ppid) %>% 
  inner_join(cd4_seq_corrected_only_trb,
             by=c('cdr3aa_trb','vj_trb','Participant_PPID','OIT')) %>% 
  select(cdr3aa_trb,vj_trb,Participant_PPID,cluster_id,OIT,cell_subset,Cell_type) %>% 
  unique() %>% 
  group_by(cluster_id) %>% 
  dplyr::mutate(n_clones_in_subset=n())

#barplot with proportions of cell subsets
subsets_frequncy_subsets_barplot<-subsets_frequncy_subsets %>% 
  group_by(cluster_id,cell_subset) %>% 
  dplyr::summarise(Normalized_n_clonest=n()/dplyr::first(n_clones_in_subset)) %>% 
  ungroup() 

subsets_frequncy_subsets_barplot %>%   
  mutate(cluster_id=as.factor(cluster_id)) %>% 
  ggplot(aes(x=cluster_id,y=Normalized_n_clonest,fill=cell_subset))+
  geom_col()+
  scale_fill_manual(values=cbPalette_alpha)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title='Proportion of cell subsets in different trajectory clusters ')+
  ylab('normalized n clones')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/cell_subsets_in_3_PCA_clusters_OIT.png',
       width = 6,height = 4,device = 'png',dpi = 400)

subsets_frequncy_subsets_by_ppid<- subsets_frequncy_subsets %>% 
  group_by(cluster_id,Participant_PPID,cell_subset) %>% 
  dplyr::summarise(Normalized_n_clonest_by_patient=n()/dplyr::first(n_clones_in_subset))

subsets_frequncy_subsets_by_ppid %>% 
  ungroup() %>% 
  mutate(cluster_id=as.factor(cluster_id)) %>%
  ggplot(aes(x=cluster_id,y=Normalized_n_clonest_by_patient,color=cluster_id))+
  geom_boxplot(outlier.shape = NA)+
  #scale_color_manual(values=c('#009999','#993300','#66CC99'))+
  geom_point(size=2.5,alpha=1)+
  theme_bw()+
  facet_wrap(~cell_subset,nrow = 2,scale='free_y')+
  stat_compare_means()+
  labs(title='Proportion of cell subsets in different trajectory clusters',
       y='Normalized n of clones by patient',
       x='PCA trajectory cluster id')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  theme(legend.position="none")

#Minervina and Pogoreliy PCA trajectories----

#artem's way of pr in 3 tp----
pr_clones_3tp<-cd4_seq_corrected_only_trb |> 
  group_by(Participant_PPID,Group_dose_outcome,Cell_type,cell_subset,OIT,cdr3aa_trb,vj_trb) |> 
  summarise(n_timepoints=n_distinct(timepoint_single_cell),
            timepoints=paste(timepoint_single_cell,collapse=',')) |> 
  filter(n_timepoints>=3) |> 
  filter(str_detect(timepoints,'BL')) |> #at least one of this is BL
  select(-n_timepoints,-timepoints) |> 
  unique()

#bl
bl_3_timepoints_pr<-left_join(pr_clones_3tp |> ungroup(),bl01_repertoire_not_downsampled |> 
             inner_join(second_run_files_deepSeq_new %>% 
                          select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
           by=c('cdr3aa_trb','vj_trb','Participant_PPID')) |> 
  select(-Cell_type,-cell_subset) |> 
  unique() |> 
  mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT) |> 
  summarise(sum_freq=sum(freq_bulk)) |> 
  mutate(week_common='BL01') |> 
  filter(Participant_PPID %in% (nClones_ncReads_realpath |> 
                                  filter(week=='BL01'|week=='BL02'|week=='BL03') |> 
           pull(Participant_PPID)))

#w0
w0_3_timepoints_pr<-left_join(pr_clones_3tp |> ungroup(),w0_not_downsampled,
          by=c('cdr3aa_trb','vj_trb','Participant_PPID')) |> 
  select(-Cell_type,-cell_subset) |> 
  unique() |> 
  mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT) |> 
  summarise(sum_freq=sum(freq_bulk)) |> 
  mutate(week_common='W000') |> 
  filter(Participant_PPID %in% w0_not_downsampled$Participant_PPID)

#uniting w0 and bl
w0_bl_united<-bl_3_timepoints_pr |> 
  dplyr::rename(week=week_common) |> 
  left_join(nClones_ncReads_realpath |> 
              mutate(week=case_when((Participant_PPID=='P034' & week=='BL02') ~'BL01',
                     (Participant_PPID=='P074' & week=='BL02') ~'BL01',
                     TRUE~week)),by=c('Participant_PPID','week')) |> 
  bind_rows(w0_3_timepoints_pr |> 
              dplyr::rename(week=week_common) |> 
              left_join(nClones_ncReads_realpath,by=c('Participant_PPID','week'))) |> 
  mutate(week_common2='W000') |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT,week_common2) |> 
  mutate(mean_sum_freq=mean(sum_freq))
  dplyr::rename(week_common=week_common2,
                sum_freq=mean_sum_freq)
  

#w52

w52_3_timepoints_pr<-left_join(pr_clones_3tp |> ungroup(), w52_bulk_not_downsampled |> 
                                inner_join(metadata_3d_w0_w52 %>% 
                                             select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                              by=c('cdr3aa_trb','vj_trb','Participant_PPID')) |> 
  select(-Cell_type,-cell_subset) |> 
  unique() |> 
  mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT) |> 
  summarise(sum_freq=sum(freq_bulk)) |> 
  mutate(week_common='W052')

#w104
w104_3_timepoints_pr<-left_join(pr_clones_3tp |> ungroup(), w104_bulk_corrected_freq_not_downsampled,
                               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) |> 
  select(-Cell_type,-cell_subset) |> 
  unique() |> 
  mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT) |> 
  summarise(sum_freq=sum(freq_bulk)) |> 
  mutate(week_common='W104')

pr_bl_w0_w52_w104_not_downsampled_3timepoints<-w0_bl_united %>% 
  bind_rows(w52_3_timepoints_pr) %>% 
  bind_rows(w104_3_timepoints_pr)

stat_pr_0_bl_w0_w52_w104_not_downsampled<-pr_bl_w0_w52_w104_not_downsampled_3timepoints %>% ungroup |> 
  wilcox_test(sum_freq ~ week_common, paired = F,comparisons = list(c('W052','W000'),
                                                                    c('W104','W052'),
                                                                    c('W000','W104'))) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

pr_bl_w0_w52_w104_not_downsampled_3timepoints %>% 
  ggplot(aes(x=week_common,y=sum_freq))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7)+
  stat_pvalue_manual(data=stat_pr_0_bl_w0_w52_w104_not_downsampled,
                     tip.length = 0.001,
                     y.position = 'y.position',
                     label = 'p.adj.signif',hide.ns = TRUE)+
  labs(title='Summed frequency of PR clones that were found in at least 3 sc timepoints',
       y='Summed frequency (by Patient Id)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/3tp_pr_bl_w0_w52_w104.png',
       width =8,height = 4,device = 'png',dpi = 400) 

#SU vs DS
w117_3_timepoints_pr<-left_join(pr_clones_3tp |> ungroup(), w117_bulk_not_downsampled,
                                by=c('cdr3aa_trb','vj_trb','Participant_PPID')) |> 
  select(-Cell_type,-cell_subset) |> 
  unique() |> 
  mutate(freq_bulk=ifelse(is.na(freq_bulk),0,freq_bulk)) |> 
  group_by(Participant_PPID,Group_dose_outcome,OIT) |> 
  summarise(sum_freq=sum(freq_bulk)) |> 
  mutate(week_common='W117')

pr_w104_w117_not_downsampled_3timepoints<-w117_3_timepoints_pr |> 
  bind_rows(w104_3_timepoints_pr)

stat_pr_pr_w104_w117_not_downsampled<-pr_w104_w117_not_downsampled_3timepoints %>% ungroup |> 
  group_by(Group_dose_outcome) |> 
  wilcox_test(sum_freq ~ week_common, paired = TRUE,comparisons = c('W104','W117')) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_significance() %>% 
  add_y_position()

pr_w104_w117_not_downsampled_3timepoints %>% 
  filter(week_common=='W117') |> 
  ggplot(aes(x=Group_dose_outcome,y=sum_freq))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(alpha=0.7)+
  labs(title='Summed frequency of PR clones that were found in at least 3 sc timepoints',
       y='Summed frequency (by Patient Id)')+
  theme_classic()+
  stat_compare_means()+
  #facet_wrap(~Group_dose_outcome,scales = 'free')+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

###subsets in different clusters----
decreasing_their_adundance<-log10_bar_error %>% 
  filter(cluster_id=='4')

increasing_their_adundance<-log10_bar_error %>% 
  filter(cluster_id=='2')

increasing_their_adundance_w104<-log10_bar_error %>% 
  filter(cluster_id=='6')

inc_decr_their_adundance_w104<-log10_bar_error %>% 
  filter(cluster_id=='5')

n_clones_per_pat_and_subset<-log10_bar_error %>% 
  inner_join(cd4_seq_corrected_only_trb,
             by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  select(-week,-freq_bulk,-Count_single_cell,
         -n_cdr3nt_variants_sc) %>% 
  unique() %>%
  mutate(cell_subset=case_when(cell_subset=='PR_TNa_1'~'PR_Naive',
                               cell_subset=='PR_TNa_2'~'PR_Naive',
                               cell_subset=='PR_TNa_3'~'PR_Naive',
                               cell_subset=='PR_TNa_4'~'PR_Naive',
                               cell_subset=='PR_Teff_Me_act1'~'PR_Teff_memory_activ',
                               cell_subset=='PR_Teff_Me_act2'~'PR_Teff_memory_activ',
                               cell_subset=='PR_Teff_Me_act3'~'PR_Teff_memory_activ',
                               cell_subset=='PR_Treg_act1'~'PR_Treg_activ',
                               cell_subset=='PR_Treg_act2'~'PR_Treg_activ',
                               TRUE~cell_subset)) %>%
  mutate(large_subset=case_when(cell_subset=='PR_Naive'~'PR_Naive',
                                cell_subset=='PR_Teff_memory_activ'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Th2a'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Th22Th17'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Th1_CTL'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Th1Th17'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Tfh'~'PR_Teff',
                                cell_subset=='PR_Teff_Me_Tfh13'~'PR_Teff',
                                cell_subset=='PR_Treg_activ'~'PR_Treg_activ',
                                )) %>% 
  group_by(large_subset,Participant_PPID,cluster_id) %>% 
  dplyr::summarise(n=n()) %>% 
  ungroup() %>% 
  inner_join(final_metadata %>% 
               select('Week 117.ChallengeResult',Participant_PPID) %>% unique(),by='Participant_PPID') %>% 
  dplyr::rename(W117_FC_result='Week 117.ChallengeResult') %>% 
  complete(cluster_id,Participant_PPID,large_subset,fill = list(n=0)) %>% 
 mutate(W117_FC_result= case_when(Participant_PPID=='P012'~'FAILURE',
                                  Participant_PPID=='P025'~'SUCCESS',
                                  Participant_PPID=='P034'~'SUCCESS',
                                  Participant_PPID=='P036'~'SUCCESS',
                                   TRUE~W117_FC_result))

n_clones_per_pat_and_subset %>% 
    filter(W117_FC_result!='FAILURE') %>% 
    mutate(cluster_id=factor(x=cluster_id,
                             levels=c('1','2','3','4','5','6'))) %>% 
    ggplot(aes(x=large_subset,y=n,fill=large_subset))+
    scale_fill_manual(values=c("#105BCC","#659406","#C26A27"))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0),
               size=2)+
    theme_bw()+
    stat_compare_means()+
    facet_wrap(~cluster_id)+
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
    labs(title='Number of clones assigned to different trajectory 
       clusters per patient in different subsets',
         y='Number of clones',
         x='Subset',
         fill='Subset')
  
  ggsave(filename = '/data/vskatova/vskatova/single_cell_poised_cd4/analysis/PCA_clusters_n_clones_per_subset.png',
         width =9,height = 5,device = 'png',dpi = 400) 


###accurate subsets----
  
n_clones_per_pat_and_accurate_subset<-log10_bar_error %>% 
    inner_join(cd4_seq_corrected_only_trb,
               by=c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
    select(-week,-freq_bulk,-Count_single_cell,
           -n_cdr3nt_variants_sc) %>% 
    unique() %>%
    mutate(cell_subset=case_when(cell_subset=='PR_TNa_1'~'PR_Naive',
                                 cell_subset=='PR_TNa_2'~'PR_Naive',
                                 cell_subset=='PR_TNa_3'~'PR_Naive',
                                 cell_subset=='PR_TNa_4'~'PR_Naive',
                                 cell_subset=='PR_Teff_Me_act1'~'PR_Teff_memory_activ',
                                 cell_subset=='PR_Teff_Me_act2'~'PR_Teff_memory_activ',
                                 cell_subset=='PR_Teff_Me_act3'~'PR_Teff_memory_activ',
                                 cell_subset=='PR_Treg_act1'~'PR_Treg_activ',
                                 cell_subset=='PR_Treg_act2'~'PR_Treg_activ',
                                 TRUE~cell_subset)) %>% 
    group_by(cell_subset,Participant_PPID,cluster_id) %>% 
    dplyr::summarise(n=n()) %>% 
    ungroup() %>% 
    inner_join(final_metadata %>% 
                 select('Week 117.ChallengeResult',Participant_PPID) %>% unique(),by='Participant_PPID') %>% 
    dplyr::rename(W117_FC_result='Week 117.ChallengeResult') %>% 
    complete(cluster_id,Participant_PPID,cell_subset,fill = list(n=0)) %>% 
    mutate(W117_FC_result= case_when(Participant_PPID=='P012'~'FAILURE',
                                     Participant_PPID=='P025'~'SUCCESS',
                                     Participant_PPID=='P034'~'SUCCESS',
                                     Participant_PPID=='P036'~'SUCCESS',
                                     TRUE~W117_FC_result))  

n_clones_per_pat_and_accurate_subset %>% 
  filter(W117_FC_result!='FAILURE') %>% 
  mutate(cluster_id=factor(x=cluster_id,
                           levels=c('1','2','3','4','5','6'))) %>% 
  ggplot(aes(x=cell_subset,y=n,fill=cell_subset))+
  scale_fill_manual(values=mi_dark)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0),
             size=2)+
  theme_bw()+
  stat_compare_means()+
  facet_wrap(~cluster_id)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  labs(title='Number of clones assigned to different trajectory 
       clusters per patient in different subsets',
       y='Number of clones',
       x='Subset',
       fill='Subset')

ggsave(filename = '/data/vskatova/vskatova/single_cell_poised_cd4/analysis/PCA_clusters_n_clones_per_subset.png',
       width =9,height = 5,device = 'png',dpi = 400)   
  






  


##between SU vs DS----
#to get w117 repertoire
w117_repertoires<-read_tsv(deep_seq_run1_metadata %>% #to group repertoires by cdr3aa_vdj
                                  filter(Participant_PPID %in% th2_clones_bulk_freq$Participant_PPID) %>% #to filter what patients have sc data
                                  filter(week=='W117'|week=='W117B' | week=='W117A') %>% 
                                  pull(real_path),id = 'real_path') %>% 
  inner_join(deep_seq_run1_metadata,by='real_path') %>% 
  dplyr::rename(v_trb=allVHitsWithScore,
                j_trb=allJHitsWithScore,
                cdr3aa_trb=aaSeqCDR3,
                cdr3nt_trb=nSeqCDR3,
                count=readCount,
                freq=readFraction) %>% 
  select(-targetQualities,-allDHitsWithScore,-allVAlignments,-allDAlignments,-allJAlignments,-allCAlignments,
         -minQualCDR3,-refPoints,-targetSequences) %>% 
  mutate(v_trb=str_remove(v_trb,'\\*\\d+.*')) %>% 
  mutate(j_trb=str_remove(j_trb,'\\*\\d+.*')) %>% 
  group_by(real_path,cdr3aa_trb,v_trb,j_trb,sample_name,Participant_PPID,
           week,Group,Group2,Group3,age,week104_FC_result,week117_FC_result,
           week143_FC_result,week156_FC_result) %>% 
  summarise(count_bulk=sum(count),
            freq_bulk=sum(freq),
            n_cdr3nt_bulk=n_distinct(cdr3nt_trb)) %>% 
  ungroup() %>% 
  mutate(vj_trb=paste0(v_trb,'_',j_trb))

week117a_w117b_mean<-w117_repertoires %>% #to make mean clone freq because of the replicates
  select(cdr3aa_trb,vj_trb,Participant_PPID,week,freq_bulk) %>% 
  mutate(week=case_when(Participant_PPID=='P114'~'W117', #to change because there is no w117a/117b sample
                        TRUE~week)) %>% 
  filter(week=='W117A' | week=='W117B') %>% #to take only weeks with replicates
  dplyr::rename(freq_bulk_w117=freq_bulk) %>% 
  pivot_wider(names_from =week,values_from = freq_bulk_w117,
              id_cols = c('cdr3aa_trb','vj_trb','Participant_PPID')) %>% 
  mutate(W117A=ifelse(is.na(W117A),0,W117A),
         W117B=ifelse(is.na(W117B),0,W117B)) %>% 
  mutate(freq_bulk_w117=rowMeans(select(.,c('W117A','W117B')))) %>% 
  select(-W117A,-W117B)

w117_repertoire_clones<-w117_repertoires %>%  #these are clones with w117 without replicates
  select(cdr3aa_trb,vj_trb,Participant_PPID,week,freq_bulk) %>% 
  mutate(week=case_when(Participant_PPID=='P114'~'W117',
                        TRUE~week)) %>% 
  filter(week=='W117') %>% 
  dplyr::rename(freq_bulk_w117=freq_bulk) %>% 
  bind_rows(week117a_w117b_mean)

Th2_clones_tracking_w117<-left_join(th2_clones_bulk_freq %>% 
                                 select(cdr3aa_trb,vj_trb,Participant_PPID,OIT,timepoint_single_cell,
                                        cell_subset,week117_FC_result,freq_bulk,week) %>% 
                                 dplyr::rename(freq_bulk_baseline=freq_bulk),
                               w117_repertoire_clones,
                               by=c('cdr3aa_trb','vj_trb','Participant_PPID'))  %>%  
  filter(Participant_PPID %in% (first_run_files_deepSeq %>% 
           filter(week=='W117'|week=='W117A'|week=='W117B') %>% pull(Participant_PPID))) %>% 
  filter (Participant_PPID %in% (second_run_files_deepSeq %>% 
                                   filter(week=='BL01' | week=='BL02' | week=='BL03' ) %>% 
                                   pull(Participant_PPID))) %>% 
  mutate(freq_bulk_w117=ifelse(is.na(freq_bulk_w117),0,freq_bulk_w117),
         freq_bulk_baseline=ifelse(is.na(freq_bulk_baseline),0,freq_bulk_baseline)) %>% 
  filter(!(freq_bulk_baseline==0 & freq_bulk_w117==0)) %>% 
  mutate(week.x=ifelse(is.na(week.x),'BL','BL'),
         week.y=ifelse(is.na(week.y),'W117','W117')) %>% 
  select(-cell_subset,-timepoint_single_cell)

th2_clones_w117<-Th2_clones_tracking_w117 %>% 
  select('week.y' ,'freq_bulk_w117','cdr3aa_trb',
         'vj_trb','Participant_PPID','OIT','week117_FC_result') %>% 
  dplyr::rename(freq_bulk=freq_bulk_w117,
                week=week.y)

th2_clones_bl_w117<-Th2_clones_tracking_w117 %>% 
  dplyr::rename(freq_bulk=freq_bulk_baseline,
                week=week.x) %>% 
  select(-week.y,-freq_bulk_w117) %>% 
  bind_rows(th2_clones_w117)

th2_clones_bl_w117 %>% #graph of th2a clones trajectories by clones
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               select(Participant_PPID,Group_dose_outcome) %>% unique(),
             by='Participant_PPID') %>% 
  filter(Group_dose_outcome=='Pea_0_fa' | Group_dose_outcome=='Pea_0_suc') %>% 
  mutate(cdr3aa_trb_vj=paste0(cdr3aa_trb,'_',vj_trb)) %>% 
  ggplot(aes(x=week,y=freq_bulk))+
  geom_line(aes(group=cdr3aa_trb_vj),alpha=0.5)+
  geom_point(alpha=0.5)+
  facet_wrap(~Group_dose_outcome,scales = 'free_y')+
  stat_compare_means(paired = TRUE,method ='wilcox.test')

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th2_changing_freq_by_clones_w117.png',
       width = 6,height = 4,device = 'png',dpi = 400)

th2_clones_bl_w117 %>% #graph of th2a clones trajectories by patient 
  inner_join(cd4_seq_corrected_only_trb %>% ungroup() %>% 
               select(Participant_PPID,Group_dose_outcome) %>% unique(),
             by='Participant_PPID') %>% 
  filter(Group_dose_outcome=='Pea_0_fa' | Group_dose_outcome=='Pea_0_suc') %>% 
  group_by(Participant_PPID,Group_dose_outcome,week) %>% 
  summarise(bulk_freq_sum=sum(freq_bulk)) %>% 
  ggplot(aes(x=week,y=bulk_freq_sum))+
  geom_line(aes(group=Participant_PPID),alpha=0.5)+
  geom_point(alpha=0.5)+
  facet_wrap(~Group_dose_outcome,scales = 'free_y')+
  stat_compare_means(paired = TRUE,method ='wilcox.test' )

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/th2_changing_freq_by_pat_id_w117_su_ds.png',
       width = 6,height = 4,device = 'png',dpi = 400)




#shared top ----

find_n_shared_clones_between_2samples<-function(data_tissues,data_pbmc,list_of_ppid){

    pair_comparison<-function(ppid){
    
    table_pbmc<-data_pbmc %>% 
      filter(Participant_PPID==ppid) %>% 
      filter(count>1)
    
    table_tissues<-data_tissues %>% 
      filter(Participant_PPID==ppid) %>% 
      filter(count>1)
    
    contamination<-inner_join(table_pbmc,table_tissues,
                              by=c('cdr3nt_trb','v_trb','j_trb','Participant_PPID'))#to exclude contamination
  
    clean_table_pbmc<-anti_join(table_pbmc,contamination %>% 
                                  dplyr::rename(Participant_PPID=Participant_PPID.x),by=c('cdr3nt_trb','v_trb','j_trb','Participant_PPID')) %>% 
      group_by(Participant_PPID,cdr3aa_trb,v_trb,j_trb,Group3,real_path) %>% 
      summarise(freq=sum(freq),
                count=sum(count)) %>% ungroup()
    
    clean_table_tissues<-anti_join(table_tissues,contamination %>% 
                                     dplyr::rename(Participant_PPID=Participant_PPID.y),by=c('cdr3nt_trb','v_trb','j_trb','Participant_PPID')) %>% 
      group_by(Participant_PPID,cdr3aa_trb,v_trb,j_trb,Group3,real_path) %>% 
      summarise(freq=sum(freq),
                count=sum(count)) %>% ungroup()
    
    shared_clones<-inner_join(clean_table_pbmc,clean_table_tissues,by=c('cdr3aa_trb','v_trb','j_trb')) %>% 
      group_by(Participant_PPID.x,Participant_PPID.y,age.x,age.y,Group.x,Group.y,real_path.x,real_path.y) %>% 
      summarise(n_shared_clones=n(),
                weighted_n_shared_clones=n_shared_clones/(as.numeric(nrow(clean_table_pbmc))*
                                                            as.numeric(nrow(clean_table_tissues)))) %>% 
      mutate(size_1st_rep_no_cont=nrow(clean_table_pbmc),size_2nd_rep_no_singles=nrow(clean_table_tissues)) %>% 
      dplyr::rename(ppid_1=Participant_PPID.x,ppid_2=Participant_PPID.y) %>% 
      mutate(age_group=paste0(age.x,'_',age.y),
             group_comparison=paste0(Group.x,'_',Group.y))
    return(shared_clones)
    }
  
  
    lapply(list_of_pairs,pair_comparison) %>% bind_rows() %>% 
    mutate(group_comparison=case_when(group_comparison=='B_A'~'A_B',
                                      group_comparison=='C_A'~'A_C',
                                      group_comparison=='C_B'~'B_C',
                                      TRUE~group_comparison)) %>% 
    mutate(OIT=case_when(group_comparison=='A_B'~'OIT_OIT',
                         group_comparison=='A_C'~'OIT_placebo',
                         group_comparison=='B_C'~'OIT_placebo',
                         group_comparison=='C_C'~'placebo_placebo',
                         TRUE~'OIT_OIT'))
}
##shared clones between pbmc and tissues----

find_n_shared_clones_between_tissue_and_pbmc<-function(data_tissues,data_pbmc,list_of_ppid){

  pbmc_tissue_comparison<-function(ppid){
  
    list_of_tissues<-data_tissues |> 
      filter(Participant_PPID==ppid) |> 
      pull(tissue) |> unique()
 print(list_of_tissues)
  lapply(list_of_tissues,function(tissue_part){
  print(tissue_part)
     
   table_pbmc<-data_pbmc %>% 
     filter(Participant_PPID==ppid)
   
   table_tissues<-data_tissues %>% 
     filter(Participant_PPID==ppid,
            tissue==tissue_part)
   
   shared_clones<-inner_join(table_pbmc,table_tissues,
                             by=c('cdr3aa_trb','v_trb','j_trb','Participant_PPID')) |> 
     group_by(tissue,Participant_PPID,week) |> 
     summarise(n_shared_clones_by_tissue=n(),
               weighted_n_shared_clones=n_shared_clones_by_tissue/(as.numeric(nrow(table_pbmc))*
                                                           as.numeric(nrow(table_tissues)))) %>% 
     mutate()
   return(shared_clones)
  }) |> bind_rows()
      }
  lapply(list_of_ppid,pbmc_tissue_comparison)
} |> bind_rows()
## baseline clones----
baseline_allPat<-process_mixcr_clones(list_samples = second_run_files_deepSeq_new %>% 
                                        filter(week=='BL01' | week=='BL02' | week=='BL03' ) %>% 
                                        pull(real_path))
baseline_allPat %<>% 
  inner_join(second_run_files_deepSeq_new %>% 
               select(real_path,Participant_PPID,age,Group) %>% unique(),by='real_path')

shared_clones_bl<-find_n_shared_clones(data=baseline_allPat,list_of_pairs = baseline_allPat %>% pull(Participant_PPID) %>% unique() %>% 
                       combn(.,2,c,simplify=F))

baseline_allPat_corr<-baseline_allPat %>% pull(Participant_PPID) %>% unique() %>% 
  combn(.,2,c,simplify=F) %>% as.data.frame() %>% 
  pivot_longer(cols=starts_with('c'),names_to = 'patient_comparison') %>% 
  mutate(ppid_1=str_extract(patient_comparison,'c..P\\d+'),
         ppid_1=str_remove(ppid_1,'c..'),
         ppid_2=str_extract(patient_comparison,'....P\\d+'),
         ppid_2=str_remove(ppid_2,'....')) %>% 
  select(-value,-patient_comparison) %>% 
  left_join(shared_clones_bl,by=c('ppid_1','ppid_2')) %>% unique() %>% 
  mutate(age.x=ifelse(is.na(age.x),'Juvenile',age.x),
         age.y=ifelse(is.na(age.y),'Juvenile',age.y),
         n_shared_clones=ifelse(is.na(n_shared_clones),0,n_shared_clones),
         weighted_n_shared_clones=ifelse(is.na(weighted_n_shared_clones),0,weighted_n_shared_clones),
         Group.x=ifelse(is.na(Group.x),'A',Group.x),
         age_group=ifelse(is.na(age_group),'Juvenile_Juvenile',age_group),
         Group.y=ifelse(is.na(Group.y),'A',Group.y)) %>% 
  mutate(Group.y=case_when(ppid_2=='P053'~'B',
                           ppid_2=='P050'~'C',
                           TRUE~Group.y),
         group_comparison=paste0(Group.x,'_',Group.y)) %>% 
  mutate(OIT=ifelse(is.na(OIT),case_when(group_comparison=='A_C'~'OIT_placebo',
                       TRUE~'OIT_OIT'),OIT))
  
baseline_allPat_corr %>% 
  ggplot(aes(x=age_group,y=weighted_n_shared_clones))+
  geom_quasirandom(alpha=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.6)+
  geom_text(data=baseline_allPat_corr %>% group_by(age_group) %>% summarise(n=n()),
            aes(x=age_group,y=0.000000015,label=n))+
 scale_y_log10(oob = scales::squish_infinite)+
  stat_compare_means()+
  xlab('Age Group')+
  ylab('Weighted n of shared clones (log10)')+
  labs(title='Number of shared clones at the baseline')+
  theme_bw()
  
#to save the df
write_tsv(baseline_allPat_corr,file = '/data/vskatova/vskatova/poised_trb/data_sharedCL_bl.tsv')
#to save the plot
ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/sharedCL_baseline_byAgeGroup.png',
       width =5,height =4 ,device = 'png',dpi = 400)

##week104 clones----
w104_allPat<-process_mixcr_clones(list_samples = deep_seq_run1_metadata %>% 
                                        filter(week=='W104' | week=='W104A' | week=='W104B') %>% 
                                        pull(real_path))

w104_allPat %<>% inner_join(deep_seq_run1_metadata %>% select(real_path,Participant_PPID,Group,age) %>% unique(),
                                                             by='real_path')

shared_clones_W104<-find_n_shared_clones(data=w104_allPat,
                                        list_of_pairs = w104_allPat %>% pull(Participant_PPID) %>% unique() %>% 
                                         combn(.,2,c,simplify=F))
shared_clones_W104 %<>% ungroup()

shared_clones_w104_test<-dunn_test(data = shared_clones_W104 %>% 
  mutate(weighted_n_shared_clones_new=weighted_n_shared_clones*1e6) ,formula = weighted_n_shared_clones ~ OIT) %>% 
  add_y_position(formula = weighted_n_shared_clones_new ~ OIT,scales = 'free_y') %>% 
  mutate(y.position = y.position /1e6)

shared_clones_W104 %>% 
  #filter(age_group=='Juvenile_Juvenile') %>% 
  ggplot(aes(x=age_group,y=weighted_n_shared_clones,color=factor(age_group)))+
  geom_quasirandom(pch = 21, alpha=0.3,dodge.width = 0.8)+
  geom_boxplot(outlier.shape = NA,alpha=0.6)+
  scale_color_manual(values=c("#009999", "#FFCC99","#5F31CC"))+
  #geom_text(data=shared_clones_W104 %>% group_by(OIT) %>% summarise(n=n()),
            #aes(x=OIT,y=0.00000001,label=n))+
  xlab('Age Group')+
  ylab('Weighted n of shared clones (log10)')+
  labs(title='Number of shared clones at week 104')+
  facet_grid(~OIT,scales = 'free_y')+
  theme_bw()

baseline_allPat_corr %>% 
  mutate(week_common='BL') %>% 
  bind_rows(shared_clones_W104 %>% mutate(week_common='W104')) %>% 
  filter(age_group=='Juvenile_Juvenile') %>% 
  ggplot(aes(x=week_common,y=weighted_n_shared_clones))+
  facet_wrap(~OIT,scales = 'free_y',nrow =  1)+
  geom_quasirandom()+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means()

##pbmc and tissues----
###week 0----
w0_not_downsampled<-process_mixcr_clones(list_samples = metadata_3d_w0_w52 %>% filter(week=='W000') %>% 
                                           pull(real_path))

w0_not_downsampled<-grouping_mixcr_clones(file_mixcr_clones = w0_not_downsampled )

w0_not_downsampled<-w0_not_downsampled |> 
  inner_join(metadata_3d_w0_w52 |> 
               select(Participant_PPID,real_path,Group3) |> unique(),by='real_path')

week0_tissues<-process_mixcr_clones(list_samples = final_metadata %>% filter(tissue!='PBMC') %>% #upload tissue repertoires week 0 
                                      filter(week=='W000') %>% 
                                      mutate(sample=str_replace(sample,'.txt','clones_trb_prod.tsv')) |> 
                                      mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr_V4/',sample)) |> 
                                      pull(real_path))
week0_tissues<-grouping_mixcr_clones(file_mixcr_clones = week0_tissues )

week0_tissues<- week0_tissues |>  inner_join(final_metadata %>% 
                               mutate(sample=str_replace(sample,'.txt','clones_trb_prod.tsv')) |> 
                               mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr_V4/',sample)) |> 
                               select(real_path,Participant_PPID,week,tissue,Group3) %>% unique(),
                               by='real_path')

w0_pbmc_tissues_shared<-find_n_shared_clones_between_tissue_and_pbmc(data_tissues =week0_tissues,data_pbmc = w0_not_downsampled,
                                             list_of_ppid = w0_not_downsampled |>
                                               filter(Participant_PPID %in% week0_tissues$Participant_PPID) |>  
                                               pull(Participant_PPID) |>  unique())
###baseline----

bl_w0_pbmc_tissue_shared<-find_n_shared_clones_between_tissue_and_pbmc(data_pbmc=bl01_repertoire_not_downsampled |> 
  left_join(second_run_files_deepSeq_new |> 
              select(Participant_PPID,real_path,Group3),by='real_path'),data_tissues =week0_tissues,
              list_of_ppid = bl01_repertoire_not_downsampled |>
    left_join(second_run_files_deepSeq_new |> 
                select(Participant_PPID,real_path,Group3),by='real_path') |> 
    filter(Participant_PPID %in% week0_tissues$Participant_PPID) |> pull(Participant_PPID) |> unique()) 
  
bl_w0_pbmc_tissue_mean<-bl_w0_pbmc_tissue_shared |> 
  dplyr::rename(week_com=week) |> 
  mutate(week='BL') |> 
  bind_rows(w0_pbmc_tissues_shared |> mutate(week_com='W000')) |> 
  group_by(Participant_PPID,tissue,week_com) |> 
  summarise(weighted_n_shared_clones_mean=mean(weighted_n_shared_clones)) 

###week 104----
week104_tissues<-process_mixcr_clones(list_samples = final_metadata %>% filter(tissue!='PBMC') %>% #upload tissue repertoires week 0 
                                        filter(week=='W104') %>% 
                                        mutate(sample=str_replace(sample,'.txt','clones_trb_prod.tsv')) |> 
                                        mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr_V4/',sample)) |> 
                                        pull(real_path))

week104_tissues<-grouping_mixcr_clones(week104_tissues)
week104_tissues %<>%  inner_join(final_metadata %>% 
                                   mutate(sample=str_replace(sample,'.txt','clones_trb_prod.tsv')) |> 
                                   mutate(real_path=paste0('/data/vskatova/vskatova/poised_trb/poised_clonsets/mixcr_V4/',sample)) |> 
                                   select(real_path,Participant_PPID,week,tissue,Group3) %>% unique(),
                                 by='real_path')

w104_pbmc_tissues_shared<-find_n_shared_clones_between_tissue_and_pbmc(data_tissues = week104_tissues,
                                                                       data_pbmc = w104_bulk_corrected_freq_not_downsampled |> 
                                                                         select(-week) |> 
                                                                         separate(vj_trb,into=c('v_trb','j_trb'),sep = '_'),
                                                                     list_of_ppid = week104_tissues |>
                                                                       filter(Participant_PPID %in% week104_tissues$Participant_PPID) |>  
                                                                       pull(Participant_PPID) |>  unique())
#w0-w104 comparison
bl_w0_pbmc_tissue_mean |> 
  dplyr::rename(weighted_n_shared_clones=weighted_n_shared_clones_mean,
                week=week_com) |> 
  bind_rows(w104_pbmc_tissues_shared) |>
  filter(tissue!='Duodenum_Stomach') |> 
  inner_join(final_metadata |> select(Participant_PPID,Group2) |> unique(),
             by='Participant_PPID') |> 
  group_by(tissue,Group2,Participant_PPID) |> 
  mutate(n_weeks=n_distinct(week)) |> 
  filter(n_weeks==2) |> 
         #Group2=='AB') |> View()
  ggplot(aes(x=week,y=weighted_n_shared_clones))+
  geom_boxplot(aes(color=tissue))+
  geom_point(aes(color=tissue))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_line(aes(group=Participant_PPID,color=tissue))+
  stat_compare_means(paired = TRUE)+
  facet_wrap(~tissue+Group2,scales = 'free')+
  labs( y='Weighted n shared clones with PBMC (by ppid) ')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/N_shared_clones_tissue_pbmc_w0_w104.png',
       width =10,height =7 ,device = 'png',dpi = 400)

#tissues comparison between each other
bl_w0_pbmc_tissue_mean |> 
  dplyr::rename(weighted_n_shared_clones=weighted_n_shared_clones_mean,
                week=week_com) |> 
  bind_rows(w104_pbmc_tissues_shared) |>
  filter(tissue!='Duodenum_Stomach') |> 
  inner_join(final_metadata |> select(Participant_PPID,Group2) |> unique(),
             by='Participant_PPID') |> 
  mutate(tissue=factor(x=tissue,
                       levels=c('prox_esophagus','mid_esophagus','dist_esophagus','stomach','duodenum'))) %>%
  mutate(OIT=ifelse(Group2=='AB','OIT','Placebo')) |> 
  ggplot(aes(x=tissue,y=weighted_n_shared_clones,color=tissue))+
  geom_boxplot(alpha=0.5)+
  geom_point(alpha=0.5)+
  stat_compare_means()+
  facet_wrap(~week+OIT,scales = 'free')+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs( y='Weighted n shared clones with PBMC (by ppid) ')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/N_shared_clones_tissue_pbmc.png',
       width =12,height =7 ,device = 'png',dpi = 400)
 
#expanded clones(sc count > 1)----
сd4_seq_corrected_only_trb_expanded<- сd4_seq_corrected_only_trb |> 
  filter(Count_single_cell>1)

##bl in bl, w104 in w104, w117 in w117----
###large subsets----
#bl
pr_exp_bl_in_bulk_bl_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= bl01_repertoire_not_downsampled |> 
                                                                                 inner_join(second_run_files_deepSeq_new %>% 
                                                                                              select(real_path,Participant_PPID,week) %>% unique(),by='real_path'),
                                                                               data_sc = сd4_seq_corrected_only_trb_expanded, 
                                                                          week_sc = 'BL',week_bulk = 'BL01')
#week 104
pr_exp_w104_in_bulk_w104_large_subset<-find_sc_clones_in_bulk_by_large_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled,
                                                                                   data_sc = сd4_seq_corrected_only_trb_expanded, 
                                                                                  week_sc = 'W104',  week_bulk = 'W104')



#bl-w104
bl_w104_large_subsets_exp_clones<-pr_exp_bl_in_bulk_bl_large_subset |> 
  bind_rows(pr_exp_w104_in_bulk_w104_large_subset) |> 
  group_by(Participant_PPID,Cell_type) |> 
  mutate(n_tp=n_distinct(week_common)) |> 
  filter(n_tp==2) |>
  inner_join(cd4_seq_corrected_only_trb |> ungroup() |> 
               select(Participant_PPID,Group_dose_outcome) |> unique(),by='Participant_PPID') |> 
  filter(Cell_type!='PR_CD154+Teff_Na') |> 
  inner_join(cd4_seq_corrected |> 
               select(OIT,Participant_PPID) |> unique(),by='Participant_PPID')

#statistics
stat_bl_w104_exp<-bl_w104_large_subsets_exp_clones %>% 
  group_by(Cell_type,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, 
              paired = T,comparisons =list(c('BL01','W104'))) %>% 
  add_significance(p.col = 'p') |> 
  add_y_position() |> 
  select(-y.position) |> 
  unique()

bl_w104_large_subsets_exp_clones |> 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset))+
  geom_boxplot(outlier.shape = NA,aes(color=OIT))+
  geom_line(aes(group=Participant_PPID),alpha=0.2)+
  geom_point(alpha=0.7,size=2,aes(color=OIT))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  facet_wrap(~OIT+Cell_type,scales = 'free_y',nrow = 2 )+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))+
  labs(y='Summed freq (by Patient Id)')+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/expanded_clones_large_subs_bl_w104_OIT_pl.png',
       width = 6,height = 5,device = 'png',dpi = 400)

###Th1/Th2a/Th2conv----

pr_exp_bL_bulk_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= bl01_repertoire_not_downsampled %>% 
                                                            inner_join(second_run_files_deepSeq_new %>% 
                                                                         select(real_path,Participant_PPID,week) %>% unique(),by='real_path'), 
                                                          data_sc = сd4_seq_corrected_only_trb_expanded, 
                                                          week_sc = 'BL',week_bulk = 'BL01')

#w104  
pr_exp_w104_bulk_by_subset<-find_sc_clones_in_bulk_by_subset(data_bulk= w104_bulk_corrected_freq_not_downsampled, 
                                                                       data_sc = сd4_seq_corrected_only_trb_expanded,
                                                                       week_sc = 'W104',week_bulk = 'W104')

#w0,w52 and w104 together
exp_clones_w0_w104_bulk<-pr_exp_bL_bulk_by_subset %>% 
  bind_rows(pr_exp_w104_bulk_by_subset) %>%
  left_join(cd4_seq_corrected_only_trb %>% ungroup() %>% select(Participant_PPID,OIT) %>% unique(),by='Participant_PPID') %>% 
  group_by(Participant_PPID,OIT,cell_subset) %>% 
  mutate(n_tp=n_distinct(week_common)) %>% 
  filter(n_tp==2) |> 
  filter(cell_subset=='PR_Teff_Me_Th2a'|cell_subset=='PR_Th2conv'|cell_subset=='PR_Teff_Me_Th1_CTL')

#statistics
stat_exp_pr_clones_w0_w104_bulk<-exp_clones_w0_w104_bulk %>% 
  group_by(cell_subset,OIT) %>% 
  wilcox_test(sum_freq_in_bulk_by_subset ~ week_common, paired = T,
              comparisons =list(c('BL01','W104'))) %>% 
  add_significance(p.col = 'p') %>% 
  add_y_position()

exp_clones_w0_w104_bulk %>% 
  ggplot(aes(x=week_common,y=sum_freq_in_bulk_by_subset,color=OIT))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=Participant_PPID))+
  geom_point(alpha=0.7,size=2)+
  facet_wrap(~cell_subset+OIT,scales = 'free_y',nrow =  3 )+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  stat_pvalue_manual(data=stat_exp_pr_clones_w0_w104_bulk,
                     label = 'p.signif',hide.ns = TRUE)+
  labs( y='Summed frequency (by patient ID)')+
  theme_classic()+
  theme(strip.background = element_blank())+
  xlab(element_blank())+
  theme(legend.position="none")+ 
  theme(strip.text.x = element_text(size = 12))

ggsave(filename = '/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/expanded_w0_w52_w104_bulk.png',
       width =5,height = 6,device = 'png',dpi = 400)

#edgeR processing----

patients_more_than_3tp<-metadata_3d_w0_w52 |> 
  bind_rows(first_run_deep_seq_metadata,
            second_run_files_deepSeq_new) |> 
  select(Participant_PPID,week,real_path) |> 
  unique() |>
  group_by(Participant_PPID) |> 
  summarise(n_tp=n_distinct(week),
            weeks=paste(week,collapse = ',')) |>
  filter(n_tp>4)

metadata_edgeR<-metadata_3d_w0_w52 |> 
  bind_rows(first_run_deep_seq_metadata,
            second_run_files_deepSeq_new) |> 
  select(Participant_PPID,week,real_path) |> 
  unique() |>
  filter(Participant_PPID %in% patients_more_than_3tp$Participant_PPID) |> 
  mutate(week=case_when(week=='BL01'~'W000',
                        week=='BL02'~'W000',
                        week=='BL03'~'W000',
                        week=='BL04'~'W000',
                        week=='W104A'~'W104',
                        week=='W104B'~'W104',
                        week=='W117A'~'W117',
                        week=='W117B'~'W117',
                        week=='W36'~'W052',
                        week=='W34'~'W052',
                        week=='W065'~'W052',
                        week=='W24'~'W052', TRUE~week))
  
edgeR_clone_search<-function(metadata){

  pat_list<-metadata %>% 
     pull(Participant_PPID) %>% unique() #лист пациентов
    
cloneCounts_wide<-lapply(pat_list,function(pat){
  
  metadata_filtered_patID<-metadata %>%
    filter(Participant_PPID==pat)
 
   weeks_list<-metadata_filtered_patID |> pull(week) # лист недель
 print(weeks_list)
  paths_list<-metadata_filtered_patID %>% pull(real_path) #лист путей файлов
 
   mixcr_clonsets<-process_mixcr_clones(paths_list) #общая таблица с клонами
 
   grouped_mixcr_table<-grouping_mixcr_clones(mixcr_clonsets) # группированные каунты клонов
  
   table_counts_cdr3aa<- grouped_mixcr_table |> 
    mutate(aaSeqCDR3_V_J=paste(cdr3aa_trb,
                               v_trb,
                               j_trb, sep="_")) %>% 
    select(real_path,aaSeqCDR3_V_J,count_bulk) %>% 
    pivot_wider(id_cols = aaSeqCDR3_V_J,
                values_from = count_bulk,
                names_from = real_path,
                values_fill = 0) %>% 
    filter(rowMeans(select(.,-1))>1)
  
  table_edgeR <- table_counts_cdr3aa[,-1] %>% as.data.frame()
  rownames(table_edgeR) <- table_counts_cdr3aa %>% pull(aaSeqCDR3_V_J)
  print(weeks_list)
  y <- DGEList(counts = table_edgeR,group = weeks_list) 
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  
time_pairs<-combn(weeks_list,m=2,FUN = c,simplify = F)
print(time_pairs)
  lapply(time_pairs,function(time_pair){
    print(time_pair)
    test_vector<-c(time_pair |> unlist())
    de <- exactTest(y,pair = test_vector,dispersion = "trended")
    result <- topTags(de,n=1000,p.value = 1, sort.by="PValue")$table
    result<- as.data.table(result,keep.rownames = TRUE)
    result_filtered<-result %>% filter(PValue<0.05)
    print(nrow(result_filtered))
    path<-paste0('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/edgeR/',
                 pat,'_',time_pair[1],'_',time_pair[2],'.txt')
    print(path)
    readr::write_tsv(x=result_filtered,
                     file=path,
                     quote = 'none')
    
  }) 
})  
}

edgeR_results<-edgeR_clone_search(metadata =metadata_edgeR ) #run the function

edgeR_clones<-read_tsv(dir('/data/vskatova/vskatova/poised_trb/poised_trb_deepSeq/analysis/edgeR/',
                               full.names = TRUE),id = 'real_path') |> 
  mutate(Participant_PPID=str_extract(real_path,'P\\d+'),
         week1=str_extract(real_path,'W\\d+'),
         week2=str_extract(real_path,'W\\d+.txt'),
         week2=str_remove(week2,'.txt')) 

metadata_edgeR2<-metadata_edgeR |> 
 full_join(metadata_edgeR,by='Participant_PPID') |>  
  mutate(week1_week2=paste0(week.x,'_',week.y)) |> 
  select(-week.x,-week.y,-real_path.x,-real_path.y) |> 
  filter(week1_week2=='W000_W104'|
         week1_week2=='W000_W117'| 
         week1_week2=='W000_W052'|
         week1_week2=='W052_W104'| 
         week1_week2=='W052_W117') |> unique()

nClones_ncReads_new<-nClones_ncReads_realpath |>  
  mutate(week=case_when(week=='BL01'~'W000',
                        week=='BL02'~'W000',
                        week=='BL03'~'W000',
                        week=='BL04'~'W000',
                        week=='W104A'~'W104',
                        week=='W104B'~'W104',
                        week=='W117A'~'W117',
                        week=='W117B'~'W117',
                        week=='W36'~'W052',
                        week=='W34'~'W052',
                        week=='W065'~'W052',
                        week=='W24'~'W052', TRUE~week)) |> 
  group_by(Participant_PPID,week) |> 
  summarise(nClones_mean=mean(nclones),
            nReads_mean=mean(nReads))

edgeR_clones_processed<-edgeR_clones |> 
  mutate(week1_week2=paste0(week1,'_',week2)) |>
  filter(week1_week2!='W104_W117') |> 
  mutate(week1_week2=case_when(week1_week2=='W104_W000'~'W000_W104',
                               week1_week2=='W117_W000'~'W000_W117',
                               week1_week2=='W052_W000'~'W000_W052',
                               week1_week2=='W104_W052'~'W052_W104',
                               week1_week2=='W117_W052'~'W052_W117',
                               TRUE~week1_week2
                               )) |> 
  select(-logFC,-logCPM,-real_path,-week1,-week2) |> unique() 

##difference in number of edgeR clones between groups or week pairs----
edgeR_normalized<-edgeR_clones_processed |> 
  #filter(FDR<1) |> 
  group_by(week1_week2,Participant_PPID) |> 
  summarise(n_clones_by_week_pair=n()) |> 
  right_join(metadata_edgeR2,by=c('Participant_PPID','week1_week2')) |> 
  inner_join(final_metadata |> select(Group,Participant_PPID) |> unique(),by='Participant_PPID') |> 
  separate(col = week1_week2,into = c('week1','week2'),sep = '_') |> 
  left_join(nClones_ncReads_new |> select(nClones_mean,Participant_PPID,week) |> 
              dplyr::rename(week1=week),by=c('Participant_PPID','week1')) |> 
  left_join(nClones_ncReads_new |> select(nClones_mean,Participant_PPID,week) |> 
              dplyr::rename(week2=week),by=c('Participant_PPID','week2')) |> 
  unique() |> 
  dplyr::rename(nclones_1=nClones_mean.x,
                nclones_2=nClones_mean.y) |> ungroup() |> 
  mutate(nclones_1=as.numeric(nclones_1),
         nclones_2=as.numeric(nclones_2)) |> 
  mutate(n_clones_by_week_pair=ifelse(is.na(n_clones_by_week_pair),0,n_clones_by_week_pair)) |> 
  mutate(normalized_n_clones_by_week=(n_clones_by_week_pair/(nclones_1*nclones_2))) |> 
  mutate(week1_week2=paste0(week1,'_',week2))

edgeR_normalized |> 
 # mutate(normalized_n_clones_by_week=ifelse(normalized_n_clones_by_week== 0,
                                         #   min(normalized_n_clones_by_week[normalized_n_clones_by_week>0])/2,
                                         #   normalized_n_clones_by_week)) %>% 
  filter(normalized_n_clones_by_week>0) |> 
  ggplot(aes(x=Group,y=normalized_n_clones_by_week))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  stat_compare_means()+
  facet_wrap(~week1_week2,scales = 'free')+
  scale_y_log10()
  #ylim(0,1e-11)

##intersection with sc clones----

edgeR_clones_in_SC_data<-edgeR_clones_processed |> 
  inner_join(cd4_seq_corrected_only_trb |> mutate(rn=paste0(cdr3aa_trb,'_',vj_trb)),
             by=c('Participant_PPID','rn'))
  
##motif search----

edgeR_clones_processed2<-edgeR_clones_processed |> 
  dplyr::rename(cdr3aa_vj_trb=rn)

matrixClones<-stringdistmatrix(edgeR_clones_processed2 %>% pull(cdr3aa_vj_trb),
                               edgeR_clones_processed2 %>% pull(cdr3aa_vj_trb),method = 'hamming') %>% as.matrix() %>% 
  set_rownames( edgeR_clones_processed2 %>% unite('x',cdr3aa_vj_trb,Participant_PPID,week1_week2) |>  pull(x))  %>%
  set_colnames(edgeR_clones_processed2  %>% unite('y',cdr3aa_vj_trb,Participant_PPID,week1_week2) |>  pull(y)) 

New_table<-matrixClones %>%  
  reshape2::melt() %>%
  set_colnames(c("seq1","seq2","dist")) %>%
  filter(dist<2)

New_table_new<- New_table %>% separate(col=seq1,into=c('cdr3aa','V_trb','J_trb','patients'),
                                       sep='_',remove=TRUE) %>% separate(col=seq2,into=c('CDR3aa_pog','V_pog','J_pog'),
                                                                         sep='_',remove=TRUE) %>% 
  filter(V_alice==V_pog & J_alice ==J_pog)

Innerclones<-test_new_2 %>% mutate(CDR3aaPog=CDR3.amino.acid.sequence %in% New_table_new$CDR3aa_alice) %>% 
  filter(CDR3aaPog==TRUE)