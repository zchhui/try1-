
rm(list = ls()) 
options(stringsAsFactors = F)
f='GSE100795_eSet.Rdata'

library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE100795', destdir=".",
                 AnnotGPL = T, ##ע���ļ�   
                 getGPL = T)   ##ƽ̨�ļ�   
  save(gset,file=f)            ##���浽����
}#����'GSE100795'�������
load('GSE100795_eSet.Rdata')   
class(gset)
gset[[1]]            ##������ȡgset
class(gset[[1]])
dat1 <- exprs(gset[[1]])##ʹ�ú���exprs��ȡ�����������
pd1 <- pData(gset[[1]]) ##ʹ�ú���pData��ȡ����������Ϣ
fd1 <- fData(gset[[1]])##ƽ̨ע����Ϣ

write.table(pd1,"D:\\Rworkplace\\GSE100795\\GSE100795_gsegroup.txt",sep="\t",quote=F)
write.table(dat1,"D:\\Rworkplace\\GSE100795\\GSE100795_raw.txt",sep="\t",quote=F)


##���ػ���ID
library(GEOquery)
GPL23126 <- getGEO('GSE100795', destdir=".")
GPL23126_anno <- table(GPL23126)
geneid <- fread("D:\\Rworkplace\\GSE100795\\GSE100795_family.txt",skip ='ID')
library(dplyr)
library(tidyr)
probe2symbol_df <- geneid %>%
  select(ID,gene_assignment) %>%
  filter(gene_assignment != "---") %>%
  separate(gene_assignment,c("drop","symbol"),sep="//") %>%
  select(-drop)
probe2symbol_df <- na.omit(probe2symbol_df)#ɾ��NAֵ

write.table(probe2symbol_df,"D:\\Rworkplace\\GSE100795\\GSE100795_geneid.txt",sep="\t",quote=F)
