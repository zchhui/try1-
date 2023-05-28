
rm(list = ls()) 
options(stringsAsFactors = F)
data <- read.table("D:\\Rworkplace2\\GSE100795\\GSE100795.txt",header=T,stringsAsFactors=F)
probe_geneid <- read.table("D:\\Rworkplace2\\GSE100795\\GSE100795_geneid.txt",header=T,stringsAsFactors=F)
names(probe_geneid) <- c("Id","symbol")

#ƥ�����id�ͻ�������
data <- as.data.frame(data)#��data������ݿ�
data$Id<- rownames(data)#��Id���ӵ�data���������µ�һ��
data <- merge(data,probe_geneid,by='Id')#merge����������data��Id��probe-geneid��ƥ��

#���̽����һ�����򣬺ϲ���һ��ȡ��ƽ��ֵ
library(limma)
data <- avereps(data[,-c(1,34)],ID=data$symbol)#ɾ����һ�к͵�53�У�����������Ϊ��������Ȼ��ȡƽ��ֵ
write.table(data,"D:\\Rworkplace\\GSE100795\\GSE\\experset.txt");



#####################
#####################Limma�Ҳ������
rm(list = ls()) 
options(stringsAsFactors = F)
dats <- read.table("D:\\Rworkplace2\\GSE100795\\GSE\\experset.txt")
dat = dats[apply(dats,1,function(x) sum(x > 0.5) >= 0.9*(ncol(dats)) ), ]


#����
group_list <- factor(c(rep('leanP',2),rep('none',12),rep('leanA',2),rep('obeseP',2),rep('none',12),rep('obeseA',2)))
table(group_list)


library(limma)
design <- model.matrix(~0+group_list)#����һ������ľ���
colnames(design)=levels(group_list)
rownames(design)=colnames(dat)
head(design)

#��������ȽϾ���
contrast.matrix1<-makeContrasts("obeseP-leanP",levels=design)

##step1��lmFitΪÿһ���������һϵ�е��������������ģ��
fit <- lmFit(dat,design)
##step2��eBayes������һ��΢��������ģ����ϣ�ͨ����Ҷ˹������һ����ͬ��ֵ�������������tͳ�������������fͳ������΢�ֱ���ʽ�Ķ������ʡ�
fit2 <- contrasts.fit(fit, contrast.matrix1)
fit2 <- eBayes(fit2)  
##step3��toptable������ģ���������ȡ��������ǰ�Ļ����
options(digits = 4)#����ȫ�ֵ�������Чλ��Ϊ4
tempOutput1 <- topTable(fit2, coef=1, n=Inf)
tempOutput1 <- na.omit(tempOutput1) #�Ƴ�NAֵ
head(tempOutput1)
a <- tempOutput1[which(tempOutput1$P.Value<0.05&abs(tempOutput1$logFC+10^(-10)) >= log2(1.5)),]
a1 <- tempOutput1[which(tempOutput1$P.Value<0.05),]


contrast.matrix2 <- makeContrasts("obeseA-leanA",levels=design)#��leanP-obeseP,leanA-obeseA����в�������Ƚ�

fit <- lmFit(dat,design)
fit2 <- contrasts.fit(fit, contrast.matrix2)
fit2 <- eBayes(fit2)  
options(digits = 4)
tempOutput2 <- topTable(fit2, coef=1, n=Inf)
tempOutput2 <- na.omit(tempOutput2) 
head(tempOutput2)

b <- tempOutput2[which(tempOutput2$P.Value<0.05&abs(tempOutput2$logFC+10^(-10)) >= log2(1.5)),]
b1 <- tempOutput2[which(tempOutput2$P.Value<0.05),]


write.table(a,"D:\\Rworkplace2\\GSE100795\\GSE\\res\\p0.05_fc1.5\\limma\\pro_obese&lean_p0.05_fc1.5.txt",sep="\t",quote=F,);
write.table(b,"D:\\Rworkplace2\\GSE100795\\GSE\\res\\p0.05_fc1.5\\limma\\adi_obese&lean_p0.05_fc1.5.txt",sep="\t",quote=F,);
write.table(a1,"D:\\Rworkplace2\\GSE100795\\GSE\\res\\p0.05\\limma\\pro_obese&lean_p0.05.txt",sep="\t",quote=F,);
write.table(b1,"D:\\Rworkplace2\\GSE100795\\GSE\\res\\p0.05\\limma\\adi_obese&lean_p0.05.txt",sep="\t",quote=F,);

write.table(dat,"D:\\Rworkplace2\\GSE100795\\GSE\\scale_experset.txt",sep="\t",quote=F,);
