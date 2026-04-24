
phyT.noContam <- transform_sample_counts(phyR.noContam, function(x) x/sum(x))
otu_table(phyT.noContam)[is.na(otu_table(phyT.noContam))] <- 0
sample_data(phyT.noContam)$number <- factor(sample_data(phyT.noContam)$number, levels=c("1","4","7","10","14"))
sample_data(phyT.noContam)$sampleSite <- factor(sample_data(phyT.noContam)$sampleSite, levels=c("Jejunum","Ileum","Caecum","Colon","Fecal"))
sample_data(phyT.noContam)$mutant <- with(sample_data(phyT.noContam), ifelse(genotype%in%c("WT1","WT2"), "WT", "dtssBC")) 


otu.noContam <- data.frame(otu_table(phyT.noContam))
sd.noContam <- data.frame(sample_data(phyT.noContam))

otu.noContam2 <- rbind(otu.noContam, data.frame(tax_table(phyT.noContam))[,c("Spp")])
dir.create("output") # create directory for output files to go
write.csv(otu.noContam2,"output/table_abundance.csv") # Save otu table as .csv
write.csv(sd.noContam, "output/table_metadata.csv") # Save otu table as .csv

ps.noContam <- psmelt(phyT.noContam)
ps.noContam$Spp <- with(ps.noContam, paste(Genus, Species))
data.plot <- subset(ps.noContam, Control=="sample" & !sampleID%in%c("8317_087","8317_076"))
data.plotBact <- subset(data.plot, Genus=="Bacteroides")
data.plotBact$number <- factor(data.plotBact$number, levels=c("1","4","7","10","14"))
data.plotBact$sampleSite <- factor(data.plotBact$sampleSite, levels=c("Jejunum","Ileum","Caecum","Colon","Fecal"))

pBar.a <- ggplot(data=subset(data.plotBact, community=="Total" & number=="14"), aes(x=sampleSite, y=Abundance, fill=Spp)) +
  geom_bar(position = "stack", stat="identity") +  #geom_jitter(width=0.2, height=0.01) +
  facet_grid(genotype+number~community+sampleRvalue, scales="free_x", space="free_x") +
  theme_bw() +
  #xlab("Culture Result") +
  #ylab("Percent Bacteria") +
  scale_fill_manual(values=pal.CB) +
  theme(text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),axis.title.x =element_blank(),panel.spacing = unit(0, "lines"))
pBar.a
ggsave(pBar.a, filename="plots/pBar_total14.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file

pBar.b <- ggplot(data=subset(data.plotBact, community=="Transconjugants" & sampleSite=="Fecal"), aes(x=number, y=Abundance, fill=Spp)) +
  geom_bar(position = "stack", stat="identity") +  #geom_jitter(width=0.2, height=0.01) +
  facet_grid(genotype~community+sampleRvalue, scales="free_x", space="free_x") +
  theme_bw() +
  #xlab("Culture Result") +
  #ylab("Percent Bacteria") +
  scale_fill_manual(values=pal.CB) +
  theme(text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),axis.title.x =element_blank(),panel.spacing = unit(0, "lines"))
pBar.b
ggsave(pBar.b, filename="plots/pBar_totalFecal.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file


# alpha diversity
data.alpha <- unique(data.plotBact[,c("sampleID","genotype","sampleSite","sampleRvalue","number","community","Observed","Shannon","Simpson")])

pAlpha.14 <- ggplot(data=subset(data.alpha, community!="unk" & number=="14"), aes(x=sampleSite, y=Simpson, color=community)) +
  geom_jitter(size=3) +
  facet_grid(.~mutant, scales="free_x", space="free_x") +
  theme_bw()
pAlpha.14
ggsave(pAlpha.14, filename="plots/pSimpson_sampleSite.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file

pAlpha.fecal <- ggplot(data=subset(data.alpha, community!="unk" & sampleSite=="Fecal"), aes(x=number, y=Simpson, color=community)) +
  geom_jitter(size=3) +
  facet_grid(.~mutant, scales="free_x", space="free_x") +
  theme_bw()
pAlpha.fecal
ggsave(pAlpha.fecal, filename="plots/pSimpson_fecal.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file


# beta diversity?

phy.ord <- subset_taxa(subset_samples(phyT.noContam, genotype!="unk" & !is.na(sampleRvalue)), Genus=="Bacteroides")
phy.ord <- subset_taxa(phy.ord, taxa_sums(phy.ord)>0)
phy.ord <- subset_samples(phy.ord, sample_sums(phy.ord)>0)

ord.a = ordinate(subset_samples(phy.ord,  sampleSite=="Fecal"), "NMDS", "bray") 

ord.a <- plot_ordination(subset_samples(phy.ord, sampleSite=="Fecal"), ord.a, color = "number", shape="sampleRvalue") + 
  geom_point(size = 3) +
  facet_grid(community~mutant) +
  theme_bw() +
  labs(title="Fecal") +
  theme(legend.position="bottom", text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),axis.title.x =element_blank(),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=c("#E92028","#F46F22","#FFE000","#90C73B","#0988BC","#5D2E91"))
ord.a
ggsave(ord.a, filename="plots/nmdsBray_fecal.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file

ord.b = ordinate(subset_samples(phy.ord,  number=="14"), "NMDS", "bray") 

ord.b <- plot_ordination(subset_samples(phy.ord, number=="14"), ord.b, color = "sampleSite", shape="sampleRvalue") + 
  geom_point(size = 3) +
  facet_grid(community~mutant) +
  theme_bw() +
  labs(title="Fecal") +
  theme(legend.position="bottom", text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),axis.title.x =element_blank(),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=c("#E92028","#F46F22","#FFE000","#90C73B","#0988BC","#5D2E91"))
ord.b
ggsave(ord.b, filename="plots/nmdsBray_14.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file

