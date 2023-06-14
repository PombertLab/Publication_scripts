#!/usr/bin/env Rscript

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17853474.fastq.gz.50604.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17853474.fastq.gz.50604.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17853474.fastq.gz.50604.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17853474.fastq.gz.50604.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17853474.fastq.gz.50604.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17853474.fastq.gz.50604.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17858634.fastq.gz.50602.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17858634.fastq.gz.50602.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17858634.fastq.gz.50602.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17858634.fastq.gz.50602.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17858634.fastq.gz.50602.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17858634.fastq.gz.50602.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17865590.fastq.gz.50506.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17865590.fastq.gz.50506.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17865590.fastq.gz.50506.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17865590.fastq.gz.50506.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR17865590.fastq.gz.50506.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR17865590.fastq.gz.50506.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR23560420.fastq.gz.50451.fasta.minimap2.both.cc.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR23560420.fastq.gz.50451.fasta.minimap2.both.cc.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR23560420.fastq.gz.50451.fasta.minimap2.both.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR23560420.fastq.gz.50451.fasta.minimap2.both.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

matrix <- as.matrix(read.table("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/SRR23560420.fastq.gz.50451.fasta.minimap2.both.sub.sorted.10.90.tsv", header=FALSE, sep="\t"))
pdf("/home/jpombert/Analyses/Encephalitozoon/GenDiv/for_R_plots/PDFs/SRR23560420.fastq.gz.50451.fasta.minimap2.both.sub.sorted.10.90.pdf")
library(MASS)
colnames(matrix)[1] <- "allelic frequency"
fit <- fitdistr(matrix, "normal")
class(fit)
para <- fit$estimate
hist(matrix, breaks = 50, prob = TRUE, xlim=c(0,100), ylim=c(0,0.12), xlab = "allelic frequencies")
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
while (!is.null(dev.list()))  dev.off()

